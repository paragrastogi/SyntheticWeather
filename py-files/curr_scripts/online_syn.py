"""
Create Synthetic Weather based on some recorded data.
The algorithm works by creating synthetic time series over
short periods based on short histories. These short series
may be comined to obtain a longer series.
Script originally written by Parag Rastogi. Started: July 2017

Description of algorithm:
    1. Load data.
    2. Scale data using a standard scaler (subtract mean and divide by std).
    3. Enter model-fitting loop:
        a. Select 14 days of history (less when beginning).
        b. Use history to train model for next day.
        c. Sample from next day to obtain synthetic "predictions".
        d. Once the "predictions" are obtained for every day of the year,
           we are left with synthetic time series.
    4. Un-scale the data using the same scaler as (2) above.
    5. Clean / post-process the data if needed, e.g., oscillation of solar
       values around sunrise and sunset.

Features to be implemented:
    1. Load climate change model information and add that to the model.
       The outputs of climate models are usually in the form of mean changes.
       I've added them to my models in the past - need to figure out best way
       to do this now.
    2. Ability to get recorded or typical data automatically from some web
       service if the user provides coordinates or WMO station number.
    3. Ability to download climate data the same way.
    4. Ability to learn a transfer function between urban and rural stations
       for quick-and-not-too-dirty estimates of urban heat island effects.
       So long as time series of a year-ish are available from the locations
       of interest, the transfer _should_ be able to proceed without requiring
       too much information about the local urban canopy.
"""

import os.path
import random as random

import matplotlib.pyplot as plt
# import matplotlib.mlab as mlab
import numpy as np
import pandas as pd

# import scipy as sp
# from scipy import signal

import default_colours as colours

# To time the script.
import time

# Lets you copy objects if you need two
# distinct objects instead of the usual pointers.
import copy

# from sklearn.neural_network import MLPRegressor
# from sklearn.gaussian_process import GaussianProcessRegressor

# from sklearn.model_selection import GridSearchCV
# from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import KFold

# from skopt import gp_minimize
# from sklearn.datasets import fetch_mldata
from sklearn.preprocessing import StandardScaler

# import statsmodels as sm
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.graphics.tsaplots import plot_pacf
# Not needed if you don't want to plot autorrrelation plots.
# Could also use pandas:
# from pandas.tools.plotting import autocorrelation_plot


# These custom functions load and clean recorded data.
# For now, we are only concerned with ncdc and nsrdb.
import get_weather_data as gw
#from get_weather_data import load_actual
#from get_weather_data import load_typical

# Custom functions to calculate error metrics.
# import losses.
from losses import rmseloss
from losses import maeloss

# import les petites fonctionnes utiles.
from petites import yfinder
from petites import gp_funcs

# For storing data by pickling.
import pickle

# importlib is only needed if libraries need to be re-imported mid-script,
# which would usually happen only during debugging.
# import importlib as im
# from importlib import reload

__author__ = 'Parag Rastogi'

get_ipython().magic('matplotlib inline')

# Where do you want to save figures?
figpath = 'figures'
# ... and the csv outputs?
path_fldr_csv = 'pred_means'
# Check if these folders exist.
# If they don't, create them.
if not os.path.isdir(path_fldr_csv):
    os.mkdir(path_fldr_csv)
if not os.path.isdir(figpath):
    os.mkdir(figpath)
# %%

# Make a standard scaler to scale variables to \mu = 0 and \sigma = 1.
scaler = StandardScaler()

# Seed random number generators.
seed = 8760
np.random.seed(seed)
random.seed = seed

# These vectors should also be updated in the script called petites.py .
# This is the master tuple of column names, which should not be modified.
column_names = ('day', 'hour', 'tdb', 'tdp', 'rh', 'ghi', 'dni', 'dhi',
                'wspd', 'wdr')
# This is the master tuple of time variable names, which should also not
# be modified.
date_cols = ('day', 'hour')
dc = len(date_cols)

# %%

# This bit of code loads data from pre-processed pickle files or
# processes CSV files. This is the 'starter' data that the script needs.

# Load data about the cities. This is just for this example.
citytab = pd.read_csv(os.path.join('..', 'CityData.csv'),
                      dtype=dict(WMO=str, StCode=str))
# The following inputs would be input by the user.
# Station code.
stcode = citytab.StCode[22]
# Longitude.
stlong = 46.25
# Latitude.
stlat = 6.13
# Altitude.
stalt = 416.0

# See accompanying script "gw".
typicaldata, actualdata = gw.get_weather(stcode, citytab)

T = typicaldata.shape[0]

# %%

# Extract the relevant data from the loaded 'starter' data.
# This step won't be needed if only the relevant variables are in the
# input matrix.
# Conversely, the user might have to identify which time series are which.
xin_un = (np.vstack([typicaldata.index.dayofyear, typicaldata.hour.values,
                    typicaldata.tdb.values, typicaldata.tdp.values,
                    typicaldata.rh.values, typicaldata.ghi.values,
                    typicaldata.dni.values, typicaldata.dhi.values,
                    typicaldata.wspd.values, typicaldata.wdr.values])).T

# Number of data points read in.
N = xin_un.shape[0]
# Number of variables being considered =
# number of columns minus time variables.
D = xin_un.shape[1]

# The column order is as follows:
# day, hour, tdb, tdp, rh, ghi, dni, dhi, wspd, wdr

# %%

# Back to the generator.

# Change the value of variables *l_start* and *l_end* to alter
# the time for which synthetic time series should be generated.
# Generally leave to 0 and 365.

# Re-scale the hourly values (Normalize).
s_h = scaler.fit(xin_un)
xin_norm = scaler.transform(xin_un)

scaled_days = xin_norm[:, 0]
scaled_hours = xin_norm[:, 1]

start_time = time.monotonic()

# Length of history when "predicting" the next day's hourly profiles.
histlim = int(14*24)
# Start of time loop.
l_start = int(0*24)
# End of time loop.
l_end = int(365*24)
# Step size for time loop.
l_step = int(48)

# Pre-allocate a list of lists to store gp models for
# each day and variable.
listicle = [None]*D
gp_list = list()
for n in range(0, int(l_end/24)+1, 1):
    gp_list.append(copy.copy(listicle))

# Number of synthetic series to be created.
n_samples = int(histlim/l_step)

# To check the change in values from the previous step.
mae = np.zeros(len(column_names))

# Array to store the day-ahead 'predictions'.
pred_means = np.zeros([xin_norm.shape[0],
                       xin_norm.shape[1], n_samples])
# pred_stds = np.zeros_like(xin_norm_diff)

# Pre-allocate predictions matrices.
y_sample = np.zeros([l_step, n_samples])
y_std = np.zeros([l_step, n_samples])

# Hours of the year.
hh = typicaldata.hour + (typicaldata.index.dayofyear-1)*24

# Simply copy the input for the first l_start+l_step hours.
pred_means[0:l_start+l_step, :, :] = np.tile(np.atleast_3d(
    xin_norm[0:l_start+l_step, :]), n_samples)

# Loop over each day, starting from the second day.
for h in range(l_start+l_step, l_end+l_step, l_step):

    # Train initial model on the first day.
    # Take the last $histlim$ hours into account.
    slicer = range(l_start, h)

    if len(slicer) > histlim:
        b = int(len(slicer) - histlim)
        slicer = slicer[b:h+1]
        del b

    m = np.squeeze(typicaldata.month.values[h == hh])
    d = int(h/24)

    print("Month %d, day %d" % (m, d))

    # Concatenate the daily means with the hourly values, leaving
    # out the day variable from the daily means array (first column).
    train_slice = xin_norm[slicer, :]

    # Array to store stdev of predictions.
    # Array shape is l_step x number of variables - 2,
    # (since time variables do not have predictions).
    ystd_slice = np.zeros((l_step, xin_norm.shape[1]-2))

    # Cycle through all columns.
    for c, colname in enumerate(column_names):

        if colname in date_cols:
            # gp_h.append(None)
            continue

        else:
            # Find current output column.
            find_y = yfinder(colname)

            # Separate the training data (history) into x & y,
            # and reshape them to be 2D.
            x_train = train_slice[:, np.logical_not(find_y)]
            y_train = train_slice[:, find_y]

            # No test values.

        # End of columns loop.

    # Index for the next day's data.
    tomorrows_slice = range(h, h+l_step)

    # If tomorrow's slice is beyond the 31st of December,
    # then break this loop.
    if (h+l_step) > T:
        break

    # The next day's values ($histlim$), which will be updated by the
    # Gibbs sampling loop.
    tomorrows_vals = np.zeros((l_step, xin_norm.shape[1], n_samples))

    for c, colname in enumerate(column_names):

        if colname in date_cols:
            # Use original values of time variables.
            tomorrows_vals[:, c, :] = np.tile(np.atleast_2d(
                xin_norm[tomorrows_slice, c]).T, n_samples)
        else:
            # Find the current column in the overall array.
            find_y = yfinder(colname)

            # Separate the training data (history) into x & y.
            x_train = train_slice[:, np.logical_not(find_y)]
            y_train = train_slice[:, find_y]

            # Call
            gp_temp, ypred, ystd = gp_funcs(colname, x_train, y_train)
            # The second and third outputs are dummies,
            # though calling them stupid would be mean.
            # Delete them immediately.
            del ypred, ystd

            x_query = x_train[-l_step:, :]

            x_pred = np.squeeze(gp_temp.sample_y(x_query,
                    n_samples=n_samples), axis = 1)

    # Use the GP of each day to predict that day,
    # but add the predictions from using the GPs of
    # surrounding days as well.

            tomorrows_vals[:, c, :] = x_pred

#            # Day in recent history which will be used as a query point for GP.
#            dd = 0
#
#            for xq in range(0, x_train.shape[0], l_step):
#                x_pred = x_train[xq:xq+l_step, :]
#
#                # Predict the current variable for the next day at x_pred.
#                tomorrows_vals[:, c, dd] = np.reshape(
#                        gp_temp.predict(x_pred),
#                        tomorrows_vals[:, c, dd].shape)
#                dd += 1
#            del dd

            gp_list[d][c] = gp_temp

    pred_means[tomorrows_slice, :, :] = tomorrows_vals

    # End of hourly loop.

end_time = time.monotonic()
print("Time taken to train month-wise hourly models was %6.2f seconds."
      % (end_time - start_time))

# Save list to pickle file.
with open('gp_list', 'wb') as fp:
    pickle.dump(gp_list, fp)

# Un-scale the series.
# Create a ndarray like pred_means.
xout_un = np.zeros_like(pred_means)
# Inverse transform all columns together, for each sample.
for n in range(0, n_samples):
    xout_un[:, :, n] = s_h.inverse_transform(pred_means[:, :, n])

# Synthetic solar data requires post-processing.
for c, colname in enumerate(column_names):

    if colname in ['ghi', 'dni', 'dhi']:

        for n in range(0, n_samples):
            # Using the source data - check to see if there
            # should be sunlight at a given hour. If not,
            # then set corresponding synthetic value to zero.

            temp = xout_un[:, c, n]
            temp[xin_un[:, c]==0] = 0
            xout_un[:, c, n] = temp

            # If there is a negative value (usually at sunrise
            # or sunset), interpolate.
            temp = xout_un[:, c, n]
            temp[temp<0] = np.nan
            nans = np.isnan(temp)
            temp[nans] = np.interp(xin_un[nans, 1],
                xin_un[~nans, 1], temp[~nans])
            xout_un[:, c, n] = temp

# A potential improvement would be to calculate sunrise and sunset
# independently since that is an almost deterministic calculation.

# Save synthetic time series.
for n in range(0, n_samples):
    filepath = os.path.join(path_fldr_csv, 'ts_syn_%s.csv' % n)
    np.savetxt(filepath, np.squeeze(xout_un[:, :, n]), '%6.2f',
               delimiter=',', header=','.join(column_names))

# The generation part ends here - the rest is just plotting various things.

#
# %%
#
# =============================================================================
#
# # This bit is only for plotting (to see how generated values correspond
# # to summary statistics) and has nothing to do with the generator.
#
# # Calculate the iqr and median to plot later.
# actualq3 = gw.weather_stats(actualdata, ['month'], 'q3')
# actualmed = gw.weather_stats(actualdata, ['month'], 'med')
# actualq1 = gw.weather_stats(actualdata, ['month'], 'q1')
#
# temp = np.repeat(actualq1.values,
#                  int(typicaldata.shape[0]/actualq1.shape[0]), axis=0)
# temp2 = pd.DataFrame(data=temp, columns=actualq1.columns)
# actualq1 = temp2
# del (temp, temp2)
#
# temp = np.repeat(actualq3.values,
#                  int(typicaldata.shape[0]/actualq3.shape[0]), axis=0)
# temp2 = pd.DataFrame(data=temp, columns=actualq3.columns)
# actualq3 = temp2
# del (temp, temp2)
#
# temp = np.repeat(actualmed.values,
#                  int(typicaldata.shape[0]/actualmed.shape[0]), axis=0)
# temp2 = pd.DataFrame(data=temp, columns=actualmed.columns)
# actualmed = temp2
# del (temp, temp2)

# ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
#                facecolor='w', edgecolor='k').add_subplot(111)
# plt.plot(actualq1.tdb.values)
# plt.plot(actualq3.tdb.values)
# plt.show

# =============================================================================

# %%

plotrange = range(l_start, l_end)
p = 2
n = range(0, n_samples+1)

# Line plot.
ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                facecolor='w', edgecolor='k').add_subplot(111)

p1 = plt.plot(plotrange, xin_un[plotrange, p],
              color=colours.orange, linewidth=1.5, zorder=1)
p2 = plt.plot(plotrange, xout_un[plotrange, p, :], linewidth=0.5,
              alpha=1)
#p31 = plt.plot(plotrange, actualq1.tdb.values, linewidth=2,
#               alpha=1, color=colours.blackest, zorder = 3, marker='o')
#p32 = plt.plot(plotrange, actualmed.tdb.values, linewidth=2,
#               alpha=1, color=colours.blackest, zorder = 4, marker='o')
#p33 = plt.plot(plotrange, actualq3.tdb.values, linewidth=2,
#               alpha=1, color=colours.blackest, zorder = 5, marker='o')

# ax.set_xticks(range(int(l_start), int(l_end)+l_step, l_step*4)) color=colours.orange,
ax.grid()
plt.xlim(plotrange[0], plotrange[-1])
plt.xlabel('Hour')
plt.title(column_names[p])
plt.ylabel('[degC]') # Change according to plotted variable.
plt.legend(['Recorded', 'Synthetic']) # , 'q1', 'q2','q3'])

figname = os.path.join(figpath, 'line_{0}.pdf'.format('tdb'))
#plt.savefig(figname)
plt.show()

# %%

p = 2
n = 0  # random.randint(0, n_samples)
plotrange = range(l_start, l_end)

# ACF plot.
ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                facecolor='w', edgecolor='k').add_subplot(111)

p1 = plot_acf(xin_un[plotrange, p], ax=ax, lags=72, alpha=0.05)
p1 = plot_acf(xout_un[plotrange, p, n], ax=ax, lags=72, alpha=0.05)

plt.xlabel('Lags [Hours]')
plt.ylabel('ACF')
plt.title('TDB')
plt.legend(['', 'Recorded', '', 'Synthetic'])

# PACF plot.
ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                facecolor='w', edgecolor='k').add_subplot(111)

p1 = plot_pacf(xin_un[plotrange, p], ax=ax, lags=72, alpha=0.05)
p1 = plot_pacf(xout_un[plotrange, p, n], ax=ax, lags=72, alpha=0.05)

plt.xlabel('Lags [Hours]')
plt.ylabel('PACF')
plt.title('TDB')
plt.legend(['', 'Recorded', '', 'Synthetic'])

ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                facecolor='w', edgecolor='k').add_subplot(111)

p1 = plt.scatter(xin_un[plotrange, 2], xin_un[plotrange, 3], color = colours.blackest)
p1 = plt.scatter(xout_un[plotrange, 2, n], xout_un[plotrange, 3, n],  color=colours.reddest)

plt.xlabel('Lags [Hours]')
plt.ylabel('Corr.')
plt.title('TDB vs TDP')
plt.legend(['Recorded', 'Synthetic'])
# %%

# Repeat the RGB spec for grey n_sample times to plot n_sample synthetic series.
grey_reps = np.repeat(np.atleast_2d(colours.grey).T,n_samples, axis=1).T

# Histogram.
ax = plt.figure(num=None, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k').add_subplot(111)


p1 = plt.hist(xin_un[plotrange, p], normed=True, cumulative=False,
              bins=50, histtype='step', color=colours.blackest, linewidth=2,
              zorder=1)
p2 = plt.hist(xout_un[plotrange, p, :], normed=True, cumulative=False,
              bins=50, histtype='step', color=grey_reps)
# =============================================================================
# p3 = plt.hist(actualdata.tdb.values[np.logical_not(np.isnan(
#         actualdata.tdb.values))], normed=True, cumulative=False,
#               bins=50, histtype='step', zorder=2, linewidth=2, color=colours.orange)
# =============================================================================

plt.xlabel('Temp.')
plt.ylabel('CDF')
plt.legend(['Recorded', 'Synthetic'], loc='upper right')

figname = os.path.join(figpath, 'hist_{0}.pdf'.format('tdb'))
plt.savefig(figname)
plt.show()
# %% TDP

plotrange = range(l_start, 48)
p = 5

ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                facecolor='w', edgecolor='k').add_subplot(111)

p1 = plt.plot(plotrange, xin_un[plotrange, p],
              color=colours.orange, linewidth=0.5, zorder=n_samples+1)
p2 = plt.plot(plotrange, xout_un[plotrange, p, :], linewidth=0.75,
              alpha=1, color=colours.grey, zorder=1)
#p2 = plt.plot(plotrange, actualstats.tdb, linewidth=0.75,
#              alpha=1, color=colours.lgrey)

# ax.set_xticks(range(int(l_start), int(l_end)+l_step, l_step*4))
# labels = [item.get_text() for item in ax.get_xticklabels()]
# labels = range(int(l_start/24), int(l_end/24)+1, 4)
# ax.set_xticklabels(labels)
ax.grid()
# plt.xlim(l_start, l_end)
plt.xlabel('Hour')
plt.ylabel('[degC]')
plt.title(column_names[p])
plt.legend(['Recorded', 'Synthetic'])

figname = os.path.join(figpath, 'line_{0}.pdf'.format('tdp'))
plt.savefig(figname)
plt.show()

# %%

# Histogram.
plt.figure(num=None, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k')

p1 = plt.hist(xin_un[plotrange, p], normed=True, cumulative=False,
              bins=50, histtype='step', color=colours.orange, linewidth=2,
              zorder=1)
p2 = plt.hist(xout_un[plotrange, p, :], normed=True, cumulative=False,
              bins=50, histtype='step', zorder=n_samples+1)

plt.xlabel('Temp.')
plt.ylabel('CDF')
plt.legend(['Recorded', 'Synthetic'], loc='upper right')

figname = os.path.join(figpath, 'hist_{0}.pdf'.format('tdp'))
plt.savefig(figname)
plt.show()
