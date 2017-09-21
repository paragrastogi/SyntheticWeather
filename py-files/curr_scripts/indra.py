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

import os
import random as random

import matplotlib.pyplot as plt
# import matplotlib.mlab as mlab
import numpy as np
import pandas as pd

# import scipy as sp
# from scipy import signal

import default_colours as colours

# To time the script.
#import time

# Lets you copy objects if you need two
# distinct objects instead of the usual pointers.
#import copy

# from sklearn.neural_network import MLPRegressor
# from sklearn.gaussian_process import GaussianProcessRegressor

# from sklearn.model_selection import GridSearchCV
# from sklearn.model_selection import RandomizedSearchCV
# from sklearn.model_selection import KFold

# from skopt import gp_minimize
# from sklearn.datasets import fetch_mldata
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import Imputer

# import statsmodels as sm
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.graphics.tsaplots import plot_pacf
# Not needed if you don't want to plot autorrrelation plots.
# Could also use pandas:
# from pandas.tools.plotting import autocorrelation_plot


# These custom functions load and clean recorded data.
# For now, we are only concerned with ncdc and nsrdb.
import get_weather_data as gw
# from get_weather_data import load_actual
# from get_weather_data import load_typical

# Custom functions to calculate error metrics.
# import losses.
from losses import rmseloss
from losses import maeloss

# import les petites fonctionnes utiles.
from petites import yfinder
from petites import fitgp
from petites import learngp
from petites import samplegp

# For storing data by pickling.
#import pickle

# importlib is only needed if libraries need to be re-imported mid-script,
# which would usually happen only during debugging.
# import importlib as im
# from importlib import reload

__author__ = 'Parag Rastogi'

# Disbale this command if you want figures in a new window.
#get_ipython().magic('matplotlib inline')
# Also disable when running in command window.

def synweather(seedfile, stcode='GLA',
               path_wthr_fldr='/usr/esru/esp-r/climate',
               outpath=os.getcwd(), n_samples=100,
               l_start=int(0), l_end=int(31*24),
               l_step=int(4*24), histlim=int(14*24),
               stlat=0.0, stlong=0.0, stalt=0.0,
               path_fldr_pickle=os.getcwd()):

    # ------------------
    # Some initialisation house work.

    # Seed random number generators.
    seed = 8760
    np.random.seed(seed)
    random.seed = seed

    # This is the master tuple of column names, which should not be modified.
    # If it is modified, then the same tuple in in the script called
    # petites.py  should also be updated.
    column_names = ('day', 'hour', 'tdb', 'tdp', 'rh', 'ghi', 'dni', 'dhi',
                    'wspd', 'wdr')
    # This is the master tuple of time variable names, which should also not
    # be modified.
    #    date_cols = ('day', 'hour')
    # dc = len(date_cols)

    # ----------------

    # The program needs to know the mode you are running in.
    # There are two modes supported:
    # 1. With starter data - seedfile = True
    # 2. Without started data - seedfile = False

    if seedfile:

        # The script needs some 'starter' data.

        currentdata = None

        try:
            currentdata, _ = gw.load_typical(
                    path_wthr_fldr, stcode, force=True)
        except:
            print('Could not find the files in %s with station ' +
                  'code %s or could not read the file I found.'
                  % path_wthr_fldr, stcode)

        # If there is no path to a weather file, then this bit of code loads
        # data from pre-processed pickle files or processes CSV files.

        if (currentdata is None):

            # Load data about the cities. This is just for this example.
            citytab = pd.read_csv(os.path.join('..', 'CityData.csv'),
                                  dtype=dict(WMO=str, StCode=str))
            # Specify the sources of the actual data -
            # please follow AMY keywords list.
            sources = ('ncdc', 'nsrdb', 'nrel_indiasolar',
                       'ms', 'WY2', 'nasa_saudi')

            # See accompanying script "gw".
            currentdata, historicaldata = gw.get_weather(
                    stcode, citytab, sources)

        # T = currentdata.shape[0]

        # Extract the relevant data from the loaded 'starter' data.
        # This step won't be needed if only the relevant variables are
        # in the input matrix. Conversely, the user might have
        # to identify the relevant time series.
        ts_curr_in = (np.vstack([
                currentdata.index.dayofyear, currentdata.hour.values,
                currentdata.tdb.values, currentdata.tdp.values,
                currentdata.rh.values, currentdata.ghi.values,
                currentdata.dni.values, currentdata.dhi.values,
                currentdata.wspd.values, currentdata.wdr.values])).T

        # ts_hist_norm = scaler.transform(ts_hist_in)

        # ts_hist_in = (np.vstack([
        #        historicaldata.index.dayofyear, historicaldata.hour.values,
        #        historicaldata.tdb.values, historicaldata.tdp.values,
        #        historicaldata.rh.values, historicaldata.ghi.values,
        #        historicaldata.dni.values, historicaldata.dhi.values,
        #        historicaldata.wspd.values, historicaldata.wdr.values])).T

        # Number of data points read in.
        # histlen = ts_curr_in.shape[0]

        # The column order is as follows:
        # day, hour, tdb, tdp, rh, ghi, dni, dhi, wspd, wdr

    # end seedfile if statement.
    #

    # Where do you want to save figures?
    figpath = 'figures'
    # ... and the csv outputs?
    # outpath = os.path.join('/home', 'rasto', 'Documents', 'WeatherData',
    #                                 'SyntheticData', 'csv_saves')

    # Check if these folders exist.
    # If they don't, create them.
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    if not os.path.isdir(figpath):
        os.mkdir(figpath)

    # %%

    # This is the part of the generator where GP models are fitted to
    # the incoming data. Naturally, it is only invoked if the function
    # was invoked with starter data.

    gp_list, month_tracker = learngp(l_start, l_end, l_step, histlim,
                       ts_curr_in)

    # %%

    picklepath = os.path.join(path_fldr_pickle, 'ts_syn_%s.p' % stcode)

    xout_un = samplegp(gp_list, l_start, l_end, l_step, histlim, n_samples,
              ts_curr_in, month_tracker, picklepath, outpath)

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
    # actualq3 = gw.weather_stats(historicaldata, ['month'], 'q3')
    # actualmed = gw.weather_stats(historicaldata, ['month'], 'med')
    # actualq1 = gw.weather_stats(historicaldata, ['month'], 'q1')
    #
    # temp = np.repeat(actualq1.values,
    #                  int(currentdata.shape[0]/actualq1.shape[0]), axis=0)
    # temp2 = pd.DataFrame(data=temp, columns=actualq1.columns)
    # actualq1 = temp2
    # del (temp, temp2)
    #
    # temp = np.repeat(actualq3.values,
    #                  int(currentdata.shape[0]/actualq3.shape[0]), axis=0)
    # temp2 = pd.DataFrame(data=temp, columns=actualq3.columns)
    # actualq3 = temp2
    # del (temp, temp2)
    #
    # temp = np.repeat(actualmed.values,
    #                  int(currentdata.shape[0]/actualmed.shape[0]), axis=0)
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
    p = 3
    n = range(0, n_samples+1)

    # Line plot.
    ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                    facecolor='w', edgecolor='k').add_subplot(111)

    p1 = plt.plot(plotrange, ts_curr_in[plotrange, p],
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

    p1 = plot_acf(ts_curr_in[plotrange, p], ax=ax, lags=72, alpha=0.05)
    p1 = plot_acf(xout_un[plotrange, p, n], ax=ax, lags=72, alpha=0.05)

    plt.xlabel('Lags [Hours]')
    plt.ylabel('ACF')
    plt.title('TDB')
    plt.legend(['', 'Recorded', '', 'Synthetic'])

    # PACF plot.
    ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                    facecolor='w', edgecolor='k').add_subplot(111)

    p1 = plot_pacf(ts_curr_in[plotrange, p], ax=ax, lags=72, alpha=0.05)
    p1 = plot_pacf(xout_un[plotrange, p, n], ax=ax, lags=72, alpha=0.05)

    plt.xlabel('Lags [Hours]')
    plt.ylabel('PACF')
    plt.title('TDB')
    plt.legend(['', 'Recorded', '', 'Synthetic'])

    ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                    facecolor='w', edgecolor='k').add_subplot(111)

    p1 = plt.scatter(ts_curr_in[plotrange, 2], ts_curr_in[plotrange, 3], color = colours.blackest)
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


    p1 = plt.hist(ts_curr_in[plotrange, p], normed=True, cumulative=False,
                  bins=50, histtype='step', color=colours.blackest, linewidth=2,
                  zorder=1)
    p2 = plt.hist(xout_un[plotrange, p, :], normed=True, cumulative=False,
                  bins=50, histtype='step', color=grey_reps)
    # =============================================================================
    # p3 = plt.hist(historicaldata.tdb.values[np.logical_not(np.isnan(
    #         historicaldata.tdb.values))], normed=True, cumulative=False,
    #               bins=50, histtype='step', zorder=2, linewidth=2, color=colours.orange)
    # =============================================================================

    plt.xlabel('Temp.')
    plt.ylabel('CDF')
    plt.legend(['Recorded', 'Synthetic'], loc='upper right')

    figname = os.path.join(figpath, 'hist_{0}.pdf'.format('tdb'))
    plt.savefig(figname)
    plt.show()
    # %% TDP

    plotrange = range(l_start, l_end)
    p = 5

    ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                    facecolor='w', edgecolor='k').add_subplot(111)

    p1 = plt.plot(plotrange, ts_curr_in[plotrange, p],
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

    p1 = plt.hist(ts_curr_in[plotrange, p], normed=True, cumulative=False,
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
