#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 11:47:58 2017

@author: rasto

Translating, as faithfully as possible, the resampling method first proposed
my thesis (Rastogi, 2016, EPFL).
"""

import pickle
import copy

import numpy as np
import pandas as pd

from scipy.optimize import curve_fit

import fourier
from ts_models import select_models
# Useful small functions like solarcleaner.
import petites as petite

# Number of variables resampled - TDB and RH.
NUM_VARS = 2

# "Standard" length of output year.
STD_LEN_OUT = 8760

# Number of nearest neighbours (of each synthetic day in array of recorded
# days using daily mean temperature) from which to choose solar radiation.
# NUM_NBOURS = 10

# This is the master tuple of column names, which should not
# be modified.
# COLUMNS = ('year', 'month', 'day', 'hour', 'tdb', 'tdp', 'rh',
#            'ghi', 'dni', 'dhi', 'wspd', 'wdr')


def resampling(xy_train, train=True, n_samples=10,
               picklepath='test_out.npy', counter=0,
               arma_params=None, bounds=None,
               cc_data=None, cc_scenario=None):

    # Reassign defaults if incoming list params are None
    # (i.e., nothing passed.)
    if arma_params is None:
        arma_params = [2, 2, 1, 1, 24]

    if bounds is None:
        bounds = [0.01, 99.9]

    # if train is True, then the model must be fit to incoming data.
    if train:

        ffit, selmdl, xout = trainer(xy_train, n_samples, arma_params,
                                     bounds, cc_data, cc_scenario)

        # Save the outputs as a pickle.

        with open(picklepath, 'wb') as f:
            pickle.dump(xout, f)

        sample = None

    else:

        # The non-training call is only meant to return a sample.
        ffit = None
        selmdl = None

        sample = sampler(picklepath, counter)

    # Return all three outputs.

    return ffit, selmdl, sample


def trainer(xy_train, n_samples, arma_params, bounds, cc_data, cc_scenario):
    '''Train the model with this function.'''

    # Save a copy of all data to calculate quantiles later.
    xy_train_all = xy_train

    # Get the unique years. Remove the last element since that is 2223,
    # which exists because Python interprets the last hour of a typical
    # year (2222) as the first hour of the next year.
    all_years = np.unique(xy_train.index.year)[:-1]

    xy_train = pd.DataFrame()

    while xy_train.shape[0] < 8760:

        select_year = all_years[np.random.randint(0, len(all_years), 1)[0]]

        # Keep only that one year of data.
        xy_train = xy_train_all[str(select_year) + '-01-01':
                                str(select_year) + '-12-31']

        xy_train = petite.remove_leap_day(xy_train)

        if xy_train.shape[0] > 8760:
            xy_train = xy_train.iloc[0:STD_LEN_OUT, :]

    x_calc_params = np.arange(0, xy_train_all.shape[0])
    x_fit_models = np.arange(0, STD_LEN_OUT)

    # Fit fourier functions to the tdb and rh series.

    # The curve_fit function outputs two things: parameters of the fit and
    # the estimated covariance. We only use the first.
    # Inputs are the function to fit (fourier in this case),
    # xdata, and ydata.
    # Use all the data available to calculate these parameters.
    params = [curve_fit(fourier.fit_tdb, x_calc_params, xy_train_all['tdb']),
              curve_fit(fourier.fit_rh, x_calc_params, xy_train_all['rh'])]

    # Call the fourier fit function with the calculated
    # parameters to get the values of the fourier fit at each time step
    ffit = [fourier.fit('tdb', x_fit_models, *params[0][0]),
            fourier.fit('rh', x_fit_models, *params[1][0])]

    # ARIMA models cannot be fit to overly long time series.
    # Select a random year from all the data available.
    # unique_years = np.unique(xy_train.index.year)
    # selected_year = unique_years[np.random.randint(0, len(unique_years))]

    #  Filtered_xy = xy_train[xy_train.index.year==selected_year]
    # Filtered_DM = DeMeaned[DeMeaned.index.year==selected_year]

    # Now subtract the low- and high-frequency fourier fits
    # (whichever is applicable) from the raw values to get the
    # 'de-meaned' values (values from which the mean has
    # been removed).

    DeMeaned = pd.concat([x - y for x, y in
                          zip([xy_train["tdb"], xy_train["rh"]], ffit)],
                         axis=1)
    DeMeaned.index = xy_train.index

    # Fit ARIMA models.

    selmdl = list()
    resid = np.zeros([DeMeaned["tdb"].shape[0], NUM_VARS])

    for idx, ser in enumerate(DeMeaned):
        mdl_temp, resid[:, idx] = select_models(
            arma_params, DeMeaned[ser])
        selmdl.append(mdl_temp)

    print("Done with fitting models to TDB and RH.\r\n")

    print(("Simulating the learnt model to get synthetic noise series."
           "This might take some time.\r\n"))

    # num_years = int(xy_train.shape[0] / STD_LEN_OUT)

    resampled = np.zeros([STD_LEN_OUT, NUM_VARS, n_samples])

    for midx, mdl in enumerate(selmdl):
        for sample_num in range(0, n_samples):
            resampled[:, midx, sample_num] = mdl.simulate(
                nsimulations=STD_LEN_OUT)
        # End n for loop.
    # End mdl for loop.

    # Add the resampled time series back to the fourier series.

    if cc_data is None:
        # First make the xout array using all variables. Variables other
        # than RH and TDB are just repeated from the incoming files.
        xout = list()

        # Clean the generated temperature values using extreme percentiles
        # as proxies for 'outliers'.

        # Clean the RH values using phyiscal limits (0-100).

        # Add the fourier fits from the training data to the
        # resampled/resimulated ARMA model outputs.

        for nidx in range(0, n_samples):

            # Copy the master datatable of all values.
            xout_temp = copy.deepcopy(xy_train)

            for idx, var in enumerate(DeMeaned[["tdb", "rh"]]):

                # Replace only var (tdb or rh).
                # Also send it to the quantile cleaner.
                xout_temp[var] = petite.quantilecleaner(
                    (resampled[:, idx, nidx] + ffit[idx]), xy_train_all,
                    var, bounds=bounds)

            xout.append(xout_temp)

    else:

        # ccindex = cc_data[cc_scenario].index.get_level_values(1)
        cc_models = set(cc_data[cc_scenario].index.get_level_values(0))
        xout = list()  # ([xy_train] * n_samples)

#        for idx, var in enumerate(DeMeaned[["tdb", "rh"]]):
        for model in cc_models:

            this_cc_out = cc_data[cc_scenario].loc[model]
            gcm_years = np.unique(cc_data[cc_scenario].loc[model].index.year)

            for yidx, future_year in enumerate(gcm_years):

                # Select only this year of cc model outputs.
                cctable = this_cc_out[str(future_year) + '-01-01':
                                      str(future_year) + '-12-31']
                # leap_idx = [idx for idx, x in enumerate(cctable.index)
                # if x == pd.to_datetime(str(future_year) + '-02-29 12:00:00')]
                # cctable = cctable.drop(cctable.index[leap_idx])
                cctable = petite.remove_leap_day(cctable)

                if cctable.shape[0] < 365:
                    continue

                for nidx in range(0, n_samples):

                    xout_temp = copy.deepcopy(xy_train)

                    for idx, var in enumerate(["tdb", "rh"]):

                        tas = cctable["tas"].values

                        if var == "tdb":
                            ccvar = np.repeat(tas, [24], axis=0)
                        elif var == "rh":
                            huss = cctable["huss"].values
                            ps = cctable["ps"].values

                            # Convert specific humifity to humidity ratio.
                            w = -huss / (huss - 1)

                            # Convert humidity ratio (w) to
                            # Relative Humidity (RH).
                            rh = petite.w2rh(w, tas, ps)

                            # Is there some way to replace the fourier fit at
                            # a finer grain instead of repeating the daily
                            # mean value 24 times?
                            ccvar = np.repeat(rh, [24], axis=0)

                        xout_temp[var] = resampled[:, idx, nidx] + ccvar
                        future_index = pd.DatetimeIndex(
                            start=str(future_year) + "-01-01 00:00:00",
                            end=str(future_year) + "-12-31 23:00:00",
                            freq='1H')
                        # Remove leap days.
                        future_index = future_index[
                            ~((future_index.month == 2) &
                              (future_index.day == 29))]
                        xout_temp.index = future_index

                    xout.append(xout_temp)

    # End for loop.

    # End loop over samples.

    # Calculate TDP.

# %%
    for idx, df in enumerate(xout):

        df["tdp"] = petite.calc_tdp(df["tdb"], df["rh"])
        xout[idx] = df

    # Calculate daily means of temperature.
    mean_list = list()
    for df in xout:
        df_dm = df['tdb'].resample('1D').mean()

        if len(df_dm) > 365:
            df_dm = petite.remove_leap_day(df_dm)

        mean_list.append(df_dm)

    tdb_dailymeans = dict(syn=mean_list,
                          rec=xy_train['tdb'].resample('1D').mean())

    sol_idx = [x for x, y in enumerate(xy_train)
               if y in ['ghi', 'dhi', 'dni']]

    # Number of nearest neighbours to keep when varying solar quantities.
    nn_top = 10

    for this_month in range(1, 13):

        # This month's indices.
        idx_this_month = tdb_dailymeans['rec'].index.month == this_month

        nearest_nbours = list()

        # Cycle through each array of daily means.
        for _, syn_means in enumerate(tdb_dailymeans['syn']):
            # Find the nearest neighbour to each element of the syn_means
            # array. Nearest neighbours must be in the same month to preserve
            # length of day.

            # First find 10 nearest neighbours.
            nearest_nbours_temp = [(np.abs(
                x - tdb_dailymeans['rec'][idx_this_month])).argsort()[:nn_top]
                for x in syn_means[idx_this_month]]
            # This works by sorting the values in ascending order, taking the
            # indices, and picking the first nn_top.

            # Pick one of the ten nearest neighbours.
            nearest_nbours.append([x[np.random.randint(
                0, nn_top, size=1)] for x in nearest_nbours_temp])

        # Find the solar data for this month.
        idx_solar_this_month = np.asarray(
            [xy_train.iloc[xy_train.index.month == this_month, x]
             for x in sol_idx]).T
        # Reshape into day-sized blocks.
        idx_solar_this_month = np.reshape(idx_solar_this_month, [-1, 24, 3])

        solar_samples = np.zeros([len(nearest_nbours),
                                  np.sum(idx_this_month),
                                  24, len(sol_idx)])

        for ng_idx, n_group in enumerate(nearest_nbours):
            for sg_idx, n_subgroup in enumerate(n_group):
                solar_samples[ng_idx, sg_idx, :, :] = np.array(
                    [idx_solar_this_month[x, :, :].tolist()
                     for x in n_subgroup][0])

        # Send the solar columns to cleaning.
        for nidx in range(0, n_samples):

            # Put the solar samples back in to xout.
            for sidx, solcol in enumerate(sol_idx):
                xout[nidx].iloc[:, sidx].mask(
                    xout[nidx].index.month == this_month,
                    other=pd.Series(np.reshape(
                        solar_samples[nidx, :, :, sidx], [-1])))

    # End month loop.

    for nidx in range(0, n_samples):
        xout[nidx].iloc[:, solcol] = petite.solarcleaner(
            xout[nidx].iloc[:, solcol], xy_train.iloc[:, solcol])

    # End nidx loop.

    return ffit, selmdl, xout


def sampler(picklepath, counter):
    '''Only opens the pickle of saved samples and returns ONE sample.'''

    try:

        xout = pickle.load(open(picklepath, 'rb'))
        sample = xout[counter]

    except AttributeError:

        print("I could not open the pickle file with samples. " +
              "Please check it exists at {0}.".format(picklepath))
        sample = None

    return sample
