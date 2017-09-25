#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 16:23:50 2017

"""
__author__ = 'Parag Rastogi'

import os
import time
import copy
# import pickle
import random
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.preprocessing import StandardScaler
import pandas as pd

# This is the master tuple of column names, which should not be modified.
column_names = ('month', 'day', 'hour', 'tdb', 'tdp', 'rh',
                'ghi', 'dni', 'dhi', 'wspd', 'wdr')

# This is the master tuple of time variable names,
# which should also not be modified.
date_cols = ('month', 'day', 'hour')
dc = len(date_cols)

# Make a scaler - either standard (to scale variables to \mu = 0 and \sigma = 1) or robust (with median and iqr).
scaler = StandardScaler()

# Alternate (more complicated) kernels:
# With seasonal periodicity (ExpSineSquared) and allow it to
# decay away from exact periodicity (RBF).
# Allow for a rising/falling trend with a separate RBF kernel and noise
# with a WhiteKernel.


def fitgp(series, xtrain, ytrain, xtest=np.empty(0), *args):

    if len(args) > 0:
        gp_in = args[0]

    # Initial values for hyperparameters.
    hypguess = np.random.random(xtrain.shape[1])

    # n_restarts_optimizer selects "result from run with minimal
    # (negative) log-marginal likelihood", line 232-236 in
    # github.com/scikit-learn/scikit-learn/blob/master/sklearn/
    # gaussian_process/gpr.py#L206 .
    n_restarts = 5

    # WhiteKernel represents signal noise ($\sigma_n$).
    if series in ['tdb_dm', 'ghi_dm', 'tdp_dm']:
        kd = RBF(length_scale=hypguess) + \
            WhiteKernel(noise_level=1, noise_level_bounds=[0.9, 1.1])
    elif series in column_names:
        kd = RBF(length_scale=hypguess) + \
            WhiteKernel(noise_level=0.1, noise_level_bounds=[1e-2, 0.1])

    if len(args) == 0:
        # No incoming regressor object, so learn hyper-parameters.
        gp_out = GaussianProcessRegressor(kernel=kd, alpha=0,
                                          n_restarts_optimizer=n_restarts,
                                          optimizer='fmin_l_bfgs_b',
                                          normalize_y=False)

        gp_out.fit(xtrain, ytrain)
    else:
        # Do not retrain kernel.
        gp_out = gp_in
        gp_out.fit(xtrain, ytrain)

    #     print("GPML kernel: %s" % gp.kernel_)
    #     print("Log-marginal-likelihood: %.3f"
    #           % gp.log_marginal_likelihood(gp.kernel_.theta))

    if xtest.size > 0:
        y_pred, y_std = gp_out.predict(xtest, return_std=True)

    elif xtest.size == 0:
        y_pred = np.zeros_like(ytrain)
        y_std = np.zeros_like(ytrain)

    # For some reason, y_std comes out as 1D [365, ] array.
    # Convert it to 2d and transpose (column vector).
    if len(y_std.shape) == 1:
        y_std = (np.atleast_2d(y_std)).T

    return gp_out, y_pred, y_std

    # %%


def learngp(l_start, l_end, l_step, histlim,
            ts_curr_in, path_fldr_pickle=os.getcwd()):

    start_time = time.monotonic()

    # Re-scale the hourly values (Normalize).
    s_c = scaler.fit(ts_curr_in)
    # s_h = scaler.fit(ts_hist_in)
    ts_curr_norm = s_c.transform(ts_curr_in)

    # Number of variables being considered =
    # number of columns minus time variables.
    metvars = ts_curr_norm.shape[1]

    # Pre-allocate a list of lists to store gp models for
    # each day and variable.
    listicle = [None]*int(((l_end-l_start)/l_step)+1)
    gp_list = list()
    for n in range(0, metvars, 1):
        gp_list.append(copy.copy(listicle))

    # To check the change in values from the previous step.
    # mae = np.zeros(len(column_names))

    # Pre-allocate predictions matrices.
    # y_sample = np.zeros([l_step, n_samples])
    # y_std = np.zeros([l_step, n_samples])

    # Hours of the year - calculated from hours of the day and day of the year.
    hh = ts_curr_in[:, 2] + (ts_curr_in[:, 1] - 1)*24

    d = 0

    month_tracker = np.zeros(len(range(l_start+l_step, l_end+l_step, l_step)))

    # Loop over each day, starting from the second day.
    for h in range(l_start+l_step, l_end+l_step, l_step):

        # Train initial model on the first day.
        # Take the last $histlim$ hours into account.
        slicer = np.arange(l_start, h)

        if len(slicer) > histlim:
            b = int(len(slicer) - histlim)
            slicer = slicer[b:h+1]
            del b

        date_slice = [(mslice, dslice, hslice)
                      for (mslice, dslice, hslice) in
                      zip(ts_curr_in[:, 0], ts_curr_in[:, 1], hh)
                      if hslice in slicer]

        # If tomorrow's slice is beyond the 31st of December,
        # then break this loop.
        if (h+l_step) > ts_curr_norm.shape[0]:
            break

        # Output the current month and day.
        month_tracker[d] = int(np.unique(np.squeeze(
                ts_curr_in[h == hh, 0])))
        # Change the unique function here to accomodate
        # multi-year records.

        print("Model no. {0}, month {1}, day {2}".format(
                d, month_tracker[d], int(h/24)))

        # Concatenate the daily means with the hourly values, leaving
        # out the day variable from the daily means array (first column).
        train_slice = ts_curr_norm[slicer, :]
        # USE dateslice here instead of slicer.

        # Cycle through all columns.
        for c, colname in enumerate(column_names):

            if colname in date_cols:
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

        for c, colname in enumerate(column_names):

            if colname in date_cols:
                # Use original values of time variables.
                # tomorrows_vals[:, c, :] = np.tile(np.atleast_2d(
                #     ts_curr_norm[tomorrows_slice, c]).T, n_samples)
                continue

            else:
                # Find the current column in the overall array.
                find_y = yfinder(colname)

                # Separate the training data (history) into x & y.
                x_train = train_slice[:, np.logical_not(find_y)]
                y_train = train_slice[:, find_y]

                # Call
                gp_temp, _, _ = fitgp(colname, x_train, y_train)
                # The second and third outputs are dummies,
                # though calling them stupid would be mean.

                gp_list[c][d] = gp_temp

        # Update 'training day' iterator.
        d += 1

        # End 'h' loop.

    end_time = time.monotonic()
    print("Time taken to train models was %6.2f seconds."
          % (end_time - start_time))

    return gp_list, month_tracker


def samplegp(gp_list, l_start, l_end, l_step, histlim, n_samples,
             ts_curr_in, month_tracker,
             picklepath=os.getcwd(), outpath=os.getcwd()):

    start_time = time.monotonic()

    # Re-scale the hourly values (Normalize).
    s_c = scaler.fit(ts_curr_in)
    # s_h = scaler.fit(ts_hist_in)
    ts_curr_norm = s_c.transform(ts_curr_in)

    # Array to store the day-ahead 'predictions'.
    xout_norm = np.zeros([ts_curr_norm.shape[0],
                          ts_curr_norm.shape[1], n_samples])
    # pred_stds = np.zeros_like(ts_curr_norm)

    # Simply copy the input for the first l_start+l_step hours.
    xout_norm[0:l_start+l_step, :, :] = np.tile(np.atleast_3d(
        ts_curr_norm[0:l_start+l_step, :]), n_samples)

    # Month tracker index.
    d = 0

    # Loop over each day, starting from the second day.
    for h in range(l_start+l_step, l_end+l_step, l_step):

        # Train initial model on the first day.
        # Take the last $histlim$ hours into account.
        slicer = range(l_start, h)

        if len(slicer) > histlim:
            b = int(len(slicer) - histlim)
            slicer = slicer[b:h+1]
            del b

        train_slice = ts_curr_norm[slicer, :]

        # Index for the next day's data.
        tomorrows_slice = range(h, h+l_step)

        # Find the current month and date.
        # month_tracker[m] = np.squeeze(currentdata.month.values[h == hh])

        print("Sampling for month %d, day %d" % (month_tracker[d], int(h/24)))

        # The 'predicted' next day's values (size = l_step).
        tomorrows_vals = np.zeros((l_step, ts_curr_norm.shape[1],
                                   n_samples))

        for c, colname in enumerate(column_names):

            if colname in date_cols:
                # Use original values of time variables.
                tomorrows_vals[:, c, :] = np.tile(
                        np.atleast_2d(ts_curr_norm[tomorrows_slice, c]).T,
                        tomorrows_vals.shape[-1])

            else:

                # Find the current column in the overall array.
                find_y = yfinder(colname)

                # Separate the training data (history) into x & y.
                x_train = train_slice[:, np.logical_not(find_y)]
                y_train = train_slice[:, find_y]

                # The model will be 'queried' at these points.
                # The query points do not include the variable
                # which is the current output.
                x_query = x_train[-l_step:, :]

                x_pred = np.zeros([l_step, n_samples])

                month_gps = [gp_list[c][a] for a, b in
                             enumerate(month_tracker[d] ==
                                       month_tracker) if b]
                for s in range(0, x_pred.shape[1]):
                    # Select a random gp model from this month's models.
                    gp_sel = random.choice(month_gps)

                    # Take a sample from that model.
                    x_pred[:, s] = np.squeeze(gp_sel.sample_y(
                                x_query, n_samples=1))

                    # Use the GP of each day to predict that day,
                    # but add the predictions from using the GPs of
                    # surrounding days as well.

                # Store the values for 'tomorrow' for this variable.

                tomorrows_vals[:, c, :] = x_pred

        xout_norm[tomorrows_slice, :, :] = tomorrows_vals

        d += 1

        # End of hourly loop.

    # Un-scale the series.
    # Create a ndarray like xout_norm.
    xout_un = np.zeros_like(xout_norm)
    # Inverse transform all columns together, for each sample.
    for n in range(0, n_samples):
        xout_un[:, :, n] = s_c.inverse_transform(xout_norm[:, :, n])

    # Synthetic solar data requires post-processing.
    for c, colname in enumerate(column_names):

        if colname in ['ghi', 'dni', 'dhi']:

            for n in range(0, n_samples):
                # Using the source data - check to see if there
                # should be sunlight at a given hour. If not,
                # then set corresponding synthetic value to zero.

                temp = xout_un[:, c, n]
                temp[ts_curr_in[:, c] == 0] = 0
                xout_un[:, c, n] = temp

                # If there is a negative value (usually at sunrise
                # or sunset), interpolate.
                temp = xout_un[:, c, n]
                temp[temp < 0] = np.nan
                nans = np.isnan(temp)
                temp[nans] = np.interp(ts_curr_in[nans, 1],
                    ts_curr_in[~nans, 1], temp[~nans])
                xout_un[:, c, n] = temp

            # End loop over samples.

        # End colname if statement.

    # A potential improvement would be to calculate sunrise and sunset
    # independently since that is an almost deterministic calculation.

    pd.to_pickle(xout_un, picklepath)

    end_time = time.monotonic()
    print("Time taken to sample models was %6.2f seconds."
          % (end_time - start_time))

    return xout_un, column_names


def yfinder(colname):

    # Find the current column in the hourly values matrix.
    findy = np.zeros(len(column_names), dtype=bool)
    col = [i for i, x in enumerate(column_names) if x == colname]
    findy[col] = True
    del col

    # Tack on a vector of False for the daily means that will be
    # added to the input matrix below.
    # findy = np.hstack([findy, np.zeros(len(column_names)-2, dtype=bool)])
    # The right-half of this logical index is all False since the daily
    # means for the next day are always part of the input.
    # That is, they are NEVER selected as outputs.

    return findy