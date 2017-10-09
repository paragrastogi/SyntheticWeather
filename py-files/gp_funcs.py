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

from petites import solarcleaner

# This is the master tuple of column names, which should not be modified.
column_names = ('year', 'month', 'day', 'hour', 'tdb', 'tdp', 'rh',
                'ghi', 'dni', 'dhi', 'wspd', 'wdr')

# This is the master tuple of time variable names,
# which should also not be modified.
date_cols = ('year', 'month', 'day', 'hour')
dc = len(date_cols)
midx = 1  # Month is in the second column - will be needed later.

# Make a scaler - either standard (with \mu = 0 and
# \sigma = 1) or robust (with median and iqr).
scaler = StandardScaler()

# Alternate (more complicated) kernels:
# With seasonal periodicity (ExpSineSquared) and allow it to
# decay away from exact periodicity (RBF).
# Allow for a rising/falling trend with a separate RBF kernel and noise
# with a WhiteKernel.


def fitgp(series, xtrain, ytrain, xtest=np.empty(0), gp_in=None):

    # Initial values for hyperparameters.
    hypguess = np.random.random(xtrain.shape[1])

    # n_restarts_optimizer selects "result from run with minimal
    # (negative) log-marginal likelihood", line 232-236 in
    # github.com/scikit-learn/scikit-learn/blob/master/sklearn/
    # gaussian_process/gpr.py#L206 .
    n_restarts = 5

    # WhiteKernel represents signal noise ($\sigma_n$).
    if series in ['tdb_dm', 'ghi_dm', 'tdp_dm']:
        # Monthly data.
        kd = RBF(length_scale=hypguess) + \
            WhiteKernel(noise_level=1, noise_level_bounds=[0.9, 1.1])
    elif series in column_names:
        # Hourly data.
        kd = RBF(length_scale=hypguess) + \
            WhiteKernel(noise_level=0.1, noise_level_bounds=[1e-2, 0.1])

    if gp_in is None:
        # No incoming regressor object, so learn hyper-parameters.
        gp_out = GaussianProcessRegressor(
                kernel=kd, alpha=0, n_restarts_optimizer=n_restarts,
                optimizer='fmin_l_bfgs_b', normalize_y=False)

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

# ----------- END fitgp function. -----------


def learngp(l_start, l_end, l_step, histlim,
            masterdata, path_fldr_pickle=os.getcwd()):

    start_time = time.monotonic()

    # Re-scale the hourly values (Normalize).
    s_c = scaler.fit(masterdata)
    # s_h = scaler.fit(ts_hist_in)
    ts_norm = s_c.transform(masterdata)

    # Number of variables being considered =
    # number of columns minus time variables.
    metvars = ts_norm.shape[1]

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

    # Hours of the year - calculated from hours of the day and
    # day of the year.
    hh = masterdata[:, 3] + (masterdata[:, 2] - 1)*24

    d = 0

    mtrack = np.zeros(len(range(l_start+l_step, l_end+l_step, l_step)))

    # Loop over each day, starting from the second day.
    for h in range(l_start+l_step, l_end+l_step, l_step):

        # Train initial model on the first day.
        # Take the last $histlim$ hours into account.
        xslice = np.arange(l_start, h)

        if len(xslice) > histlim:
            b = int(len(xslice) - histlim)
            xslice = xslice[b:h+1]
            del b

        # The time slice from which the function will pick up the
        # 'y' it needs to 'predict'.
        yslice = np.arange(h, h+xslice.shape[0])

        # This bit of code is dormant right now but should be revived when
        # the function is updated to learn from more than one year of data.
        # Then the train slice would be selected by month, day, and hour,
        # with the same timestamps used from each year.
        #        date_slice = [(mslice, dslice, hslice)
        #                      for (mslice, dslice, hslice) in
        #                      zip(masterdata[:, 0], masterdata[:, 1], hh)
        #                      if hslice in xslice]

        # If the train or predict slice is beyond the 31st of December,
        # then break this loop.
        if (h+xslice.shape[0]) > ts_norm.shape[0]:
            break

        # Output the current month and day.
        mtrack[d] = int(np.unique(np.squeeze(
                masterdata[h == hh, midx])))
        # Change the unique function here to accomodate
        # multi-year records.

        print("Model no. {0}, month {1}, day {2}".format(
                d, mtrack[d], int(h/24)))

        # Concatenate the daily means with the hourly values, leaving
        # out the day variable from the daily means array (first column).
        train_slice = ts_norm[xslice, :]
        train_y = ts_norm[yslice, :]
        # USE dateslice here instead of xslice.

        for c, colname in enumerate(column_names):

            if colname in date_cols:
                # Use original values of time variables.
                # tomorrows_vals[:, c, :] = np.tile(np.atleast_2d(
                #     ts_norm[tomorrows_slice, c]).T, n_samples)
                continue

            else:
                # Find the current column in the overall array.
                find_y = yfinder(colname)

                # Separate the training data (history) into x & y.
                x_train = train_slice
                y_train = train_y[:, find_y]

                # Remove year from training data (first column).
                x_train = x_train[:, 1:]

                # Call the fitgp function.
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

    return gp_list, mtrack, s_c

# ----------- END learngp function. -----------


def samplegp(gp_list, l_start, l_end, l_step, histlim, n_samples,
             masterdata, mtrack, s_c, picklepath="./xxx.npy"):

    start_time = time.monotonic()

    # Re-scale the hourly values (Normalize).
    # s_c = scaler.fit(masterdata)
    # s_h = scaler.fit(ts_hist_in)
    ts_norm = s_c.transform(masterdata)

    # Array to store the day-ahead 'predictions'.
    xout_norm = np.zeros([ts_norm.shape[0],
                          ts_norm.shape[1], n_samples])

    # Simply copy the input for the first l_start+l_step hours.
    xout_norm[0:l_start+l_step, :, :] = np.tile(np.atleast_3d(
        ts_norm[0:l_start+l_step, :]), n_samples)

    # Month tracker index.
    d = 0

    # Loop over each day, starting from the second day.
    for h in range(l_start+l_step, l_end+l_step, l_step):

        # Train initial model on the first day.
        # Take the last $histlim$ hours into account.
        xslice = range(l_start, h)

        if len(xslice) > histlim:
            b = int(len(xslice) - histlim)
            xslice = xslice[b:h+1]
            del b

        pred_slice = ts_norm[xslice, :]

        if (h+pred_slice.shape[0]) > ts_norm.shape[0]:
            break

        # Index for the next day's data.
        tomorrows_slice = range(h, h+pred_slice.shape[0])

        # Find the current month and date.
        # mtrack[m] = np.squeeze(currentdata.month.values[h == hh])

        print("Sampling for month %d, day %d" % (mtrack[d], int(h/24)))

        # The 'predicted' next day's values (size = l_step).
        tomorrows_vals = np.zeros((pred_slice.shape[0], ts_norm.shape[1],
                                   n_samples))

        for c, colname in enumerate(column_names):

            if colname in date_cols:
                # Use original values of time variables.
                tomorrows_vals[:, c, :] = np.tile(
                        np.atleast_2d(ts_norm[tomorrows_slice, c]).T,
                        tomorrows_vals.shape[-1])

            else:

                # Find the current column in the overall array.
                # find_y = yfinder(colname)

                # Separate the training data (history) into x & y.
                # x_train = pred_slice

                # Remove year from training data (first column).
                x_query = pred_slice[:, 1:]

                # The model will be 'queried' at these points.
                # The query points do not include the variable
                # which is the current output.
                # x_query = x_train[-l_step:, :]

                # Pre-allocate the array that will contain predicted/sampled y.
                y_pred = np.zeros([x_query.shape[0], n_samples])

                month_gps = [gp_list[c][a] for a, b in
                             enumerate(mtrack[d] ==
                                       mtrack) if b]
                for s in range(0, y_pred.shape[1]):
                    # Select a random gp model from this month's models.
                    gp_sel = random.choice(month_gps)

                    # Take a sample from that model.
                    y_pred[:, s] = np.squeeze(gp_sel.sample_y(
                                x_query, n_samples=1))

                    # Use the GP of each day to predict that day,
                    # but add the predictions from using the GPs of
                    # surrounding days as well.

                # Store the values for 'tomorrow' for this variable.

                tomorrows_vals[:, c, :] = y_pred

        xout_norm[tomorrows_slice, :, :] = tomorrows_vals

        d += 1

        # End of hourly loop.

    # Un-scale the series.
    # Create a ndarray like xout_norm.
    xout = np.zeros_like(xout_norm)
    # Inverse transform all columns together, for each sample.
    for n in range(0, n_samples):
        xout[:, :, n] = s_c.inverse_transform(xout_norm[:, :, n])

    # Synthetic solar data requires post-processing.
    for c, colname in enumerate(column_names):

        if colname in ['ghi', 'dni', 'dhi']:

            for n in range(0, n_samples):
                xout[:, c, n] = solarcleaner(
                        xout[:, c, n], masterdata[:, c], masterdata[:, 3])

            # End loop over samples.

        # End colname if statement.

    # End colname for loop.


    end_time = time.monotonic()
    print("Time taken to sample models was %6.2f seconds."
          % (end_time - start_time))

    # Save the outputs as a pickle.
    np.save(picklepath, xout, allow_pickle=True)

    return xout

# ----------- END samplegp function. -----------


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

# ----------- END yfinder function. -----------
