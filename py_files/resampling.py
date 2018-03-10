#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 11:47:58 2017

@author: rasto

Translating, as faithfully as possible, the resampling method first proposed
my thesis (Rastogi, 2016, EPFL).
"""

import numpy as np
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
NUM_NBOURS = 10

# This is the master tuple of column names, which should not
# be modified.
COLUMNS = ('year', 'month', 'day', 'hour', 'tdb', 'tdp', 'rh',
           'ghi', 'dni', 'dhi', 'wspd', 'wdr')


def resampling(xy_train, train=True, n_samples=10,
               picklepath='test_out.npy', counter=0,
               arma_params=None, bounds=None):

    # Reassign defaults if incoming list params are None
    # (i.e., nothing passed.)
    if arma_params is None:
        arma_params = [2, 2, 1, 1, 24]

    if bounds is None:
        bounds = [0.01, 99.9]

    # if train is True, then the model must be fit to incoming data.
    if train:

        ffit, selmdl, xout = trainer(
            xy_train, n_samples, arma_params, bounds)

        # Save the outputs as a pickle.
        np.save(picklepath, xout, allow_pickle=True)

        sample = None

    else:

        # The non-training call is only meant to return a sample.
        ffit = None
        selmdl = None

        sample = sampler(picklepath, counter)

    # Return all three outputs.

    return ffit, selmdl, sample


def trainer(xy_train, n_samples, arma_params, bounds):
    '''Train the model with this function.'''

    fit_idx = np.arange(0, xy_train.shape[0])

    # For now, I am only considering these two variables.
    tdb_idx = [x for (x, y) in enumerate(COLUMNS) if y == "tdb"][0]
    rh_idx = [x for (x, y) in enumerate(COLUMNS) if y == "rh"][0]

    sol_idx = [x for (x, y) in enumerate(COLUMNS) if y == "ghi"] + \
        [x for (x, y) in enumerate(COLUMNS) if y == "dni"] + \
        [x for (x, y) in enumerate(COLUMNS) if y == "dhi"]

    # Later, we will calculate tdp values from tdb and rh, so store the index.
    tdp_idx = [x for (x, y) in enumerate(COLUMNS) if y == "tdp"]

    # othervars = (np.arange(0, xy_train.shape[1]-2)).tolist()

    # Fit fourier functions to the tdb and rh series.

    # The curve_fit function outputs two things: parameters of the fit and
    # the estimated covariance. We only use the first.
    # Inputs are the function to fit (fourier in this case),
    # xdata, and ydata.
    params = [
        curve_fit(fourier.fit_tdb, fit_idx, xy_train[:, tdb_idx]),
        curve_fit(fourier.fit_rh, fit_idx, xy_train[:, rh_idx])
        ]

    # Call the fourier fit function with the calculated
    # parameters to get the values of the fourier fit at each time step
    ffit = [fourier.fit('tdb', fit_idx, *params[0][0]),
            fourier.fit('rh', fit_idx, *params[1][0])]

    # Now subtract the low- and high-frequency fourier fits
    # (whichever is applicable) from the raw values to get the
    # 'de-meaned' values (values from which the mean has
    # been removed).
    DeMeaned = [x-y for x, y in zip([xy_train[:, tdb_idx],
                                     xy_train[:, rh_idx]], ffit)]

    # Fit ARIMA models.

    selmdl = list()
    resid = np.zeros([DeMeaned[0].shape[0], NUM_VARS])

    for idx, ser in enumerate(DeMeaned):
        mdl_temp, resid[:, idx] = select_models(
            arma_params, ser)
        selmdl.append(mdl_temp)

    print("Done with fitting models to TDB and RH.\r\n")

    resampled = np.zeros([STD_LEN_OUT, NUM_VARS, n_samples])

    print(("Simulating the learnt model to get synthetic noise series."
           "This might take some time.\r\n"))
    for mdl_idx, mdl in enumerate(selmdl):
        for sample_num in range(0, n_samples):
            resampled[:, mdl_idx, sample_num] = mdl.simulate(
                nsimulations=STD_LEN_OUT)
        # End n for loop.
    # End mdl for loop.

    # Add the resampled time series back to the fourier series.

    # First make the xout array using all variables. Variables other than RH
    # and TDB are just repeated from the incoming files.
    xout = np.tile(np.atleast_3d(xy_train), n_samples)

    for train_idx, xout_idx in enumerate([rh_idx, tdb_idx]):

        # Add the fourier fits from the training data to the
        # resampled/resimulated ARMA model outputs.
        xout[:, xout_idx, :] = (
            resampled[:, train_idx, :] +
            np.reshape(np.tile(ffit[train_idx], n_samples),
                       (STD_LEN_OUT, -1))
            )

    # End for rh_idx, tdp_idx loop.

    # Clean the tdb and rh values.

    print('I am going to clean all synthetic values.\r\n')

    # All synthetic data requires post-processing.

    # Clean the generated temperature values using extreme percentiles
    # as proxies for 'outliers'.

    # Clean the RH values using phyiscal limits (0-100).

    for sample_num in range(0, n_samples):

        xout[:, tdb_idx, sample_num] = petite.quantilecleaner(
            xout[:, tdb_idx, sample_num], xy_train, bounds=bounds)

        xout[:, rh_idx, sample_num] = petite.quantilecleaner(
            xout[:, rh_idx, sample_num], xy_train, bounds=bounds)

        xout[:, tdp_idx, sample_num] = np.resize(
            petite.calc_tdp(xout[:, tdb_idx, sample_num],
                            xout[:, rh_idx, sample_num]),
            xout[:, tdp_idx, sample_num].shape)

        # End loop over samples.

    # Calculate daily means of temperature.
    tdb_syn_dailymeans = np.mean(np.reshape(np.squeeze(xout[:, tdb_idx, :]),
                                            [-1, 24, n_samples]), axis=1)
    tdb_dailymeans = np.mean(np.reshape(np.squeeze(xout[:, tdb_idx]),
                                        [-1, 24]), axis=1)

    # Cycle through each array of daily means.
    for idx, syn_means in enumerate(tdb_syn_dailymeans.T):
        # Find the nearest neighbour to each element of the syn_means array.
        # Nearest neighbours must be in the same month to preserve length
        # of day.

        for month in range(0, 12):

            # THIS ISN'T WORKING.
            this_month = tdb_syn_dailymeans[:, 0] == month

            # First find 10 nearest neighbours.
            nearest_nbour_temp = [
                (np.abs(x-tdb_dailymeans[this_month])).argsort(
                )[:NUM_NBOURS] for x in syn_means[this_month]]
            # This works by taking the indices from sorting the values,
            # in ascending order, and then picking the first four of those.

            # Now take one random nearest neighbour.
            nearest_nbour = [x[np.random.randint(0, NUM_NBOURS)]
                             for x in nearest_nbour_temp]

            print(len(nearest_nbour))
            print(len(nearest_nbour[0]))
            print(np.sum(this_month))
    
    # for train_idx, xout_idx in enumerate(sol_idx):

    # # Add the fourier fits from the training data to the
    # # resampled/resimulated ARMA model outputs.
    # xout[:, xout_idx, :] = (
    #     resampled[:, train_idx, :] +
    #     np.reshape(np.tile(ffit[train_idx], n_samples),
    #                (STD_LEN_OUT, -1))
    #     )

    for sample_num in range(0, n_samples):

        for solcols in sol_idx:
            xout[:, solcols, sample_num] = petite.solarcleaner(
                xout[:, solcols, sample_num], xy_train[:, solcols])

        # End sol_idx loop

    # End for sample_num loop.

    return ffit, selmdl, xout


def sampler(picklepath, counter):
    '''Only opens the pickle of saved samples and returns ONE sample.'''

    try:

        xout = np.load(picklepath)

        sample = xout[:, :, counter]

    except AttributeError:

        print("I could not open the pickle file with samples. " +
              "Please check it exists at {0}.".format(picklepath))
        sample = None

    return sample
