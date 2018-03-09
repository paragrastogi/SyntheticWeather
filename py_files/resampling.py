#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 11:47:58 2017

@author: rasto

Translating, as faithfully as possible, the resampling method first proposed
my thesis (Rastogi, 2016, EPFL).
"""

import numpy as np


def resampling(xy_train, train=True, n_samples=10,
               picklepath='test_out.npy', counter=0,
               arma_params=[2, 2, 2, 2, 24], lb=0.01, ub=99.9):

    # Number of variables affected by this script.
    # Only two for now.
    NUM_VARS = 2

    # %%

    # if train is True, then the model must be fit to incoming data.
    if train:

        ffit, selmdl, xout = trainer(
            xy_train, n_samples, arma_params, NUM_VARS)

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


def trainer(xy_train, n_samples, arma_params, NUM_VARS, lb, ub):
    '''Train the model with this function.'''

    from scipy.optimize import curve_fit

    import fourier
    from ts_models import select_models
    from petites import solarcleaner
    from petites import rhcleaner
    from petites import calc_tdp
    from petites import quantilecleaner

    # Progress bar!
    # from tqdm import tqdm

    # "Standard" length of output year.
    STD_LEN_OUT = 8760

    # This is the master tuple of column names, which should not
    # be modified.
    column_names = ('year', 'month', 'day', 'hour', 'tdb', 'tdp', 'rh',
                    'ghi', 'dni', 'dhi', 'wspd', 'wdr')

    # This is the master tuple of time variable names,
    # which should also not be modified.
    # date_cols = ('year', 'month', 'day', 'hour')
    # dc = len(date_cols)
    # midx = 1  # Month is in the second column - will be needed later.

    # For now, I am only considering these two variables.
    fit_idx = np.arange(0, xy_train.shape[0])
    tdb_idx = [x for (x, y) in enumerate(column_names) if y == "tdb"][0]
    rh_idx = [x for (x, y) in enumerate(column_names) if y == "rh"][0]
    sol_idx = [x for (x, y) in enumerate(column_names) if y == "ghi"] + \
        [x for (x, y) in enumerate(column_names) if y == "dni"] + \
        [x for (x, y) in enumerate(column_names) if y == "dhi"]

    # Later, we will calculate tdp values from tdb and rh, so store the index.
    tdp_idx = [x for (x, y) in enumerate(column_names) if y == "tdp"]

    # varidx = [x for (x, y) in enumerate(column_names) if y == "tdb"] + \
    #     [x for (x, y) in enumerate(column_names) if y == "rh"]
    othervars = (np.arange(0, xy_train.shape[1])).tolist()

    # Fit fourier functions to the tdb and rh series.

    # The curve_fit function outputs two things:
    params = [
        curve_fit(fourier.fit_tdb, fit_idx, xy_train[:, tdb_idx]),
        curve_fit(fourier.fit_rh, fit_idx, xy_train[:, rh_idx])
        ]

    # Call the fourier fit function with the calculated
    # parameters to get the
    # values of the fourier fit at each time step
    ffit = [fourier.fit('tdb', fit_idx, *params[0][0]),
            fourier.fit('rh', fit_idx, *params[1][0])]

    # Now subtract the low- and high-frequency fourier fits
    # (whichever is applicable) from the raw values to get the
    # 'de-meaned' values (values from which the mean has
    # been removed).
    demeaned = [xy_train[:, tdb_idx] - ffit[0],
                xy_train[:, rh_idx] - ffit[1]]

    # Fit ARIMA models.

    selmdl = list()
    resid = np.zeros([demeaned[0].shape[0], NUM_VARS])

    for idx, ser in enumerate(demeaned):
        mdl_temp, resid[:, idx] = select_models(
            arma_params, ser)
        selmdl.append(mdl_temp)

    print("Done with fitting models to TDB and RH.")

    resampled = np.zeros([STD_LEN_OUT, NUM_VARS, n_samples])

    print("Simulating the learnt model to get synthetic noise series.")
    print("This might take some time.\r\n")
    for mdl_idx, mdl in enumerate(selmdl):
        for n in range(0, n_samples):
            resampled[:, mdl_idx, n] = mdl.simulate(
                nsimulations=STD_LEN_OUT)
        # End n for loop.
    # End mdl for loop.

    # Add the resampled time series back to the fourier series.
    xout = np.zeros([STD_LEN_OUT, xy_train.shape[1],
                     resampled.shape[-1]])

    v = 0
    vv = 0
    for idx in range(0, xy_train.shape[1]):
        if idx in [rh_idx, tdb_idx]:
            xout[:, idx, :] = (resampled[:, v, :] +
                               np.repeat(np.reshape(ffit[v], [-1, 1]),
                               resampled.shape[-1], axis=1))
            v += 1

        elif idx in othervars:
            xout[:, idx, :] = np.repeat(np.reshape(
                xy_train[:, idx], [-1, 1]),
                resampled.shape[-1], axis=1)
            vv += 1

        # End if idx conditional.
    # End for idx loop.

    # Clean the tdb and rh values.

    print('I am going to clean all synthetic values.')

    # All synthetic data requires post-processing.

    # Clean the generated temperature values using extreme percentiles
    # as proxies for 'outliers'.

    # Clean the RH values using phyiscal limits (0-100).

    for n in range(0, n_samples):

        for solcols in sol_idx:
            xout[:, solcols, n] = solarcleaner(
                xout[:, solcols, n], xy_train[:, solcols])

        xout[:, tdb_idx, n] = quantilecleaner(xout[:, tdb_idx, n],
            lb=lb, ub=ub)

        xout[:, rh_idx, n] = rhcleaner(xout[:, rh_idx, n])

        xout[:, tdp_idx, n] = np.resize(
            calc_tdp(xout[:, tdb_idx, n], xout[:, rh_idx, n]),
            xout[:, tdp_idx, n].shape)

        # End loop over samples.

    # End colname for loop.

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
