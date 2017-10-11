#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 11:47:58 2017

@author: rasto

Translating, as faithfully as possible, the resampling method first proposed
my thesis (Rastogi, 2016, EPFL).
"""

# import os
# import random

import numpy as np
# import pandas
from scipy.optimize import curve_fit

# from sklearn.preprocessing import StandardScaler

# from matplotlib import pyplot as plt
# from IPython import get_ipython
# get_ipython().run_line_magic('matplotlib', 'inline')

import fourier
from ts_models import select_models
# from petites import setseed
from petites import solarcleaner

# This is the master tuple of column names, which should not be modified.
column_names = ('year', 'month', 'day', 'hour', 'tdb', 'tdp', 'rh',
                'ghi', 'dni', 'dhi', 'wspd', 'wdr')

# This is the master tuple of time variable names,
# which should also not be modified.
date_cols = ('year', 'month', 'day', 'hour')
dc = len(date_cols)
midx = 1  # Month is in the second column - will be needed later.


def resampling(xy_train, counter=0, selmdl=None, ffit=None, train=True,
               sample=True, n_samples=1,
               picklepath='./xxx.npy', randseed=None):

    # Temporarily here - to be eventually fed in from main script.
    # stcode = 'gen'
    # randseed = 8760
    # xy_train, locdata, header = get_weather(
    #        stcode, "./gen_iwec.epw", "epw", outpath=stcode)
    #    picklepath = './syn_gen_8760_res'
    ###

    # Check to see if random number generation is reproducible.

    # Seed random number generators.
    #    if randseed is None:
    #        import time
    #        randseed = int(time.time())
    #    setseed(randseed)

    #    # This is the master tuple of column names, which should
    #    # not be modified.
    #    column_names = ('year', 'month', 'day', 'hour', 'tdb', 'tdp', 'rh',
    #                    'ghi', 'dni', 'dhi', 'wspd', 'wdr')
    #
    #    # This is the master tuple of time variable names,
    #    # which should also not be modified.
    #    date_cols = ('year', 'month', 'day', 'hour')
    #    dc = len(date_cols)
    #    midx = 1  # Month is in the second column - will be needed later.

    # # Make a scaler - either standard (with \mu = 0 and
    # # \sigma = 1) or robust (with median and iqr).
    # scaler = StandardScaler()

    fit_idx = np.arange(0, xy_train.shape[0])

    # %%

    # Fit fourier functions to the tdb and rh series.

    varidx = [4, 6]  # For now, only considering two variables.
    othervars = (np.arange(0, xy_train.shape[1])).tolist()

    if train:

        # The curve_fit function outputs two things:
        params = [
                curve_fit(fourier.fit_tdb, fit_idx, xy_train[:, varidx[0]]),
                curve_fit(fourier.fit_rh, fit_idx, xy_train[:, varidx[1]])
                ]
        # tdb, rh, tdb_low, tdb_high, rh_low

        # Call the fourier fit function with the calculated
        # parameters to get the
        # values of the fourier fit at each time step
        ffit = [fourier.fit('tdb', fit_idx, *params[0][0]),
                fourier.fit('rh', fit_idx, *params[1][0])]

        # Now subtract the low- and high-frequency fourier fits
        # (whichever is applicable) from the raw values to get the
        # 'de-meaned' values (values from which the mean has
        # been removed).
        demeaned = [xy_train[:, varidx[0]] - ffit[0],
                    xy_train[:, varidx[1]] - ffit[1]]

        # %%

        # Fit ARIMA models.

        # Set ranges for various model parameters.
        # Each range is one more than what we are
        # interested in because range cuts off at end-1.

        arp = range(0, 2)
        maq = range(0, 2)
        sarp = range(0, 2)
        smaq = range(0, 2)
        s = 24
        # n_samples = 50  # This should be an input to the function.

        selmdl = list()
        selmdl_type = list()
        resid = np.zeros([demeaned[0].shape[0], len(varidx)])

        for idx, ser in enumerate(demeaned):
            mdl_temp, type_temp, resid[:, idx] = select_models(
                arp, maq, sarp, smaq, s, ser)

            selmdl.append(mdl_temp)
            selmdl_type.append(type_temp)

    # %%

        resampled = np.zeros([8760, len(varidx), n_samples])

        for v in range(0, len(varidx)):
            for n in range(0, n_samples):
                resampled[:, v, n] = selmdl[v].simulate(
                        nsimulations=8760)

        # %%

        # Add the resampled time series back to the fourier series.
        xout = np.zeros([resampled.shape[0], xy_train.shape[1],
                          resampled.shape[-1]])

        v = 0
        vv = 0
        for idx in range(0, xy_train.shape[1]):
            if idx in varidx:
                xout[:, idx, :] = resampled[:, v, :] + \
                    np.resize(ffit[v], resampled[:, v, :].shape)
                v += 1
            elif idx in othervars:
                xout[:, idx, :] = np.resize(
                        xy_train[:, idx], xout[:, idx, :].shape)
                vv += 1

        # Synthetic solar data requires post-processing.
        for c, colname in enumerate(column_names):

            if colname in ['ghi', 'dni', 'dhi']:

                for n in range(0, n_samples):
                    xout[:, c, n] = solarcleaner(
                            xout[:, c, n], xy_train[:, c],
                            xy_train[:, 3])

                # End loop over samples.

            # End colname if statement.

        # End colname for loop.


        # Save the outputs as a pickle.
        np.save(picklepath, xout, allow_pickle=True)

        sample = None

    else:

        try:

            xout = np.load(picklepath)

            sample = xout[:, :, counter]

        except AttributeError:
            print("I could not open the pickle file with samples. " +
                  "Please check it exists at {0}.".format(picklepath))

#        if selmdl is None:
#            print("You did not ask me to train a model but didn't " +
#                  "supply a valid model either. Terminating with " +
#                  "None outputs.")
#            return None, None, None

    return ffit, selmdl, sample
