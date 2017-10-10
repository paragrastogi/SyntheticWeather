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
# from wfileio import get_weather
from ts_models import select_models
from petites import solarcleaner


def resampling(xy_train, selmdl=None, ffit=None, train=True,
               sample=True, n_sample=1,
               picklepath='./xxx.npy'):
    
    print(np.random.seed())

    # Temporarily here - to be eventually fed in from main script.
    # stcode = 'gen'
    # randseed = 8760
    # xy_train, locdata, header = get_weather(
    #        stcode, "./gen_iwec.epw", "epw", outpath=stcode)
    #    picklepath = './syn_gen_8760_res'
    ###

    #    # Seed random number generators.
    #    np.random.seed(randseed)
    #    random.seed = randseed

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

    numvars = 2  # For now, only considering two variables.

    if train:

        # The curve_fit function outputs two things:
        params = [
                curve_fit(fourier.fit_tdb, fit_idx, xy_train[:, 4]),
                curve_fit(fourier.fit_rh, fit_idx, xy_train[:, 6])
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
        demeaned = [xy_train[:, 4] - ffit[0], xy_train[:, 6] - ffit[1]]

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
        n_sample = 50  # This should be an input to the function.

        selmdl = list()
        selmdl_type = list()
        resid = np.zeros([demeaned[0].shape[0], numvars])

        for idx, ser in enumerate(demeaned):
            mdl_temp, type_temp, resid[:, idx] = select_models(
                arp, maq, sarp, smaq, s, ser)

            selmdl.append(mdl_temp)
            selmdl_type.append(type_temp)

    else:

        if selmdl is None:
            print("You did not ask me to train a model but didn't " +
                  "supply a valid model either. Terminating with " +
                  "None outputs.")
            return None, None, None

    # %%

    if sample:

        resampled = np.zeros([8760, numvars, n_sample])

        for v in range(0, numvars):
            for n in range(0, n_sample):
                resampled[:, v, n] = selmdl[v].simulate(
                        nsimulations=8760)

        # %%

        # Add the resampled time series back to the fourier series.
        ts_syn = np.zeros_like(resampled)

        for v in range(0, numvars):
            ts_syn[:, v, :] = resampled[:, v, :] + \
                np.resize(ffit[v], resampled[:, v, :].shape)

        ###
        ### Add solarcleaner here. ###
        ###

        # Save the outputs as a pickle.
        np.save(picklepath, ts_syn, allow_pickle=True)

    else:
        ts_syn = None

    return ffit, selmdl, ts_syn
