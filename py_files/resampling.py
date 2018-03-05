#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 11:47:58 2017

@author: rasto

Translating, as faithfully as possible, the resampling method first proposed
my thesis (Rastogi, 2016, EPFL).
"""

import numpy as np


def resampling(xy_train, counter=0, selmdl=None, ffit=None,
               train=True, sample=True, n_samples=1,
               picklepath='test_out.npy', randseed=None,
               arp_ub=2, maq_ub=2, sarp_ub=2, smaq_ub=2, s=24):

    # Set ranges for various model parameters.
    # Each range is one more than what we are
    # interested in because range cuts off at end-1.
    arp = range(0, arp_ub+1)
    maq = range(0, arp_ub+1)
    sarp = range(0, arp_ub+1)
    smaq = range(0, arp_ub+1)

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

    # %%

    # if train is True, then the model must be fit to incoming data.
    if train:

        from scipy.optimize import curve_fit

        import fourier
        from ts_models import select_models
        from petites import solarcleaner
        from petites import rhcleaner
        from petites import calc_tdp

        # Progress bar!
        # from tqdm import tqdm

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
        varidx = [x for (x, y) in enumerate(column_names) if y == "tdb"] + \
            [x for (x, y) in enumerate(column_names) if y == "rh"]
        othervars = (np.arange(0, xy_train.shape[1])).tolist()

        # Fit fourier functions to the tdb and rh series.

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

        selmdl = list()
        resid = np.zeros([demeaned[0].shape[0], len(varidx)])

        for idx, ser in enumerate(demeaned):
            mdl_temp, resid[:, idx] = select_models(
                arp, maq, sarp, smaq, s, ser)
            selmdl.append(mdl_temp)

        print("Done with fitting models to TDB and RH.")
        # %%

        resampled = np.zeros([8760, len(varidx), n_samples])

        print("Simulating the learnt model to get synthetic noise series.")
        print("This might take some time.\r\n")
        for v, mdl in enumerate(selmdl):
            for n in range(0, n_samples):
                resampled[:, v, n] = mdl.simulate(nsimulations=8760)

        # %%

        # Add the resampled time series back to the fourier series.
        xout = np.zeros([resampled.shape[0], xy_train.shape[1],
                         resampled.shape[-1]])

        v = 0
        vv = 0
        for idx in range(0, xy_train.shape[1]):
            if idx in varidx:
                xout[:, idx, :] = resampled[:, v, :] + \
                    np.repeat(np.reshape(ffit[v], [-1, 1]),
                              resampled.shape[-1], axis=1)
                v += 1

            elif idx in othervars:
                xout[:, idx, :] = np.repeat(np.reshape(
                        xy_train[:, idx], [-1, 1]),
                        resampled.shape[-1], axis=1)
                vv += 1

            # End if idx conditional.
        # End for idx loop.

        # Synthetic solar and rh data require post-processing.

        for c, colname in enumerate(column_names):

            if colname in ['ghi', 'dni', 'dhi']:

                for n in range(0, n_samples):
                    xout[:, c, n] = solarcleaner(
                            xout[:, c, n], xy_train[:, c])

                # End loop over samples.

            # End colname if statement.

        # End colname for loop.

        # Interpolate the bad rh values.
        for n in range(0, n_samples):
            rh = xout[:, varidx[1], n]
            xout[:, varidx[1], n] = rhcleaner(rh)

        # Calculate tdp values from tdb and rh.
        tdp_idx = [x for (x, y) in enumerate(column_names) if y == "tdp"]

        for n in range(0, n_samples):
            xout[:, tdp_idx, n] = np.resize(
                    calc_tdp(xout[:, varidx[0], n], xout[:, varidx[1], n]),
                    xout[:, tdp_idx, n].shape)

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
            sample = None

        # The non-training call is only meant to return a sample.
        ffit = None
        selmdl = None

    return ffit, selmdl, sample
