#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 16:23:50 2017

@author: Parag Rastogi
"""

import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
# , RationalQuadratic, ExpSineSquared, Matern


# This is the master tuple of column names, which should not be modified.
column_names = ('day', 'hour', 'tdb', 'tdp', 'rh', 'ghi', 'dni', 'dhi',
                'wspd', 'wdr')

# This is the master tuple of time variable names,
# which should also not be modified.
date_cols = ('day', 'hour')
dc = len(date_cols)


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


# Set a kernel with seasonal periodicity (ExpSineSquared) and allow it to
# decay away from exact periodicity (RBF).
# Allow for a rising/falling trend with a separate RBF kernel and noise
# with a WhiteKernel.


def gp_funcs(series, xtrain, ytrain, xtest=np.empty(0), *args):

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
