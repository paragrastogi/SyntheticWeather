#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 15:30:34 2017

@author: parag rastogi
"""

import random

import numpy as np
from scipy import interpolate
import pandas as pd


def setseed(randseed):
    '''Seed random number generators. Called as a function in main indra
    script once and only once.'''

    np.random.seed(randseed)
    random.seed = randseed

# ----------- END setseed function. -----------


def quantilecleaner(datain, xy_train, bounds=None):
    '''Generic cleaner based on quantiles. Needs a time series / dataset
       and cut-off quantiles. This function will censor the data outside
       those quantiles and interpolate the missing values using linear
       interpolation.'''

    if bounds is None:
        bounds = [0.01, 99.9]

    datain_quantiles = np.percentile(xy_train, bounds)

    dataout = pd.DataFrame(datain)

    dataout = dataout.mask(
        np.logical_or(dataout < datain_quantiles[0],
                      dataout > datain_quantiles[1]),
        other=np.NaN).interpolate(method='linear')

    # Pass back values with only one dimension.
    return np.squeeze(dataout.values)

# ----------- END quantilecleaner function. -----------


def solarcleaner(datain, master):

    '''Clean solar values by setting zeros at corresponding times in master
       to zero in the synthetic data. This is a proxy for sunrise, sunset,
       and twilight.'''

    # Using the source data - check to see if there
    # should be sunlight at a given hour. If not,
    # then set corresponding synthetic value to zero.

    datain[master == 0] = 0

    # If there is a negative value (usually at sunrise
    # or sunset), set it to zero as well.

    datain[datain < 0] = 0

    return datain

    # A potential improvement would be to calculate sunrise and sunset
    # independently since that is an almost deterministic calculation.

# ----------- END solarcleaner function. -----------


def rhcleaner(rh):

    '''RH values cannot be more than 100 or less than 0.'''

    rhout = pd.DataFrame(rh)

    rhout = rhout.mask(rhout > 100, other=99).mask(rhout < 0, other=1)

    return np.squeeze(rhout.values)

# ----------- END rhcleaner function. -----------


def wstats(datain, key, stat):

    grouped_data = datain.groupby(key)

    if stat is 'mean':
        dataout = grouped_data.mean()
    elif stat is 'sum':
        dataout = grouped_data.sum()
    elif stat is 'max':
        dataout = grouped_data.max()
    elif stat is 'min':
        dataout = grouped_data.min()
    elif stat is 'std':
        dataout = grouped_data.std()
    elif stat is 'q1':
        dataout = grouped_data.quantile(0.25)
    elif stat is 'q3':
        dataout = grouped_data.quantile(0.75)
    elif stat is 'med':
        dataout = grouped_data.median()

    return dataout

# ----------- END wstats function. -----------


def calc_tdp(tdb, rh):

    '''Calculate dew point temperature using dry bulb temperature
       and relative humidity.'''

    # Change relative humidity to fraction.
    phi = rh/100

    # Remove weird values.
    phi[phi > 1] = 1
    phi[phi < 0] = 0

    # Convert tdb to Kelvin.
    tdb_k = tdb + 273.15

    # Equations for calculating the saturation pressure
    # of water vapour, taken from ASHRAE Fundamentals 2009.
    # (Eq. 5 and 6, Psychrometrics)

    # Constants for Eq. 5, Temperature -200°C to 0°C.
    FROZEN_CONST = [-5.6745359*10**3, 6.3925247, -9.6778430*10**-3,
                    6.2215701*10**-7, 2.0747825*10**-9,
                    -9.4840240*10**-13, 4.1635019]

    # Constants for Eq. 6, Temperature 0°C to 200°C.
    LIQUID_CONST = [-5.8002206*10**3, 1.3914993, -4.8640239*10**-2,
                    4.1764768*10**-5, -1.4452093*10**-8, 6.5459673]

    # This is to distinguish between the two versions of equation 5.
    ice = tdb_k <= 273.15
    not_ice = np.logical_not(ice)

    lnp_ws = np.zeros(tdb_k.shape)

    # Eq. 5, pg 1.2
    lnp_ws[ice] = (
        FROZEN_CONST[0]/tdb_k[ice] + FROZEN_CONST[1] +
        FROZEN_CONST[2]*tdb_k[ice] + FROZEN_CONST[3]*tdb_k[ice]**2 +
        FROZEN_CONST[4]*tdb_k[ice]**3 + FROZEN_CONST[5]*tdb_k[ice]**4 +
        FROZEN_CONST[6]*np.log(tdb_k[ice]))

    # Eq. 6, pg 1.2
    lnp_ws[np.logical_not(ice)] = (
        LIQUID_CONST[0]/tdb_k[not_ice] + LIQUID_CONST[1] +
        LIQUID_CONST[2]*tdb_k[not_ice] +
        LIQUID_CONST[3]*tdb_k[not_ice]**2 +
        LIQUID_CONST[4]*tdb_k[not_ice]**3 +
        LIQUID_CONST[5]*np.log(tdb_k[not_ice]))

    # Temperature in the above formulae must be absolute,
    # i.e. in Kelvin

    # Continuing from eqs. 5 and 6
    p_ws = np.e**(lnp_ws)  # [Pa]

    # Eq. 24, pg 1.8
    p_w = (phi * p_ws) / 1000  # [kPa]

    # Constants for Eq. 39
    EQ39_CONST = [6.54, 14.526, 0.7389, 0.09486, 0.4569]

    p_w[p_w <= 0] = 1e-6
    alpha = pd.DataFrame(np.log(p_w))
    alpha = alpha.replace(
        [np.inf, -np.inf], np.NaN).interpolate(method='nearest')

    # Eq. 39
    tdp = alpha.apply(lambda x: EQ39_CONST[0] + EQ39_CONST[1]*x + EQ39_CONST[2]*(x**2) +
           EQ39_CONST[3]*(x**3) + EQ39_CONST[4]*(p_w**0.1984))

    # Eq. 40, TDP less than 0°C and greater than 93°C
    tdp_ice = tdp < 0
    tdp[tdp_ice] = 6.09 + 12.608*alpha[tdp_ice] + 0.4959*(alpha[tdp_ice]**2)

    tdp = tdp.replace(
        [np.inf, -np.inf], np.NaN).interpolate(method='nearest')

    return tdp

# ----------- END tdb2tdp function. -----------
