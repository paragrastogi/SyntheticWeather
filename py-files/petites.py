#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 15:30:34 2017

@author: parag rastogi
"""

import np
from scipy import interpolate


def wstats(data, key, stat):

    a = data.groupby(key)

    if stat is 'mean':
        b = a.mean()
    elif stat is 'sum':
        b = a.sum()
    elif stat is 'max':
        b = a.max()
    elif stat is 'min':
        b = a.min()
    elif stat is 'std':
        b = a.std()
    elif stat is 'q1':
        b = a.quantile(0.25)
    elif stat is 'q3':
        b = a.quantile(0.75)
    elif stat is 'med':
        b = a.median()

    return b

# ----------- END wstats function. -----------


def tdb2tdp(tdb, rh):

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
    c1 = -5.6745359*10**3
    c2 = 6.3925247
    c3 = -9.6778430*10**-3
    c4 = 6.2215701*10**-7
    c5 = 2.0747825*10**-9
    c6 = -9.4840240*10**-13
    c7 = 4.1635019

    # Constants for Eq. 6, Temperature 0°C to 200°C.
    c8 = -5.8002206*10**3
    c9 = +1.3914993
    c10 = -4.8640239*10**-2
    c11 = +4.1764768*10**-5
    c12 = -1.4452093*10**-8
    c13 = +6.5459673

    # This is to distinguish between the two versions of equation 5.
    ice = tdb_k <= 273.15
    not_ice = np.logical_not(ice)

    lnp_ws = np.zeros(tdb_k.shape)

    # Eq. 5, pg 1.2
    lnp_ws[ice] = (c1/tdb_k[ice] + c2 + c3*tdb_k[ice] + c4*tdb_k[ice]**2 +
                   c5*tdb_k[ice]**3 + c6*tdb_k[ice]**4 +
                   c7*np.log(tdb_k[ice]))

    # Eq. 6, pg 1.2
    lnp_ws[np.logical_not(ice)] = (
            c8/tdb_k[not_ice] + c9 +
            c10*tdb_k[not_ice] +
            c11*tdb_k[not_ice]**2 +
            c12*tdb_k[not_ice]**3 +
            c13*np.log(tdb_k[np.logical_not(ice)]))

    # Temperature in the above formulae must be absolute,
    # i.e. in Kelvin

    # Continuing from eqs. 5 and 6
    p_ws = np.e**(lnp_ws)  # [Pa]

    # Eq. 24, pg 1.8
    p_w = (phi * p_ws) / 1000  # [kPa]

    # Constants for Eq. 39
    c14 = 6.54
    c15 = 14.526
    c16 = 0.7389
    c17 = 0.09486
    c18 = 0.4569

    alpha = np.log(p_w)
    idx = np.arange(0, alpha.size)
    duds = np.logical_or(np.isinf(alpha), np.isnan(alpha))
    int_func = interpolate.interp1d(
            idx[np.logical_not(duds)], alpha[np.logical_not(duds)],
            kind='nearest', fill_value='extrapolate')
    alpha[duds] = int_func(idx[duds])

    # Eq. 39
    tdp = c14 + c15*alpha + c16*(alpha**2) + c17*(alpha**3) + \
        c18*(p_w**0.1984)

    # Eq. 40, TDP less than 0°C and greater than 93°C
    tdp_ice = np.logical_or(tdp < 0, np.isnan(tdp), np.isinf(tdp))
    tdp[tdp_ice] = 6.09 + 12.608*alpha[tdp_ice] + 0.4959*(alpha[tdp_ice]**2)

    return tdp

# ----------- END tdb2tdp function. -----------
