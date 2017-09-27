#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 14:27:53 2017

@author: rasto
"""

import os
import numpy as np
import pandas as pd
from scipy import interpolate


def tdb2tdp(tdb, rh):

    # Change relative humidity to fraction.
    phi = rh/100

    # Remove weird values.
    phi[phi>1] = 1
    phi[phi<0] = 0

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

    lnp_ws = np.zeros(tdb_k.shape)

    # Eq. 5, pg 1.2
    lnp_ws[ice] = c1/tdb_k[ice] + c2 + c3*tdb_k[ice] + c4*tdb_k[ice]**2 + \
        c5*tdb_k[ice]**3 + c6*tdb_k[ice]**4 + c7*np.log(tdb_k[ice])

    # Eq. 6, pg 1.2
    lnp_ws[np.logical_not(ice)] = c8/tdb_k[np.logical_not(ice)] + c9 + \
        c10*tdb_k[np.logical_not(ice)] + \
        c11*tdb_k[np.logical_not(ice)]**2 + \
        (c12*tdb_k[np.logical_not(ice)])**3 + \
        c13*np.log(tdb_k[np.logical_not(ice)])

    # Temperature in the above formulae must be absolute,
    # i.e. in Kelvin

    # THIS step is generating a lot of inf and nan - check against fundamentals handbook.

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
            idx, alpha, kind='nearest', fill_value='extrapolate')
    alpha[duds] = int_func(idx[duds])

    # Eq. 39
    tdp = c14 + c15*alpha + c16*(alpha**2) + c17*(alpha**3) + \
        c18*(p_w**0.1984)

    # Eq. 40, TDP less than 0°C and greater than 93°C
    tdp_ice = np.logical_or(tdp < 0, np.isnan(tdp), np.isinf(tdp))
    tdp[tdp_ice] = 6.09 + 12.608*alpha[tdp_ice] + 0.4959*(alpha[tdp_ice]**2)

    return tdp

# %%


def read_espr(filein='./che_geneva.iwec.a', filetype='ascii'):

    filein_fldr, filein_name = os.path.split(filein)
    sitename = filein_name.split(sep='.')
    sitename = sitename[0]

    filetype = 'ascii'

    if filetype is 'ascii':
        with open(filein) as f:
            content = f.readlines()

    content = [x.strip() for x in content]

    # Split the contents into a header and body.
    header = content[0:13]
    body = content[12:]
    del content

    # Find the lines with day tags.
    daylines = [[idx, line] for [idx, line] in enumerate(body)
                if 'day' in line]

    dataout = np.zeros([8760, 11])

    dcount = 0

    for idx, day in daylines:

        # Get the next 24 lines.
        daylist = np.asarray(body[idx+1:idx+25])

        # Split each line of the current daylist into separate strings.
        splitlist = [element.split() for element in daylist]

        # Convert each element to a integer, then convert the resulting
        # list to a numpy array.
        daydata = np.asarray([list(map(int, x)) for x in splitlist])

        # Today's time slice.
        dayslice = range(dcount, dcount+24, 1)

        # This will split the day-month header line on the gaps.
        splitday = day.split()

        # Month.
        dataout[dayslice, 1] = np.repeat(int(splitday[-1]), len(dayslice))

        # Day.
        dataout[dayslice, 1] = np.repeat(int(splitday[2]), len(dayslice))

        # Hour (of day).
        dataout[dayslice, 2] = np.arange(0, 24, 1)

        # tdb, input is in deci-degrees, convert to degrees.
        dataout[dayslice, 3] = daydata[:, 1]/10

        # rh, in percent.
        dataout[dayslice, 5] = daydata[:, 5]

        # ghi, in W/m2.
        dataout[dayslice, 7] = daydata[:, 0]

        # dni, in W/m2.
        dataout[dayslice, 8] = daydata[:, 2]

        # wdr, input is in deci-m/s.
        dataout[dayslice, 9] = daydata[:, 3]/10

        # wspd, clockwise deg from north.
        dataout[dayslice, 10] = daydata[:, 4]


        dcount += 24

    # tdp, calculated from tdb and rh.
    dataout[:, 4] = tdb2tdp(dataout[:, 3], dataout[:, 5])


    clmdata = pd.DataFrame(data=dataout,
            columns=['month', 'day', 'hour',
                     'tdb', 'tdp', 'rh',
                     'ghi', 'dni', 'dhi',
                     'wspd', 'wdr'])

    return clmdata



#def writeclm(fileout, filetype='ascii'):
#src = 'indra'
## filein_name and sitename taken from incoming file name.
## Set year to some ridiculous value so it's never mistaken for an actual year.
#year = 2999
#diff_rad_flag = 0
#period_st = 1
#period_end = 365
#
##espr_header = str('*CLIMATE\r\n' +
##'# ascii climate file from {0},\r\n' +
##'# defined in: {1}\r\n' +
##'# col 1: Diffuse solar on the horizontal (W/m**2)\r\n' +
##'# col 2: External dry bulb temperature   (Tenths DEG.C)\r\n' +
##'# col 3: Direct normal solar intensity   (W/m**2)\r\n' +
##'# col 4: Prevailing wind speed           (Tenths m/s) \r\n' +
##'# col 5: Wind direction     (clockwise deg from north) \r\n' +
##'# col 6: Relative humidity               (Percent)\r\n' +
##'{2}                # site name\r\n' +
##' {3},{4},{5},{6},  # year latitude, long diff, rad flag\r\n' +
##' {7},{8},    # period (julian days)\r\n').format(src, filein_name, sitename, year, lat, long, diff_rad_flag, period_st, period_end)


