# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 10:41:29 2017

@author: Parag Rastogi

This script is meant to "listen" for a request for weather files.

The script will begin by locating/recognising the location
(represented by latitude and longitude). Then, it can take
one of the following actions:
    1. Search for a pickle file containing pre-generated synthetic
    files. -> Load pickle file. -> Return as many time series as
    were asked.
    2.


"""

__author__ = 'Parag Rastogi'

import os
import numpy as np
import pandas as pd

# This bit of code loads data from pre-processed pickle files or
# processes CSV files. This is the 'starter' data that the script needs.

# Where am I to search for the pickle files?
picklepath = os.path.join('/home', 'rasto', 'Documents', 'WeatherData',
                          'SyntheticData', 'pickled')
picklelist = os.listdir(picklepath)

# Load data about the cities. This is just for this example.
citytab = pd.read_csv(os.path.join('..', 'CityData.csv'),
                      dtype=dict(WMO=str, StCode=str))
# The following inputs would be input by the user.
# Longitude.
stlong = -73.76
# Latitude.
stlat = 40.76

# Take the entered latitude and longitude, and first search for a match
# with all the significant figures.
all_lats = citytab.Latitude
all_longs = citytab.Longitude

latmatch = stlat == all_lats
longmatch = stlong == all_longs

stfind = np.zeros(latmatch.shape, dtype=np.bool)
patience = 0

while not any(stfind):

    stfind = latmatch & longmatch

    if any(stfind):
        # Station code.
        stcode = citytab.StCode[stfind]  # [41]
        # Altitude.
        stalt = citytab.Altitude[stfind]

    else:
        stlat = round(stlat, patience)
        stlong = round(stlong, patience)

        all_longs = round(all_longs, patience)
        all_lats = round(all_lats, patience)

        latmatch = stlat == all_lats
        longmatch = stlong == all_longs

        stfind = latmatch & longmatch

        stcode = citytab.StCode[stfind]
        stalt = citytab.Altitude[stfind]

    patience += 1

    if patience == 2:
        break



# Specify the sources of the actual data - please follow AMY keywords list.
# sources = ('ncdc', 'nsrdb')
# See accompanying script "gw".
# typicaldata, actualdata = gw.get_weather(stcode, citytab, sources)

# picklename = str(stlong) + ',' + str(stlat)