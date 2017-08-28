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

# Load data about the cities. This is just for this example.
citytab = pd.read_csv(os.path.join('..', 'CityData.csv'),
                      dtype=dict(WMO=str, StCode=str))
# The following inputs would be input by the user.
# Station code.
stcode = citytab.StCode[41]
# Longitude.
stlong = 40.78
# Latitude.
stlat = -73.7
# Altitude.
stalt = 47.5
# Specify the sources of the actual data - please follow AMY keywords list.
sources = ('ncdc', 'nsrdb')

# See accompanying script "gw".
#typicaldata, actualdata = gw.get_weather(stcode, citytab, sources)

picklename = str(stlong) + ',' + str(stlat)