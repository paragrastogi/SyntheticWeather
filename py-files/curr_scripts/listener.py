#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:32:27 2017

@author: rasto
"""


import argparse

# Invoke indra.
from indra import synweather

parser = argparse.ArgumentParser(description='You have invoked Indra.')

parser.add_argument('seedfile', dest='seedfile', type=bool)
parser.add_argument('--stcode=gla', dest='stcode')

parser.parse_args('-seedfile', '--stcode')

args = parser.parse_args()

print(seedfile)
print(stcode)

exit

# This command allows this script to be called from the command line.
if __name__ == '__main__':
    synweather(seedfile, stcode='gla', n_samples=10,
       path_wthr_fldr='/usr/esru/esp-r/climate',
       outpath='.',
       l_start=int(0), l_end=int(31*24),
       l_step=int(4*24), histlim=int(14*24),
       stlat=0.0, stlong=0.0, stalt=0.0,
       randomseed=8760)

#    if (len(inputlist) == 3):
#        synweather(sys.argv[1], stcode=sys.argv[2],
#                   path_wthr_fldr='/home/rasto/Documents/WeatherData/HistoricalData/gen')
#    elif (len(inputlist) == 4):
#        synweather(sys.argv[1], stcode=sys.argv[2],
#                   path_wthr_fldr=sys.argv[3])
#
#    elif (len(inputlist) == 5):
#        synweather(sys.argv[1], stcode=sys.argv[2],
#                   path_wthr_fldr=sys.argv[3])