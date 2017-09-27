#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:32:27 2017

@author: rasto
"""


import argparse

# Invoke indra.
from indra import indra

parser = argparse.ArgumentParser(
        description='This is Indra, a generator of ' +
        'synthetic weather time series.')

parser.add_argument('--data_in', type=int, choices=[0, 1], default=0,
                    help='Enter 0 for no seed data, ' +
                    'and 1 if you are passing seed data.')
parser.add_argument('--stcode', type=str, default='xxx',
                    help='Make up a station code. ' +
                    'If you are not passing seed data, and want me to ' +
                    'pick up a saved model, please use the station code' +
                    ' of the saved model.')
parser.add_argument('--n_samples', type=int, default=10,
                    help='How many samples do you want out?')
parser.add_argument('--path_wthr', type=str, help='Path to a folder' +
                    'containing the seed file or the seed file itself.' +
                    ' If you pass a path to a file, I will only use ' +
                    ' that file. If you pass a folder path, I will ' +
                    'look for files whose names contain the station code.',
                    default='/usr/esru/esp-r/climate')
parser.add_argument('--outpath', type=str, default='.',
                    help='Path to the folder where all outputs will go.' +
                    ' Default is the present working directory.')
parser.add_argument('--l_start', type=int, default=0,
                    help='Starting hour for learning or sampling loop.' +
                    ' Please key in a multiple of 24, e.g., to start ' +
                    'from day 15 of the year, input --l_start 360. ' +
                    'Default is 0.')
parser.add_argument('--l_end', type=int, default=31*24,
                    help='Ending hour for learning or sampling loop. ' +
                    'Please key in a multiple of 24, e.g., to end on ' +
                    'day 30 of the year, input --l_end 720. Default ' +
                    'is 31*24.')
parser.add_argument('--l_step', type=int, default=4*24,
                    help='Step size for learning ' +
                    'or sampling loop. Please use a multiple of 24, ' +
                    'e.g., 24 or 48. Default is 4*24.')
parser.add_argument('--histlim', type=int, default=14*24,
                    help="Maximum 'history' to consider when training " +
                    "the models. The bigger this number, i.e., the " +
                    "longer history you want me to consider, the " +
                    "slower training will be.")
parser.add_argument('--lat', type=int, default=0,
                    help='Station latitude - for now, just used to ' +
                    'identify the site. Eventually, this should be ' +
                    'used (with longitude) to interpolate or query ' +
                    'models for data.')
parser.add_argument('--long', type=int, default=0,
                    help='Station longitude.')
parser.add_argument('--alt', type=int, default=0,
                    help='Station altitude.')
parser.add_argument('--randomseed', type=int, default=8760,
                    help="Set the seed for this sampling " +
                    "run. If you don't know what this " +
                    "is, don't worry. The default is 8760.")

args = parser.parse_args()

data_in = bool(args.data_in)
stcode = args.stcode.lower()
n_samples = args.n_samples
path_wthr = args.path_wthr
outpath = args.outpath
l_start = args.l_start
l_end = args.l_end
l_step = args.l_step
histlim = args.histlim
stlat = args.lat
stlong = args.long
stalt = args.alt
randomseed = args.randomseed

if data_in:
    print('Invoking Indra for {0} with seed data.'.format(stcode))
else:
    print('Invoking Indra for {0} without seed data.'.format(stcode))

# This command allows this script to be called from the command line.
if __name__ == '__main__':
    indra(data_in, stcode=stcode,
          n_samples=n_samples,
          path_wthr=path_wthr,
          outpath=outpath,
          l_start=l_start, l_end=l_end,
          l_step=l_step, histlim=histlim,
          stlat=stlat, stlong=stlong, stalt=stalt,
          randomseed=randomseed)
