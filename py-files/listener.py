#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:32:27 2017

@author: rasto

This script is called from the command line. It only parses the arguments
and invokes Indra.
"""

# For parsing the arguments.
import argparse

# Import indra to invoke it later.
from indra import indra

parser = argparse.ArgumentParser(
        description="This is Indra, a generator of synthetic weather " +
        "time series. This function both 'learns' the structure of data " +
        "and samples from the learnt model. Both run modes need 'seed' " +
        "data, i.e., some input weather data.\r\n")

parser.add_argument("--train", type=int, choices=[0, 1], default=0,
                    help="Enter 0 for no seed data, " +
                    "and 1 if you are passing seed data.")
parser.add_argument("--stcode", type=str, default="xxx",
                    help="Make up a station code. " +
                    "If you are not passing seed data, and want me to " +
                    "pick up a saved model, please use the station code" +
                    " of the saved model.")
parser.add_argument("--n_sample", type=int, default=10,
                    help="How many samples do you want out?")
parser.add_argument("--method", type=str, default="arma",
                    help="Which method do you want to use: "
                    "Auto-Regressive Moving Average Models "
                    "(with Fourier pre-processing) [arma], "
                    "or (Experimental) Gaussian Process "
                    "Regression [gp]?")
parser.add_argument("--fpath_in", type=str, help="Path to a folder" +
                    "containing the seed file or the seed file itself." +
                    " If you pass a path to a file, I will only use " +
                    " that file. If you pass a folder path, I will " +
                    "look for files whose names contain the station code.",
                    default="/usr/esru/esp-r/climate/CHE_Geneva_IWEC")
parser.add_argument("--ftype", type=str, help="What kind of file " +
                    "are you giving me? Default is the ESP-r ascii " +
                    "format [espr]. For now, I can read EPW [epw] and " +
                    "ESP-r ascii files. If you pass a plain csv [csv] " +
                    "or python pickle [py] file, it must contain a " +
                    "table with the requisite data in the correct order. " +
                    "See file data_in_spec.txt for the format.",
                    default="espr")
parser.add_argument("--outpath", type=str, default=".",
                    help="Path to the folder where all outputs will go." +
                    " Default is the present working directory.")
parser.add_argument("--l_start", type=int, default=0,
                    help="Starting hour for learning or sampling loop." +
                    " Please key in a multiple of 24, e.g., to start " +
                    "from day 15 of the year, input --l_start 360. " +
                    "Default is 0.")
parser.add_argument("--l_end", type=int, default=31*24,
                    help="Ending hour for learning or sampling loop. " +
                    "Please key in a multiple of 24, e.g., to end on " +
                    "day 30 of the year, input --l_end 720. Default " +
                    "is 31*24.")
parser.add_argument("--l_step", type=int, default=4*24,
                    help="Step size for learning " +
                    "or sampling loop. Please use a multiple of 24, " +
                    "e.g., 24 or 48. Default is 4*24.")
parser.add_argument("--histlim", type=int, default=7*24,
                    help="Maximum 'history' to consider when training " +
                    "the models. The bigger this number, i.e., the " +
                    "longer history you want me to consider, the " +
                    "slower training will be.")
parser.add_argument("--lat", type=int, default=0,
                    help="Station latitude - for now, just used to " +
                    "identify the site. Eventually, this should be " +
                    "used (with longitude) to interpolate or query " +
                    "models for data.")
parser.add_argument("--long", type=int, default=0,
                    help="Station longitude.")
parser.add_argument("--alt", type=int, default=0,
                    help="Station altitude.")
parser.add_argument("--randseed", type=int, default=8760,
                    help="Set the seed for this sampling " +
                    "run. If you don't know what this " +
                    "is, don't worry. The default is 8760.")

args = parser.parse_args()

train = bool(args.train)
stcode = args.stcode.lower()
n_sample = args.n_sample
fpath_in = args.fpath_in
ftype = args.ftype
outpath = args.outpath
l_start = args.l_start
l_end = args.l_end
l_step = args.l_step
histlim = args.histlim
stlat = args.lat
stlong = args.long
stalt = args.alt
randseed = args.randseed
method = args.method

print("Invoking Indra for {0}.\r\n".format(stcode))

if train:
    print("You have asked to train new models. This might take a while.\r\n")
else:
    print("You have asked for samples only. This should be quick.\r\n")

# This command allows this script to be called from the command line.
# Call indra.
if __name__ == "__main__":
    indra(train, stcode=stcode,
          n_sample=n_sample,
          method=method,
          fpath_in=fpath_in,
          ftype=ftype,
          outpath=outpath,
          l_start=l_start, l_end=l_end,
          l_step=l_step, histlim=histlim,
          stlat=stlat, stlong=stlong, stalt=stalt,
          randseed=randseed)
