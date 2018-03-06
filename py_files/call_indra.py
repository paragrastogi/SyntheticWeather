#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:32:27 2017

@author: Parag Rastogi

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
                    help="Enter 0 for no seed data (sampling mode), " +
                    "or 1 if you are passing seed data (training or " +
                    "initalisation mode).")
parser.add_argument("--stcode", type=str, default="abc",
                    help="Make up a station code. " +
                    "If you are not passing seed data, and want me to " +
                    "pick up a saved model, please use the station code" +
                    " of the saved model.")
parser.add_argument("--n_samples", type=int, default=10,
                    help="How many samples do you want out?")
parser.add_argument("--fpath_in", type=str, help="Path to a weather " +
                    "file (seed file).", default="wf_in.a")
parser.add_argument("--fpath_out", type=str, help="Path to where the " +
                    "synthetic data will be written. If you ask for more " +
                    "than one sample, I will append an integer to the name.",
                    default="wf_out.a")
parser.add_argument("--ftype", type=str, help="What kind of file " +
                    "are you giving me? Default is the ESP-r ascii " +
                    "format [espr]. For now, I can read EPW [epw] and " +
                    "ESP-r ascii files. If you pass a plain csv [csv] " +
                    "or python pickle [py] file, it must contain a " +
                    "table with the requisite data in the correct order. " +
                    "See file data_in_spec.txt for the format.",
                    default="espr")
# Indra needs the data to be a numpy nd-array arranged exactly so:
# month, day of year, hour, tdb, tdp, rh, ghi, dni, dhi, wspd, wdr
parser.add_argument("--storepath", type=str, default="SyntheticWeather",
                    help="Path to the folder where all outputs will go." +
                    " Default behaviour is to create a folder in the " +
                    "present working directory called SyntheticWeather.")
parser.add_argument("--lat", type=int, default=0,
                    help="Station latitude - for now, just used to " +
                    "identify the site.")
parser.add_argument("--long", type=int, default=0,
                    help="Station longitude.")
parser.add_argument("--alt", type=int, default=0,
                    help="Station altitude.")
parser.add_argument("--randseed", type=int, default=42,
                    help="Set the seed for this sampling " +
                    "run. If you don't know what this " +
                    "is, don't worry. The default is 42.")
parser.add_argument("--arp_ub", type=int, default=4,
                    help="Upper limit of the number of AR terms " +
                    "to use in the model. If you don't know what this " +
                    "is, don't worry. The default is 4.")
parser.add_argument("--maq_ub", type=int, default=4,
                    help="Upper limit of the number of MA terms " +
                    "to use in the model. If you don't know what this " +
                    "is, don't worry. The default is 4.")
parser.add_argument("--sarp_ub", type=int, default=1,
                    help="Upper limit of the number of Seasonal AR terms " +
                    "to use in the model. If you don't know what this " +
                    "is, don't worry. The default is 1.")
parser.add_argument("--smaq_ub", type=int, default=1,
                    help="Upper limit of the number of Seasonal MA terms " +
                    "to use in the model. If you don't know what this " +
                    "is, don't worry. The default is 1.")
parser.add_argument("--seasonality", type=int, default=24,
                    help="The period of the seasonal terms (seasonality) " +
                    "to use in the model. If you don't know what this " +
                    "is, don't worry. The default is 24.")

args = parser.parse_args()

train = bool(args.train)
stcode = args.stcode.lower()
n_samples = args.n_samples
fpath_in = args.fpath_in
fpath_out = args.fpath_out
ftype = args.ftype
storepath = args.storepath
stlat = args.lat
stlong = args.long
stalt = args.alt
randseed = args.randseed
arp_ub = args.arp_ub
maq_ub = args.maq_ub
sarp_ub = args.sarp_ub,
smaq_ub = args.smaq_ub
seasonality = args.seasonality

print("\r\nInvoking indra for {0}.\r\n".format(stcode))

if storepath == "SyntheticWeather":
    storepath = storepath + stcode

if __name__ == "__main__":
    indra(train, stcode=stcode,
          n_samples=n_samples,
          fpath_in=fpath_in,
          fpath_out=fpath_out,
          ftype=ftype,
          storepath=storepath,
          stlat=stlat, stlong=stlong, stalt=stalt,
          randseed=randseed,
          arp_ub=arp_ub, maq_ub=maq_ub, sarp_ub=sarp_ub,
          smaq_ub=smaq_ub, seasonality=seasonality)
