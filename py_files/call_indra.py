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

# Define a parser.
PARSER = argparse.ArgumentParser(
    description="This is INDRA, a generator of synthetic weather " +
    "time series. This function both 'learns' the structure of data " +
    "and samples from the learnt model. Both run modes need 'seed' " +
    "data, i.e., some input weather data.\r\n", prog='INDRA',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

PARSER.add_argument("--train", type=int, choices=[0, 1], default=0,
                    help="Enter 0 for no seed data (sampling mode), " +
                    "or 1 if you are passing seed data (training or " +
                    "initalisation mode).")
PARSER.add_argument("--station_code", type=str, default="abc",
                    help="Make up a station code. " +
                    "If you are not passing seed data, and want me to " +
                    "pick up a saved model, please use the station code" +
                    " of the saved model.")
PARSER.add_argument("--n_samples", type=int, default=10,
                    help="How many samples do you want out?")
PARSER.add_argument("--path_file_in", type=str, help="Path to a weather " +
                    "file (seed file).", default="wf_in.a")
PARSER.add_argument("--path_file_out", type=str, help="Path to where the " +
                    "synthetic data will be written. If you ask for more " +
                    "than one sample, I will append an integer to the name.",
                    default="wf_out.a")
PARSER.add_argument("--file_type", type=str, default="espr",
                    help=("What kind of input weather file "
                          "are you giving me? Default is the ESP-r ascii "
                          "format [espr]. For now, I can read EPW [epw] and "
                          "ESP-r ascii files. If you pass a plain csv [csv] "
                          "or python pickle [py] file, it must contain a "
                          "table with the requisite data in the correct "
                          "order. See file data_in_spec.txt for the format."))
# Indra needs the data to be a numpy nd-array arranged exactly so:
# month, day of year, hour, tdb, tdp, rh, ghi, dni, dhi, wspd, wdr
PARSER.add_argument("--store_path", type=str, default="SyntheticWeather",
                    help="Path to the folder where all outputs will go." +
                    " Default behaviour is to create a folder in the " +
                    "present working directory called SyntheticWeather.")
PARSER.add_argument("--climate_change", type=int, choices=[0, 1], default=0,
                    help="Enter 0 to not include climate change models, or" +
                    " 1 to do so. If you want to use a CC model, you have" +
                    " to pass a path to the file containing those outputs.")
PARSER.add_argument("--path_cc_file", type=str, default="ccfile.p",
                    help="Path to the file containing CC model outputs.")
PARSER.add_argument("--station_coordinates", type=str, default="[0, 0, 0]",
                    help="Station latitude, longitude, altitude. " +
                    "Not currently used.")
PARSER.add_argument("--randseed", type=int, default=42,
                    help="Set the seed for this sampling " +
                    "run. If you don't know what this " +
                    "is, don't worry. The default is 42. Obviously.")
PARSER.add_argument("--arma_params", type=str, default="[2,2,1,1,24]",
                    help=("A list of UPPER LIMITS of the number of SARMA "
                          "terms [AR, MA, Seasonal AR, Seasonal MA, "
                          "Seasonality] to use in the model. Input should "
                          "look like a python list, i.e., [a,b,c], WITHOUT "
                          "SPACES. If you don't know what this is, " +
                          "don't worry. The default is [2,2,1,1,24]. "
                          "The default frequency of Indra is hours, so "
                          "seasonality should be declared in hours."))
PARSER.add_argument("--bounds", type=str, default="[0.01,99.9]",
                    help=("Lower and upper bound percentile values to "
                          "use for cleaning the synthetic data. Input "
                          "should look like a python list, i.e., [a,b,c], "
                          "WITHOUT SPACES. The defaults bounds are the "
                          "0.01 and 99.9 percentiles, i.e., [0.01,99.9]."))

ARGS = PARSER.parse_args()

train = bool(ARGS.train)
station_code = ARGS.station_code.lower()
n_samples = ARGS.n_samples
path_file_in = ARGS.path_file_in
path_file_out = ARGS.path_file_out
file_type = ARGS.file_type
store_path = ARGS.store_path
station_coordinates = [float(x.strip("[").strip("]"))
                       for x in ARGS.station_coordinates.split(",")]
climate_change = ARGS.climate_change
path_cc_file = ARGS.path_cc_file
randseed = ARGS.randseed
arma_params = [int(x.strip("[").strip("]"))
               for x in ARGS.arma_params.split(",")]
bounds = [float(x.strip("[").strip("]")) for x in ARGS.bounds.split(",")]

print("\r\nInvoking indra for {0}.\r\n".format(station_code))

if store_path == "SyntheticWeather":
    store_path = store_path + '_' + station_code

# Call indra using the processed arguments.
if __name__ == "__main__":
    indra(train, station_code=station_code,
          n_samples=n_samples,
          path_file_in=path_file_in,
          path_file_out=path_file_out,
          file_type=file_type,
          store_path=store_path,
          climate_change=climate_change,
          path_cc_file=path_cc_file,
          station_coordinates=station_coordinates,
          randseed=randseed,
          arma_params=arma_params,
          bounds=bounds)
