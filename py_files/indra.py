"""
Create Synthetic Weather based on some recorded data.
The algorithm works by creating synthetic time series over
short periods based on short histories. These short series
may be comined to obtain a longer series.
Script originally written by Parag Rastogi. Started: July 2017
@author = Parag Rastogi

Description of algorithm:
    1. Load data.
    2. Scale data using a standard scaler (subtract mean and divide by std).
    3. Enter model-fitting loop:
        a. Select 14 days of history (less when beginning).
        b. Use history to train model for next day.
        c. Sample from next day to obtain synthetic "predictions".
        d. Once the "predictions" are obtained for every day of the year,
           we are left with synthetic time series.
    4. Un-scale the data using the same scaler as (2) above.
    5. Clean / post-process the data if needed, e.g., oscillation of solar
       values around sunrise and sunset.

Features to be implemented:
    1. Load climate change model information and add that to the model.
       The outputs of climate models are usually in the form of mean changes.
       I"ve added them to my models in the past - need to figure out best way
       to do this now.
    2. Ability to get recorded or typical data automatically from some web
       service if the user provides coordinates or WMO station number.
    3. Ability to download climate change data the same way.
    4. Ability to learn a transfer function between urban and rural stations
       for quick-and-not-too-dirty estimates of urban heat island effects.
       So long as time series of a year-ish are available from the locations
       of interest, the transfer _should_ be able to proceed without requiring
       too much information about the local urban canopy.
"""

import os
import sys
import pickle
import time

# These custom functions load and clean recorded data.
# For now, we are only concerned with ncdc and nsrdb.
import wfileio as wf

from petites import setseed
from resampling import resampling

# Custom functions to calculate error metrics - not currently used.
# import losses.
# from losses import rmseloss
# from losses import maeloss


def indra(train=False, station_code="abc", n_samples=10,
          path_file_in="wf_in.a", path_file_out="wf_out.a",
          file_type="espr", store_path=".",
          station_coordinates=None,
          randseed=None,
          arma_params=None,
          bounds=None):

    # Reassign defaults if incoming list params are None
    # (i.e., nothing passed.)
    if arma_params is None:
        arma_params = [2, 2, 1, 1, 24]

    if bounds is None:
        bounds = [0.01, 99.9]

    if station_coordinates is None:
        station_coordinates = [0.0, 0.0, 0.0]

    # Uncomment when debugging this script to avoid having to call the
    # whole function.
    # train=False
    # station_code="lgw"
    # n_samples=100
    # store_path="lgw"
    # path_file_in= os.path.join(store_path, "GBR_London_Gatwick.a")
    # path_file_out= os.path.join(store_path, "GBR_London_Gatwick_syn.a")
    # file_type="espr"
    # randseed=None

    # ------------------
    # Some initialisation house work.

    # Convert incoming station_code to lowercase.
    station_code = station_code.lower()

    # Make a folder named using the station code in case no path to
    # folder was passed.
    if store_path == '.':
        store_path = station_code

    # Store everything in a folder named <station_code>.
    if not os.path.isdir(store_path):
        os.makedirs(store_path)

    # These will be the files where the outputs will be stored.
    path_model_save = os.path.join(store_path, "model.p")
    # Save output time series.
    picklepath = os.path.join(store_path, "syn.npy")
    path_counter_save = os.path.join(store_path, "counter.p")

    # ----------------

    # See accompanying script "wfileio".
    try:
        xy_train, locdata, header = wf.get_weather(
            station_code, path_file_in, file_type)

        print("Successfully retrieved weather data.\r\n")

    except Exception as err:
        print("Error: " + str(err))
        exc_type, _, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        print("I could not read the incoming weather file. " +
              "Terminating this run.\r\n")
        # return 0

    if train:

        # The learning/sampling functions rely on random sampling. For one
        # run, the random seed is constant/immutable; changing it during a
        # run would not make sense. This makes the runs repeatable -- keep
        # track of the seed and you can reproduce exactly the same random
        # number draws as before.

        # If the user did not specify a random seed, then the generator
        # uses the current time, in seconds since an epoch, which differs
        # between Unix and Windows. Anyhow, this is saved in the model
        # output in case the results need to be reproduced.
        if randseed is None:
            randseed = int(time.time())

        # Set the seed with either the input random seed or the one
        # assigned just before.
        setseed(randseed)

        # Train the models.
        print("Training the model. Go get a coffee or something...\r\n")

        # Call resampling with null selmdl and ffit, since those
        # haven"t been trained yet.
        ffit, selmdl, _ = resampling(
            xy_train, train=True, n_samples=n_samples,
            picklepath=picklepath,
            arma_params=arma_params,
            bounds=bounds)

        # The non-seasonal order of the model. This exists in both
        # ARIMA and SARIMAX models, so it has to exist in the output
        # of resampling.
        order = [(p.model.k_ar, 0, p.model.k_ma) for p in selmdl]
        # Also the endogenous variable.
        endog = [p.model.endog for p in selmdl]

        params = [p.params for p in selmdl]

        try:
            # Try to find the seasonal order. If it exists, save the
            # sarimax model. This should almost always be the case.
            seasonal_order = [
                (int(mdl.model.k_seasonal_ar/mdl.model.seasonal_periods),
                 0,
                 int(mdl.model.k_seasonal_ma/mdl.model.seasonal_periods),
                 mdl.model.seasonal_periods)
                for mdl in selmdl
                ]

            arma_save = dict(order=order, params=params,
                             seasonal_order=seasonal_order,
                             ffit=ffit, endog=endog,
                             randseed=randseed)

        except Exception:
            # Otherwise, ask for forgiveness and save the ARIMA model.
            arma_save = dict(order=order, params=params, endog=endog,
                             ffit=ffit, randseed=randseed)

        with open(path_model_save, "wb") as open_file:
            pickle.dump(arma_save, open_file)

        # Save counter.
        csave = dict(counter=0, n_samples=n_samples)
        with open(path_counter_save, "wb") as open_file:
            pickle.dump(csave, open_file)

        print(("I've saved the model for station '{0}'. "
               "You can now ask me for samples in folder '{1}'."
               "\r\n").format(station_code, store_path))

    else:

        # Call the functions in sampling mode.

        # The output, xout, is a numpy nd-array with the standard
        # columns ("month", "day", "hour", "tdb", "tdp", "rh",
        # "ghi", "dni", "dhi", "wspd", "wdr")

        # In this MC framework, the "year" of weather data is meaningless.
        # When the climate change models will be added, these years will
        # mean something. For now, any number will do.

        # Load counter.
        with open(path_counter_save, "rb") as open_file:
            csave = pickle.load(open_file)

        _, _, xout = resampling(
            xy_train, counter=csave["counter"], train=False, sample=True)

        # Save / write-out synthetic time series.
        wf.give_weather(xout, locdata, station_code, header, file_type=file_type,
                        s_shift=0, path_file_out=path_file_out, masterfile=path_file_in)

        # This function has been asked to give a sample, so update
        # the counter.
        csave["counter"] += 1
        if csave["counter"] >= (csave["n_samples"]-1):
            csave["counter"] = 0

        with open(path_counter_save, "wb") as open_file:
            pickle.dump(csave, open_file)
