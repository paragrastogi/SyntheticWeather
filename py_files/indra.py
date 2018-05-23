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
"""

import os
import glob
# import sys
import pickle
import time
# import copy
import pandas as pd
# import numpy as np

# These custom functions load and clean recorded data.
# For now, we are only concerned with ncdc and nsrdb.
import wfileio as wf

from petites import setseed
import resampling as resampling

# Custom functions to calculate error metrics - not currently used.
# import losses.
# from losses import rmseloss
# from losses import maeloss

WEATHER_FMTS = ["espr", "epw", "csv", "fin4"]


def indra(train=False, station_code="abc", n_samples=10,
          path_file_in="wf_in.epw", path_file_out="wf_out.epw",
          file_type="epw", store_path=".",
          climate_change=False, path_cc_file='ccfile.p',
          cc_scenario='rcp85', epoch=[2031, 2040],
          randseed=None, year=0, variant=0,
          arma_params=None,
          bounds=None):

    # Reassign defaults if incoming list params are None
    # (i.e., nothing passed.)
    if arma_params is None:
        arma_params = [2, 2, 1, 1, 24]

    if bounds is None:
        bounds = [0.01, 99.9]

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

    if isinstance(store_path, str):

        # Make a folder named using the station code in case no path to
        # folder was passed.
        if store_path == '.':
            store_path = station_code

        # Store everything in a folder named <station_code>.
        if not os.path.isdir(store_path):
            os.makedirs(store_path)

        # These will be the files where the outputs will be stored.
        path_model_save = os.path.join(
            store_path, 'model_{:d}_{:d}.p'.format(epoch[0], epoch[1]))
        # Save output time series.

        path_syn_save = os.path.join(
            store_path, 'syn_{:d}_{:d}.p'.format(epoch[0], epoch[1]))
        path_counter_save = os.path.join(
            store_path, 'counter_{:d}_{:d}.p'.format(epoch[0], epoch[1]))

    else:
        # This is for the sampling run, where a list of dataframes has
        # been passed.
        path_syn_save = store_path

    # ----------------

    if train:

        # The learning/sampling functions rely on random sampling. For one
        # run, the random seed is constant/immutable; changing it during a
        # run would not make sense. This makes the runs repeatable -- keep
        # track of the seed and you can reproduce exactly the same random
        # number draws as before.

        # If the user did not specify a random seed, then the generator
        # uses the current time, in seconds since some past year, which
        # differs between Unix and Windows. Anyhow, this is saved in the
        # model output in case the results need to be reproduced.
        if randseed is None:
            randseed = int(time.time())

        # Set the seed with either the input random seed or the one
        # assigned just before.
        setseed(randseed)

        # See accompanying script "wfileio".
        # try:
        if os.path.isfile(path_file_in):
            xy_train, locdata, header = wf.get_weather(
                station_code, path_file_in, file_type)

        elif os.path.isdir(path_file_in):

            list_wfiles = ([glob.glob(os.path.join(path_file_in, "*." + x))
                            for x in WEATHER_FMTS] +
                           [glob.glob(os.path.join(path_file_in, "*." + x.upper()))
                            for x in WEATHER_FMTS])
            list_wfiles = sum(list_wfiles, [])

            xy_list = list()

            for file in list_wfiles:
                xy_temp, locdata, header = wf.get_weather(
                    station_code, file, file.split('.')[-1])
                xy_list.append(xy_temp)

            xy_train = pd.concat(xy_list)

        print("Successfully retrieved weather data.\r\n")

        # Train the models.
        print("Training the model. Go get a coffee or something...\r\n")

        if climate_change:

            cc_data = pickle.load(open(path_cc_file, 'rb'))
            cc_data = cc_data[cc_scenario]
            cc_models = set(cc_data.index.get_level_values(0))

            # Pass only the relevant epochs to resampling.
            # For some reason, some models have repetitions and NaNs.
            # This will drop models with no data.

            temp_dict = dict()
            for model in cc_models:

                temp = cc_data.loc[model]
                temp = temp.dropna(how='any')
                # Some times there are non-unique indices, as in duplicate
                # days. Get rid of them by taking the means.
                temp = temp.groupby(temp.index).mean()
                orig_index = temp.index

                if orig_index.shape[0] > 0:
                    temp_dict[model] = temp[
                        (orig_index.year <= epoch[1]) &
                        (orig_index.year >= epoch[0])]

            # import ipdb; ipdb.set_trace()

            cc_data = pd.concat(temp_dict)

        else:
            cc_data = None

        # Hard-coded the scenario as of now - should be added as a
        # parameter later.
        # cc_scenario = 'rcp85'

        # Call resampling with null selmdl and ffit, since those
        # haven"t been trained yet.
        ffit, selmdl, _ = resampling.trainer(
            xy_train, n_samples=n_samples,
            picklepath=path_syn_save,
            arma_params=arma_params,
            bounds=bounds, cc_data=cc_data)

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
                (int(mdl.model.k_seasonal_ar / mdl.model.seasonal_periods),
                 0,
                 int(mdl.model.k_seasonal_ma / mdl.model.seasonal_periods),
                 mdl.model.seasonal_periods)
                for mdl in selmdl]

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
        csave = dict(n_samples=n_samples, randseed=randseed)
        # with open(path_counter_save, "wb") as open_file:
        pickle.dump(csave, open(path_counter_save, "wb"))

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
        # csave = pickle.load(open(path_counter_save, "rb"))

        sample = resampling.sampler(
            picklepath=path_syn_save, year=year, n=variant)

        # import ipdb; ipdb.set_trace()

        if os.path.isdir(path_file_in):

            list_wfiles = [glob.glob(os.path.join(path_file_in, "*." + x))
                           for x in WEATHER_FMTS]
            list_wfiles = sum(list_wfiles, [])

        else:
            list_wfiles = [path_file_in]

        _, locdata, header = wf.get_weather(
                    station_code, list_wfiles[0], file_type)

        # Save / write-out synthetic time series.
        wf.give_weather(sample, locdata, station_code, header,
                        file_type=file_type,
                        path_file_out=path_file_out,
                        masterfile=list_wfiles[0])
