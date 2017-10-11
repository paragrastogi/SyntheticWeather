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

# Custom functions to calculate error metrics - not currently used.
# import losses.
# from losses import rmseloss
# from losses import maeloss

# import learn and sample functions from gp_funcs.
from gp_funcs import learngp
from gp_funcs import samplegp

from petites import setseed
from resampling import resampling

# from statsmodels.tsa.arima_model import ARIMAResults
# from statsmodels.tsa.statespace.sarimax import SARIMAX


def indra(train=False, stcode="abc", n_sample=100, method="arma",
          fpath_in="wf_in.a", fpath_out="wf_out.a",
          ftype="espr", storepath=".",
          cc=False, ccpath=".",
          l_start=int(0), l_end=int(31*24),
          l_step=int(4*24), histlim=int(14*24),
          stlat=0.0, stlong=0.0, stalt=0.0,
          randseed=None):

    # Uncomment when debugging this script to avoid having to call the
    # whole function.
    train=False
    stcode="gen"
    n_sample=1
    method="arma"
    storepath="SyntheticWeather-gen"
    fpath_in= os.path.join(storepath, "che_geneva.iwec.a")
    fpath_out= os.path.join(storepath, "che_geneva.iwec_syn.a")
    ftype="espr"
    cc=False
    ccpath="."
    randseed=None

    # ------------------
    # Some initialisation house work.

    # The learning/sampling functions rely on random sampling. For one run,
    # the random seed is constant/immutable; changing it during a run
    # would not make sense. This makes the runs repeatable -- keep track of
    # the seed and you can reproduce exactly the same random number draws
    # as before.

    # If the user did not specify a random seed, then the generator
    # uses the current time, in seconds since an epoch, which differs
    # between Unix and Windows. Anyhow, this is saved in the model
    # output in case the results need to be reproduced.
    if randseed is None:
        randseed = int(time.time())

    # Set the seed with either the input random seed or the one
    # assigned just before.
    setseed(randseed)

    # Convert incoming stcode to lowercase.
    stcode = stcode.lower()

    if storepath == '.':
        storepath = stcode

    # Store everything in a folder named <stcode>.
    if not os.path.isdir(storepath):
        os.makedirs(storepath)

    # These will be the files where the outputs will be stored.
    path_model_save = os.path.join(storepath, "model.p")
    # Save output time series.
    picklepath = os.path.join(storepath, "syn.npy")
    path_counter_save = os.path.join(storepath, "counter.p")

    # ----------------

    # Load data about the cities. This isn"t usually necessary, so it
    # doesn"t matter if this file isn"t present.
    #    try:
    #        citytab = pd.read_csv(os.path.join("CityData.csv"),
    #                              dtype=dict(WMO=str, StCode=str))
    #    except:
    #        citytab = None
    #        print("I could not find CityData.csv, continuing without...")

    # See accompanying script "wfileio".
    try:
        xy_train, locdata, header = wf.get_weather(
                stcode, fpath_in, ftype)

        # The GP method needs day of year rather than day of month.
        if method == "gp":
            temp = wf.day_of_year(xy_train[:, 1], xy_train[:, 2])
            xy_train[:, 2] = temp
            del temp

        print("Successfully retrieved weather data.\r\n")

    except Exception as err:
        print("Error: " + str(err))
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        print("I could not read the incoming weather file. " +
              "Terminating this run.\r\n")
        # return 0
# %%

    if train:

        # Train the models.

        if method == "gp":

            gp_list, mtrack, scaler = learngp(
                    l_start, l_end, l_step, histlim, xy_train)

            # Save gp_list and month_tracker to pickle file.
            gp_save = dict(gp_list=gp_list, mtrack=mtrack,
                           scaler=scaler, xy_train=xy_train)

            with open(path_model_save, "wb") as fp:
                pickle.dump(gp_save, fp)

        elif method == "arma":

            # Call resampling with null selmdl and ffit, since those
            # haven"t been trained yet.
            ffit, selmdl, _, _ = resampling(
                    xy_train, selmdl=None, ffit=None,
                    train=True, sample=False, n_sample=n_sample,
                    picklepath=picklepath, randseed=randseed)

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
                        (int(p.model.k_seasonal_ar/p.model.seasonal_periods),
                         0,
                         int(p.model.k_seasonal_ma/p.model.seasonal_periods),
                         p.model.seasonal_periods) for p in selmdl]

                arma_save = dict(order=order, params=params,
                                 seasonal_order=seasonal_order,
                                 ffit=ffit, endog=endog,
                                 randseed=randseed)

            except AttributeError:
                # Otherwise, ask for forgiveness and save the ARIMA model.
                arma_save = dict(order=order, params=params, endog=endog,
                                 ffit=ffit, randseed=randseed)

            with open(path_model_save, "wb") as fp:
                pickle.dump(arma_save, fp)

            # Save counter.
            counter = 0
            with open(path_counter_save, "wb") as fp:
                pickle.dump(counter, fp)

    else:

        # Load counter.
        with open(path_counter_save, "rb") as fp:
            counter = pickle.load(fp)


        if method == "arma":

            _, _, xout = resampling(
                    xy_train, counter=counter, selmdl=None, ffit=None,
                    train=False, sample=True, n_sample=n_sample,
                    picklepath=picklepath, randseed=randseed)

        elif method == "gp":

            with open(path_model_save, "rb") as fp:
                gp_save = pickle.load(fp)

            xout = samplegp(gp_list, l_start, l_end, l_step, histlim,
                            n_sample, xy_train, mtrack, scaler,
                            picklepath=picklepath)

            for n in range(0, xout.shape[-1]):

                temp1, temp2 = wf.day_of_month(xout[:, 2, n])
                xout[:, 2, n] = temp2
                del [temp1, temp2]

        # End of if method statement.

        # Save / write-out synthetic time series.
        wf.give_weather(xout, locdata, stcode, header, ftype=ftype,
                        s_shift=0, fpath_out=fpath_out, masterfile=fpath_in)

        # This function has been asked to give a sample, so update
        # the counter.
        counter += 1
        with open(path_counter_save, "wb") as fp:
            pickle.dump(counter, fp)

#        # Load models from file.
#
#        if method == "gp":
#
#            with open(path_model_save, "rb") as fp:
#                gp_save = pickle.load(fp)
#
#            gp_list = gp_save["gp_list"]
#            mtrack = gp_save["mtrack"]
#            scaler = gp_save["scaler"]
#            xy_train = gp_save["xy_train"]
#
#        elif method == "arma":
#
#            with open(path_model_save, "rb") as fp:
#                arma_save = pickle.load(fp)
#
#            selmdl = list()
#
#            for (o, so, e, p) in zip(arma_save["order"],
#                                     arma_save["seasonal_order"],
#                                     arma_save["endog"],
#                                     arma_save["params"]):
#                mod_temp = SARIMAX(e, order=o, params=p,
#                                   seasonal_order=so,
#                                   trend=None)
#                selmdl.append(mod_temp.fit(disp=0))
#
#            ffit = list()
#
#            for f in arma_save["ffit"]:
#                ffit.append(f)

    # %%

    # Call the sampling function.

    # The output, xout, is a numpy nd-array with the standard
    # columns ("month", "day", "hour", "tdb", "tdp", "rh",
    # "ghi", "dni", "dhi", "wspd", "wdr")

    # In this MC framework, the "year" of weather data is meaningless.
    # When the climate change models will be added, these years will
    # mean something. For now, just add "2017" to every file.
