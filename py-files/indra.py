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
       I've added them to my models in the past - need to figure out best way
       to do this now.
    2. Ability to get recorded or typical data automatically from some web
       service if the user provides coordinates or WMO station number.
    3. Ability to download climate data the same way.
    4. Ability to learn a transfer function between urban and rural stations
       for quick-and-not-too-dirty estimates of urban heat island effects.
       So long as time series of a year-ish are available from the locations
       of interest, the transfer _should_ be able to proceed without requiring
       too much information about the local urban canopy.
"""

# These are the sources of measured data I am familiar with. So far, the
# get_weather script can only read the following:
# 1. NCDC, 2. NSRDB, 3. MS (MeteoSuisse)
# sources = ('ncdc', 'nsrdb', 'nrel_indiasolar', 'ms', 'WY2', 'nasa_saudi')

import os
import sys
import pickle

# import pandas as pd

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

from statsmodels.tsa.arima_model import ARIMAResults


def indra(train, stcode='gen', n_sample=10, method='arma',
          fpath_in='./che_gen_iwec.a',
          ftype='espr',
          outpath='.',
          l_start=int(0), l_end=int(31*24),
          l_step=int(4*24), histlim=int(14*24),
          stlat=0.0, stlong=0.0, stalt=0.0,
          randseed=8760):

    # ------------------
    # Some initialisation house work.

    # The learning/sampling functions rely on random sampling. For one run,
    # the random seed is constant/immutable; changing it during a run
    # would not make sense. This makes the runs repeatable -- keep track of
    # the seed and you can reproduce exactly the same random number draws
    # as before.
    setseed(randseed)

    # Convert incoming stcode to lowercase.
    stcode = stcode.lower()

    # Store everything in a folder named <stcode>.
    outpath = stcode
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    # These will be the files where the outputs will be stored.
    path_model_save = os.path.join(
            outpath, 'model_{0}_{1}.p'.format(stcode, randseed))
    path_ffit_save = os.path.join(
            outpath, 'ffit_{0}_{1}.p'.format(stcode, randseed))

    print('Storing everything in folder {0}\r\n'.format(outpath))

    # Save output time series.
    picklepath = os.path.join(
            outpath, 'syn_{0}_{1}.npy'.format(stcode, randseed))

    # ----------------

    # Load data about the cities. This isn't usually necessary, so it
    # doesn't matter if this file isn't present.
    #    try:
    #        citytab = pd.read_csv(os.path.join('CityData.csv'),
    #                              dtype=dict(WMO=str, StCode=str))
    #    except:
    #        citytab = None
    #        print("I could not find CityData.csv, continuing without...")

    # See accompanying script "wfileio".
    try:
        xy_train, locdata, header = wf.get_weather(
                stcode, fpath_in, ftype, outpath=outpath)

        # The GP method needs day of year rather than day of month.
        if method == 'gp':
            temp = wf.day_of_year(xy_train[:, 1], xy_train[:, 2])
            xy_train[:, 2] = temp
            del temp

        print('Successfully retrieved weather data.\r\n')

    except Exception as err:
        print('Error: ' + str(err))
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        print('I could not read the incoming weather file. ' +
              'Terminating this run.\r\n')
        return 0

    if train:

        # Train the models.

        if method == 'gp':

            gp_list, mtrack, scaler = learngp(
                    l_start, l_end, l_step, histlim, xy_train)

            # Save gp_list and month_tracker to pickle file.
            gp_save = dict(gp_list=gp_list, mtrack=mtrack,
                           scaler=scaler, xy_train=xy_train)

            with open(path_model_save, 'wb') as fp:
                pickle.dump(gp_save, fp)

        elif method == 'arma':

            # Call resampling with null selmdl and ffit, since those
            # haven't been trained yet.
            ffit, selmdl, _ = resampling(
                    stcode, xy_train, selmdl=None, ffit=None,
                    train=True, sample=False, n_sample=n_sample,
                    picklepath=picklepath)

            selmdl.save(path_model_save)

            with open(path_ffit_save, 'wb') as fp:
                pickle.dump(ffit, fp)

    else:

        # Load models from file.

        if method == 'gp':

            with open(path_model_save, 'rb') as fp:
                gp_save = pickle.load(fp)

            gp_list = gp_save['gp_list']
            mtrack = gp_save['mtrack']
            scaler = gp_save['scaler']
            xy_train = gp_save['xy_train']

        elif method == 'arma':

            with open(path_ffit_save, 'rb') as fp:
                ffit = pickle.load(fp)

            selmdl = ARIMAResults.load(path_model_save)

    # %%

    # Call the sampling function.

    # The output, xout, is a numpy nd-array with the standard
    # columns ('month', 'day', 'hour', 'tdb', 'tdp', 'rh',
    # 'ghi', 'dni', 'dhi', 'wspd', 'wdr')

    # In this MC framework, the 'year' of weather data is meaningless.
    # When the climate change models will be added, these years will
    # mean something. For now, just add '2017' to every file.

    if method == 'arma':

        ffit, selmdl, _ = resampling(
                stcode, xy_train, selmdl, ffit, train=False,
                sample=True, picklepath=picklepath)

    elif method == 'gp':

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
                    s_shift=0, outpath=outpath, masterfile=fpath_in)
