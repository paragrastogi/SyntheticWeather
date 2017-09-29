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
#sources = ('ncdc', 'nsrdb', 'nrel_indiasolar',
#           'ms', 'WY2', 'nasa_saudi')


def indra(train, stcode='gla', n_sample=10,
          fpath_in='/usr/esru/esp-r/climate/che_gen_iwec.a',
          ftype='espr',
          outpath='.',
          l_start=int(0), l_end=int(31*24),
          l_step=int(4*24), histlim=int(14*24),
          stlat=0.0, stlong=0.0, stalt=0.0,
          randseed=8760):

    import os
    import sys
    import random

    import numpy as np
    import pandas as pd
    import pickle

    # These custom functions load and clean recorded data.
    # For now, we are only concerned with ncdc and nsrdb.
    import get_weather as gw

    # Custom functions to calculate error metrics.
    # import losses.
    # from losses import rmseloss
    # from losses import maeloss

    # import learn and sample functions from gp_funcs.
    from gp_funcs import learngp
    from gp_funcs import samplegp

    # ------------------
    # Some initialisation house work.

    # Seed random number generators.
    np.random.seed(randseed)
    random.seed = randseed

    # Convert incoming stcode to lowercase.
    stcode = stcode.lower()

    # Store everything in a folder named <stcode>.
    outpath = stcode
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    # These will be the files where the outputs will be stored.
    path_file_list = os.path.join(
            outpath, 'model_{0}_{1}.p'.format(stcode, randseed))
    path_other_file = os.path.join(
            outpath, 'm_tracker_{0}_{1}.p'.format(stcode, randseed))
    path_seed_data = os.path.join(
            outpath, 'seed_data_{0}_{1}.p'.format(stcode, randseed))

    # Oth
#    if not os.path.isdir(os.path.join(outpath, 'pickled_data')):
#        os.makedirs(os.path.join(outpath, 'pickled_data'))
#    if not os.path.isdir(os.path.join(
#            outpath, 'csv_data')):
#        os.makedirs(os.path.join(outpath, 'csv_data'))

    print('Storing everything in folder {0}\r\n'.format(outpath))

    # ----------------

    # Load data about the cities. This isn't usually necessary, so it
    # doesn't matter if this file isn't present.
    #    try:
    #        citytab = pd.read_csv(os.path.join('CityData.csv'),
    #                              dtype=dict(WMO=str, StCode=str))
    #    except:
    #        citytab = None
    #        print("I could not find CityData.csv, continuing without...")

    # See accompanying script "gw".
    try:
        ts_in, locdata = gw.get_weather(stcode, fpath_in, ftype,
                                        outpath=outpath)

        temp = gw.day_of_year(ts_in[:, 0], ts_in[:, 1])
        ts_in[:, 1] = temp
        del temp

        print('Successfully retrieved weather data.\r\n')

    except Exception as err:
        print('Error: ' + str(err))
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        print('I could not read the incoming weather file. ' +
              'Terminating this run.')
        return 0

    # Extract the relevant data from the loaded 'starter' data.
    # This step won't be needed if only the relevant variables are
    # in the input matrix. Conversely, the user might have
    # to identify the relevant time series.

    if train:

        # Train the models and store them on drive.

        gp_list, month_tracker = learngp(l_start, l_end, l_step,
                                         histlim, ts_in)

        # Save gp_list and month_tracker to pickle file.
        with open(path_file_list, 'wb') as fp:
            pickle.dump(gp_list, fp)

        with open(path_other_file, 'wb') as fp:
            pickle.dump(month_tracker, fp)

        with open(path_seed_data, 'wb') as fp:
            pickle.dump(ts_in, fp)

    else:

        # Load models from file.

        with open(path_file_list, 'rb') as fp:
            gp_list = pickle.load(fp)
        with open(path_other_file, 'rb') as fp:
            month_tracker = pickle.load(fp)
        with open(path_seed_data, 'rb') as fp:
            ts_in = pickle.load(fp)

    # %%

    picklepath = os.path.join(
            outpath, 'syn_{0}_{1}.p'.format(stcode, randseed))

    xout_un, column_names = samplegp(
            gp_list, l_start, l_end, l_step, histlim, n_sample,
            ts_in, month_tracker)

    # Save the outputs as a pickle.
    pd.to_pickle(xout_un, picklepath)


    # Save synthetic time series.
    gw.

    # The generation part ends here - the rest is just plotting various things.
