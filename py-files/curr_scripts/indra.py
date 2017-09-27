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


def indra(seedfile, stcode='gla', n_samples=10,
          path_wthr='/usr/esru/esp-r/climate',
          outpath='.',
          l_start=int(0), l_end=int(31*24),
          l_step=int(4*24), histlim=int(14*24),
          stlat=0.0, stlong=0.0, stalt=0.0,
          randomseed=8760):

    import os
    import random

    import numpy as np
    import pandas as pd
    import pickle

    # These custom functions load and clean recorded data.
    # For now, we are only concerned with ncdc and nsrdb.
    import get_weather_data as gw

    # Custom functions to calculate error metrics.
    # import losses.
    # from losses import rmseloss
    # from losses import maeloss

    # import les petites fonctionnes utiles.
    from petites import learngp
    from petites import samplegp

    # ------------------
    # Some initialisation house work.

    # Seed random number generators.
    np.random.seed(randomseed)
    random.seed = randomseed

    # Convert incoming stcode to lowercase.
    stcode = stcode.lower()

    # Store everything in a folder named <stcode>.
    outpath = stcode
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    path_file_list = os.path.join(outpath,
                                  'gp_list_{0}.p'.format(stcode))
    path_other_file = os.path.join(outpath,
                                   'm_tracker_{0}.p'.format(stcode))
    path_seed_data = os.path.join(outpath,
                                  'seed_data_{0}.p'.format(stcode))

    if not os.path.isdir(os.path.join(outpath, 'pickled_data')):
        os.makedirs(os.path.join(outpath, 'pickled_data'))
    if not os.path.isdir(os.path.join(
            outpath, 'csv_collated_data')):
        os.makedirs(os.path.join(outpath, 'csv_collated_data'))

    print('Storing everything in folder {0}\r'.format(outpath))

    # ----------------

    # The program needs to know the mode you are running in.
    # There are two modes supported:
    # 1. With starter data - seedfile = True
    # 2. Without started data - seedfile = False

    # print('cwd = %s' % os.getcwd())

    if seedfile:

        # The script needs some 'starter' data.

        currentdata = None

        try:
            # See accompanying script "gw".
            currentdata, _ = gw.read_epw(
                    path_wthr, stcode, force=True, outpath=outpath)
        except:
            print('Could not find file {0} or could not read ' \
                  'the file I found.'.format(path_wthr))

        # If there is no path to a weather file, then this bit of code loads
        # data from pre-processed pickle files or processes CSV files.

        if (currentdata is None):

            # Load data about the cities. This is just for this example.
            try:
                citytab = pd.read_csv(os.path.join('CityData.csv'),
                                      dtype=dict(WMO=str, StCode=str))
            except:
                citytab = None
                print("I could not find CityData.csv, continuing without...")

            # Specify the sources of the actual data -
            # please follow AMY keywords list.
            sources = ('ncdc', 'nsrdb', 'nrel_indiasolar',
                       'ms', 'WY2', 'nasa_saudi')

            # See accompanying script "gw".
            currentdata, _ = gw.get_weather(
                    stcode, sources, path_wthr,
                    citytab=citytab, outpath=outpath)

        # T = currentdata.shape[0]

        # Extract the relevant data from the loaded 'starter' data.
        # This step won't be needed if only the relevant variables are
        # in the input matrix. Conversely, the user might have
        # to identify the relevant time series.
        ts_curr_in = (np.vstack([
                currentdata.month.values,
                currentdata.index.dayofyear, currentdata.hour.values,
                currentdata.tdb.values, currentdata.tdp.values,
                currentdata.rh.values, currentdata.ghi.values,
                currentdata.dni.values, currentdata.dhi.values,
                currentdata.wspd.values, currentdata.wdr.values])).T

        # ts_hist_norm = scaler.transform(ts_hist_in)

        # ts_hist_in = (np.vstack([
        #        historicaldata.index.dayofyear, historicaldata.hour.values,
        #        historicaldata.tdb.values, historicaldata.tdp.values,
        #        historicaldata.rh.values, historicaldata.ghi.values,
        #        historicaldata.dni.values, historicaldata.dhi.values,
        #        historicaldata.wspd.values, historicaldata.wdr.values])).T

        #
        # %%

        # This is the part of the generator where GP models are fitted to
        # the incoming data. Naturally, it is only invoked if the function
        # was invoked with starter data.

        gp_list, month_tracker = learngp(l_start, l_end, l_step,
                                     histlim, ts_curr_in)

        # Save gp_list and month_tracker to pickle file.
        with open(path_file_list, 'wb') as fp:
            pickle.dump(gp_list, fp)

        with open(path_other_file, 'wb') as fp:
            pickle.dump(month_tracker, fp)

        with open(path_seed_data, 'wb') as fp:
            pickle.dump(ts_curr_in, fp)

        # end seedfile if statement.

    else:

        with open(path_file_list, 'rb') as fp:
            gp_list = pickle.load(fp)
        with open(path_other_file, 'rb') as fp:
            month_tracker = pickle.load(fp)
        with open(path_seed_data, 'rb') as fp:
            ts_curr_in = pickle.load(fp)

    # %%

    picklepath = os.path.join(outpath, 'syn_{0}.p'.format(stcode))

    xout_un, column_names = samplegp(
            gp_list, l_start, l_end, l_step, histlim, n_samples,
            ts_curr_in, month_tracker, picklepath, outpath)
    # Save synthetic time series.
    for n in range(0, n_samples):
        filepath = os.path.join(
                outpath, 'syn_{0}_{1}.csv'.format(stcode, n))
        np.savetxt(filepath, np.squeeze(xout_un[:, :, n]), '%6.2f',
                   delimiter=',', header=','.join(column_names))

    # The generation part ends here - the rest is just plotting various things.
