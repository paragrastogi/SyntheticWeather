# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 15:51:51 2017

@author: rasto
"""
import os
import pandas as pd
import numpy as np

# List of values that could be nans.
nanlist = ('9900', '-9900', '9999', '99', '-99', '9999.9', '999.9', ' ', '-')
# List of keywords that identify TMY and AMY files.
keywords = dict(tmy=('nrel', 'iwec', 'ishrae', 'cwec',
                     'igdg', 'tmy3', 'meteonorm'),
                amy=('ncdc', 'nsrdb', 'nrel_indiasolar',
                     'ms', 'WY2', 'nasa_saudi'))


def read_others(stcode, sources, outpath='xxx',
                path_actual_top=os.path.join(
                        '..', '..', 'WeatherData', 'HistoricalData')):

    # Create an empty data frame to be appended to later.
    dataout = pd.DataFrame()

    if 'ncdc' in sources:
        path_ncdc = os.path.join(path_actual_top, 'ncdc')
        ncdc_filelist = [os.path.join(path_ncdc, o)
                         for o in os.listdir(path_ncdc)
                         if not os.path.isdir(os.path.join(path_ncdc, o))]
        ncdc_colnames = ['date', 'time', 'wdr', 'wspd',
                         'tdb', 'tdp', 'atmpr', 'rh']

        wdata_ncdc = pd.DataFrame()

        for f in ncdc_filelist:
            if stcode is 'all':
                print('I cannot load NCDC data for all stations -- ' +
                      'the resulting dataframe would be too big.\r\n')
                break

            if ('ncdc' in f) and (stcode in f):
                print('Processing station {0}.\r\n'.format(stcode))
                wdata_ncdc = pd.read_csv(f, delimiter=',', skiprows=2,
                                         header=None,
                                         usecols=[2, 3, 7, 10,
                                                  12, 14, 16, 18],
                                         na_filter=True,
                                         na_values=nanlist,
                                         names=ncdc_colnames,
                                         dtype=dict(date=str,
                                                    time=str))
                dater = wdata_ncdc['date'].str.extract(
                    '(\d{4})(\d{2})(\d{2})', expand=False)
                dater.columns = ['year', 'month', 'day']
                timer = wdata_ncdc['time'].str.extract(
                    '(\d{2})(\d{2})', expand=False)
                timer.columns = ['hour', 'minute']
                wdata_ncdc = pd.concat([dater, timer, wdata_ncdc],
                                       axis=1, join='inner')
                wdata_ncdc = wdata_ncdc.drop(['date', 'time'], axis=1)
                wdata_ncdc.index = pd.to_datetime(wdata_ncdc[[
                    'year', 'month', 'day', 'hour', 'minute']])
                break
            else:
                continue

        # The dataframe wdata_ncdc is renamed to dataout.
        # In case only NCDC data is found, actual data
        # returned will consist only of NCDC data.
        dataout = dataout.append(wdata_ncdc)
        del wdata_ncdc

    if 'nsrdb' in sources:
        path_nsrdb = os.path.join(path_actual_top, 'nsrdb')
        nsrdb_filelist = [os.path.join(path_nsrdb, o)
                          for o in os.listdir(path_nsrdb)
                          if not os.path.isdir(os.path.join(path_nsrdb, o))]
        nsrdb_colnames = ['date', 'time', 'GHI_metstat',
                          'DNI_metstat', 'DHI_metstat',
                          'GHI_suny', 'DNI_suny', 'DHI_suny',
                          'GHI_measure', 'DNI_measure', 'DHI_measure']

        # Create an empty data frame to be appended to later.
        wdata_nsrdb = pd.DataFrame()

        for f in nsrdb_filelist:
            if stcode is 'all':
                print('I cannot load NSRDB data for all stations -- ' +
                      'the resulting dataframe would be too big.\r\n')
                dataout = dataout.append(wdata_nsrdb)
                break

            if all(s in f for s in ('nsrdb', 'csv', stcode)):
                wdata_nsrdb = pd.read_csv(f, delimiter=',', skiprows=1,
                                          header=None,
                                          usecols=[0, 1, 6, 9, 12, 15,
                                                   17, 19, 27, 29, 31],
                                          na_filter=True,
                                          names=nsrdb_colnames,
                                          na_values=nanlist,
                                          dtype=dict(date=str, time=str))
                dater = wdata_nsrdb['date'].str.extract(
                    '(\d{4})-(\d{2})-(\d{2})', expand=False)
                dater.columns = ['year', 'month', 'day']
                timer = wdata_nsrdb['time'].str.extract(
                    '(\d{1,2}):(\d{2})', expand=False)
                timer.columns = ['hour', 'minute']
                wdata_nsrdb = pd.concat([dater, timer, wdata_nsrdb],
                                        axis=1, join='inner')
                wdata_nsrdb = wdata_nsrdb.drop(['date', 'time'], axis=1)
                wdata_nsrdb.index = pd.to_datetime(
                    wdata_nsrdb[['year', 'month', 'day', 'hour', 'minute']])

                # Little bit of code to pick best values from the three
                # solar radiation types available in the NSRDB files.
                wdata_nsrdb = nsrdb_solar_clean(wdata_nsrdb)
            else:
                continue

            dataout = dataout.append(wdata_nsrdb)

    if 'meteosuisse' in sources:
        path_ms = os.path.join(path_actual_top, 'ms')
        ms_filelist = [os.path.join(path_ms, o) for o in os.listdir(path_ms)
                       if not os.path.isdir(os.path.join(path_ms, o))]
        ms_colnames = ['datetime', 'ghi', 'atmpr', 'tdb',
                       'rh', 'tdp', 'wspd', 'wdr']

        # Create an empty data frame to be appended to later.
        wdata_ms = pd.DataFrame()

        for f in ms_filelist:
            if stcode is 'all':
                print('I cannot load MeteoSuisse data for all stations -- ' +
                      'the resulting dataframe would be too big.\r\n')
                dataout = dataout.append(wdata_ms)
                break

            if all(s in f for s in ('meteosuisse', 'txt', stcode)):
                wdata_ms = pd.read_csv(f, delimiter=';', skiprows=1,
                                       header=None,
                                       usecols=[1, 2, 3, 4, 5, 6, 7, 8],
                                       na_filter=True, na_values=nanlist,
                                       names=ms_colnames,
                                       dtype=dict(datetime=str,
                                                  ghi=np.float64,
                                                  atmpr=np.float64,
                                                  tdb=np.float64,
                                                  rh=np.float64,
                                                  tdp=np.float64,
                                                  wspd=np.float64,
                                                  wdr=np.float64))
                dater = wdata_ms['datetime'].str.extract(
                    '(\d{4})(\d{2})(\d{2})(\d{2})', expand=False)
                dater['minute'] = (np.zeros(dater.shape[0])).tolist()
                dater.columns = ['year', 'month', 'day', 'hour', 'minute']

                wdata_ms = pd.concat([dater, wdata_ms], axis=1, join='inner')
                wdata_ms = wdata_ms.drop(['datetime'], axis=1)
                wdata_ms.index = pd.to_datetime(
                    wdata_ms[['year', 'month', 'day', 'hour', 'minute']])

                # Put in dummies for DNI, DHI and WMO.
                wdata_ms['dni'] = (np.zeros(wdata_ms.shape[0])).tolist()
                wdata_ms['dhi'] = (np.zeros(wdata_ms.shape[0])).tolist()
                # wdata_ms['wmo'] = (np.zeros(wdata_ms.shape[0])).tolist()

            else:
                continue

            dataout = dataout.append(wdata_ms)

    # If blank, output blank.
    if dataout.empty:
        print('Could not locate any files with given station name. ' +
              'Returning empty table.\r\n')
        return dataout
    else:
        # delete  the minute column
        del dataout['minute']
        # Reorganise columns and return.
        dataout = dataout[['year', 'month', 'day', 'hour',
                           'tdb', 'tdp', 'rh', 'ghi', 'dni', 'dhi',
                           'wdr', 'wspd', 'atmpr']]

        # Convert time columns to numbers.
        dataout.year = pd.to_numeric(dataout.year)
        dataout.month = pd.to_numeric(dataout.month)
        dataout.day = pd.to_numeric(dataout.day)
        dataout.hour = pd.to_numeric(dataout.hour)

        tempdata = pd.DataFrame()
        for col in list(dataout):
            tempdata[col] = dataout.groupby([
                    'year', 'month', 'day', 'hour'])[col].mean()

        # If all elements of a row are nans, then discard that row.
        tempnans = tempdata.apply(np.isnan, axis=1)
        nancheck = tempnans.all(axis=1).values

        dataout = pd.DataFrame()

        for col in list(tempdata):
            dataout[col] = tempdata[col][np.logical_not(nancheck)]

        del tempdata

        # Final sanity checks.
        notnans = dataout.tdb.values[np.logical_not(
                np.isnan(dataout.tdb.values))]
        notnanidx = np.logical_not(np.isnan(dataout.tdb.values))
        notnans[(notnans > 55) | (notnans < -55)] = np.nan

        dataout.tdb.values[notnanidx] = notnans

        notnans = dataout.tdp.values[np.logical_not(
                np.isnan(dataout.tdp.values))]
        notnanidx = np.logical_not(np.isnan(dataout.tdp.values))
        notnans[(notnans > 60) | (notnans < -60)] = np.nan

        dataout.tdp.values[notnanidx] = notnans

        # Reassing original index.
        dataout.index = pd.to_datetime(
                    dataout[['year', 'month', 'day', 'hour']])

        #        pd.to_pickle(dataout, picklepath)
        #        dataout.to_csv(csvpath, na_rep='NA')

        return dataout

# ----------- END read_others function -----------


def nsrdb_solar_clean(wdata_nsrdb):

    wdata_copy = wdata_nsrdb.copy()

    metstat = wdata_copy.GHI_metstat
    measure = wdata_copy.GHI_measure

    idx = 0
    temp = wdata_copy.GHI_suny

    wdata_copy = wdata_copy.drop([
            'GHI_suny', 'GHI_metstat', 'GHI_measure'], axis=1)

    for e in temp:
        if not np.isnan(e):
            continue
        else:
            if not np.isnan(metstat[idx]):
                temp[idx] = metstat[idx]
            else:
                if not np.isnan(measure[idx]):
                    temp[idx] = measure[idx]
                else:
                    temp[idx] = float('nan')

        idx += 1

    wdata_copy = pd.concat([temp, wdata_copy], axis=1, join='inner')
    wdata_copy = wdata_copy.rename(columns={'GHI_suny': 'ghi'})

    # Do the same for the DNI columns.

    metstat = wdata_copy.DNI_metstat
    measure = wdata_copy.DNI_measure

    idx = 0
    temp = wdata_copy.DNI_suny

    wdata_copy = wdata_copy.drop([
            'DNI_suny', 'DNI_metstat', 'DNI_measure'], axis=1)

    for e in temp:
        if not np.isnan(e):
            continue
        else:
            if not np.isnan(metstat[idx]):
                temp[idx] = metstat[idx]
            else:
                if not np.isnan(measure[idx]):
                    temp[idx] = measure[idx]
                else:
                    temp[idx] = float('nan')

        idx += 1

    wdata_copy = pd.concat([temp, wdata_copy], axis=1, join='inner')
    wdata_copy = wdata_copy.rename(columns={'DNI_suny': 'dni'})

    # And for DHI.

    metstat = wdata_copy.DHI_metstat
    measure = wdata_copy.DHI_measure

    idx = 0
    temp = wdata_copy.DHI_suny

    wdata_copy = wdata_copy.drop([
            'DHI_suny', 'DHI_metstat', 'DHI_measure'], axis=1)

    for e in temp:
        if not np.isnan(e):
            continue
        else:
            if not np.isnan(metstat[idx]):
                temp[idx] = metstat[idx]
            else:
                if not np.isnan(measure[idx]):
                    temp[idx] = measure[idx]
                else:
                    temp[idx] = float('nan')

        idx += 1

    wdata_copy = pd.concat([temp, wdata_copy], axis=1, join='inner')
    wdata_copy = wdata_copy.rename(columns={'DHI_suny': 'dhi'})

    return wdata_copy

# ----------- END nsrdb_solar_clean function. -----------