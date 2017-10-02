import os
import numpy as np
import pandas as pd
import csv
from scipy import interpolate

import petites

"""
This file contains functions to:
    1. load weather data from 'typical' and 'actual' (recorded) weather
       data files.
    2. Write out synthetic weather data to EPW or ESPr weather file formats.
    3. Associated helper functions.
"""

__author__ = 'Parag Rastogi'

# ISSUES TO ADDRESS
# 1. Harmonize WMO numbers - if the incoming number is 5 digits,
# add leading zero (e.g., Geneva)
# 2. Implement something to convert GHI to DNI and DHI.
# Maybe use Erbs model like before.

# %%


def get_weather(stcode, fpath, ftype='espr', outpath='xxx'):

    # This function calls the relevant reader based on the ftype.

    # Initialise as a non-object.
    wdata = None

    locdata = None

    if os.path.isfile(fpath):
        print('Running weather file reader for station ' +
              '{0}. Expecting format {1}.\r\n'.format(stcode, ftype))
    else:
        print('I cannot find file {0}.'.format(fpath) +
              ' Returning empty dataframe.\r\n')
        return wdata, locdata

    # Load data for given station.

    if ftype == 'pickle':
        try:
            ts_in = pd.read_pickle(fpath)
            return ts_in
        except Exception as err:
            print('You asked me to read a pickle but I could not. ' +
                  'Trying all other formats.\r\n')
            print('Error: ' + str(err))
            wdata = None

    elif ftype == 'epw':

        try:
            wdata, locdata, header = read_epw(fpath)
        except Exception as err:
            print('Error: ' + str(err))
            wdata = None

    elif ftype == 'espr':

        try:
            wdata, locdata, header = read_espr(fpath)
        except Exception as err:
            print('Error: ' + str(err))
            wdata = None

    elif ftype == 'csv':

        try:
            wdata = pd.read_csv(fpath, header=0)
            wdata.columns = ['month', 'day', 'hour', 'tdb', 'tdp', 'rh',
                             'ghi', 'dni', 'dhi', 'wspd', 'wdr']
            # Location data is nonsensical, except for station code,
            # which will be reassigned later in this function.
            locdata = dict(loc=stcode, lat='00', long='00',
                           tz='00', alt='00', wmo='000000')
            header = '# Unknown incoming file format ' + \
                     '(not epw or espr)\r\n' + \
                     '# Dummy location data: ' + \
                     'loc: {0}'.format(locdata['loc']) + \
                     'lat: {0}'.format(locdata['lat']) + \
                     'long: {0}'.format(locdata['long']) + \
                     'tz: {0}'.format(locdata['tz']) + \
                     'alt: {0}'.format(locdata['alt']) + \
                     'wmo: {0}'.format(locdata['wmo']) + \
                     '\r\n'

        except Exception as err:
            print('Error: ' + str(err))
            wdata = None

    # End ftype if statement.

    # The first try didn't work for some reason.

    if wdata is None:
        print('I could not read the file you gave me with the format ' +
              'you specified. Trying all readers.\r\n')

        # Once more unto the breach.

        for fmt in wformats:
            if fmt == 'epw':
                try:
                    wdata, locdata, header = read_epw(fpath)
                except Exception as err:
                    print('Error: ' + str(err))
                    wdata = None
            elif fmt == 'espr':
                try:
                    wdata, locdata, header = read_espr(fpath)
                except Exception as err:
                    print('Error: ' + str(err))
                    wdata = None
            elif fmt == 'csv':
                try:
                    wdata = pd.read_csv(fpath, header=0)
                    wdata.columns = ['month', 'day', 'hour', 'tdb', 'tdp',
                                     'rh', 'ghi', 'dni', 'dhi', 'wspd', 'wdr']
                    locdata = dict(loc='xxx', lat='00', long='00',
                                   tz='00', alt='00', wmo='000000')
                    header = '# Unknown incoming file format ' + \
                     '(not epw or espr)\r\n' + \
                     '# Dummy location data: ' + \
                     'loc: {0}'.format(locdata['loc']) + \
                     'lat: {0}'.format(locdata['lat']) + \
                     'long: {0}'.format(locdata['long']) + \
                     'tz: {0}'.format(locdata['tz']) + \
                     'alt: {0}'.format(locdata['alt']) + \
                     'wmo: {0}'.format(locdata['wmo']) + \
                     '\r\n'

                except Exception as err:
                    print('Error: ' + str(err))
                    wdata = None

            if wdata is not None:
                break

    # End wformats for loop and if wdata is None.

    locdata['loc'] = stcode

    if wdata is None:
        print('All attempts to read your file were unsuccesful, ' +
              'returning empty table.')
        return wdata, locdata

    else:

        ts_in = (np.vstack([
                wdata.month.values,
                wdata.day.values, wdata.hour.values,
                wdata.tdb.values, wdata.tdp.values,
                wdata.rh.values, wdata.ghi.values,
                wdata.dni.values, wdata.dhi.values,
                wdata.wspd.values, wdata.wdr.values])).T
        # Note that the second column is day_of_month here but in the main
        # indra script it will be converted to day_of_year.

        return ts_in, locdata, header

    # Ignore this bit of code for now.

    # Load actual data for given station
    # force = False
    # dataout = read_others(stcode, force, sources, outpath=outpath)
    # Always use NCDC in addition to the country-specific weather source
    # (like meteosuisse), since that can help fill data.


# %%

# Number of days in each month.
m_days = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)


def day_of_year(month, day):

    month = month.astype(int) - 1
    doy = np.zeros_like(day, dtype=int)

    for m, mon in enumerate(month):
        doy[m] = (day[m] + np.sum(m_days[0:mon])).astype(int)

    return doy

# End function day_of_year


def day_of_month(day):

    month = np.zeros_like(day, dtype=int)
    dom = np.zeros_like(day, dtype=int)

    for d, doy in enumerate(day):

        rem = doy
        prev_ndays = 0

        for m, ndays in enumerate(m_days):
            # The iterator 'm' starts at zero.

            if rem <= 0:
                # Iterator has now reached the incomplete month.

                # The previous month is the correct month.
                # Not subtracting 1 because the iterator starts at 0.
                month[d] = m

                # Add the negative remainder to the previous month's days.
                dom[d] = rem + prev_ndays
                break

            # Subtract number of days in this month from the day.
            rem -= ndays
            # Store number of days from previous month.
            prev_ndays = ndays

    return month, dom

# End function day_of_month
# %%


def read_epw(fpath='./gen_iwec.epw'):

    # Uniform date index for all tmy weather data tables.
    dates = pd.date_range('1/1/2017', periods=8760, freq='H')

    # Convert the names to lowercase.
    epw_colnames = [x.lower() for x in epw_colnames]

    # List of EPW files.
    tmy_filelist = []

    # Look in folder for weather data files.
    if os.path.isdir(fpath):
        filepaths = [os.path.join(fpath, f.lower())
                     for f in (os.listdir(fpath))]
    else:
        filepaths = fpath

    if isinstance(filepaths, (list, tuple)):
        for f in filepaths:
            # Avoid AMY files.
            if 'epw' in f \
                and any(s in f.lower() for s in keywords['tmy']) \
                    and not any(s in f.lower() for s in keywords['amy']):
                tmy_filelist.append(f)
    else:
        tmy_filelist = [filepaths]

    didx = 0
    typdata = pd.DataFrame()

    for f in tmy_filelist:

        # Read table, ignoring header lines.
        wdata_typ = pd.read_csv(f, delimiter=',', skiprows=8, header=None,
                                names=epw_colnames)
        wdata_typ.index = dates

        if len(wdata_typ.columns) == 35:
            # Some files have three extra columns
            # (usually the TMY files from USDOE).
            # Delete those columns if found.
            wdata_typ = wdata_typ.drop(['unknownvar1', 'unknownvar2',
                                        'unknownvar3'], axis=1)

        # Read header and assign all metadata.
        with open(f, 'r') as epwfile:
            read_header = csv.reader(epwfile, delimiter=',')
            header = next(read_header)

        locdata = dict(loc=header[1], lat=header[6], long=header[7],
                       tz=header[8], alt=header[9], wmo=header[5])

        # Assign header information to table.
        #        wdata_typ = wdata_typ.assign(
        #                latitude=latitude, longitude=longitude,
        #                altitude=altitude, wmo=wmo, tz=tz,
        #                location=location, loccode=loccode)
        if didx == 0:
            typdata = wdata_typ
        else:
            typdata = typdata.append(wdata_typ)

        didx += 1

    if typdata.empty:

        print('Could not locate a file with given station name.' +
              ' Returning empty table.\r\n')

    return typdata, locdata, header

# ----------- END read_epw function -----------


def read_espr(fpath='./che_geneva.iwec.a'):

    # Missing functionality - reject call if path points to binary file.
    
    # Uniform date index for all tmy weather data tables.
    dates = pd.date_range('1/1/2017', periods=8760, freq='H')

    fpath_fldr, fpath_name = os.path.split(fpath)
    sitename = fpath_name.split(sep='.')
    sitename = sitename[0]

    with open(fpath) as f:
        content = f.readlines()

    content = [x.strip() for x in content]

    # Split the contents into a header and body.
    header = content[0:13]

    locline = [line for line in header if ('latitude' in line)][0].split()
    siteline = [line for line in header if ('site name' in line)][0].split()
    locdata = dict(loc=siteline[0], lat=locline[1], long=locline[2],
                   tz='00', alt='0000', wmo='000000')
    # ESP-r files do not contain timezone, altitude, or WMO number.

    body = content[12:]

    del content

    # Find the lines with day tags.
    daylines = [[idx, line] for [idx, line] in enumerate(body)
                if 'day' in line]

    dataout = np.zeros([8760, 11])

    dcount = 0

    for idx, day in daylines:

        # Get the next 24 lines.
        daylist = np.asarray(body[idx+1:idx+25])

        # Split each line of the current daylist into separate strings.
        splitlist = [element.split() for element in daylist]

        # Convert each element to a integer, then convert the resulting
        # list to a numpy array.
        daydata = np.asarray([list(map(int, x)) for x in splitlist])

        # Today's time slice.
        dayslice = range(dcount, dcount+24, 1)

        # This will split the day-month header line on the gaps.
        splitday = day.split()

        # Month.
        dataout[dayslice, 0] = np.repeat(int(splitday[-1]), len(dayslice))

        # Day of month.
        dataout[dayslice, 1] = np.repeat(int(splitday[2]), len(dayslice))

        # Hour (of day).
        dataout[dayslice, 2] = np.arange(0, 24, 1)

        # tdb, input is in deci-degrees, convert to degrees.
        dataout[dayslice, 3] = daydata[:, 1]/10

        # rh, in percent.
        dataout[dayslice, 5] = daydata[:, 5]

        # ghi, in W/m2.
        dataout[dayslice, 7] = daydata[:, 0]

        # dni, in W/m2.
        dataout[dayslice, 8] = daydata[:, 2]

        # wdr, input is in deci-m/s.
        dataout[dayslice, 9] = daydata[:, 3]/10

        # wspd, clockwise deg from north.
        dataout[dayslice, 10] = daydata[:, 4]

        dcount += 24

    # tdp, calculated from tdb and rh.
    dataout[:, 4] = petites.tdb2tdp(dataout[:, 3], dataout[:, 5])

    # wspd can have bogus values (999)
    dataout[dataout[:, 10] >= 999., 10] = np.nan
    idx = np.arange(0, dataout.shape[0])
    duds = np.logical_or(np.isinf(dataout[:, 10]), np.isnan(dataout[:, 10]))
    int_func = interpolate.interp1d(
            idx[np.logical_not(duds)], dataout[np.logical_not(duds), 10],
            kind='nearest', fill_value='extrapolate')
    dataout[duds, 10] = int_func(idx[duds])

    clmdata = pd.DataFrame(data=dataout, index=dates,
                           columns=['month', 'day', 'hour',
                                    'tdb', 'tdp', 'rh',
                                    'ghi', 'dni', 'dhi', 'wspd', 'wdr'])

    return clmdata, locdata, header

# ----------- END read_espr function -----------


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


def give_weather(ts, locdata, stcode, header, ftype, outpath='.'):
    
    print(header)
    print(ts)

    n_sample, rem = divmod(ts_out.shape[0], 8760)

    if rem != 0:
        print('Number of rows of output table given to give_weather ' +
              'is not a multiple of 8760. You might want to see why ' +
              'that is the case.')

        if ftype == 'espr':
            for n in range(0, n_sample):
                filepath = os.path.join(
                        outpath, 'syn_{0}_{1}.a'.format(stcode, n))

                np.savetxt(filepath, np.squeeze(ts[:, :, n]), '%6.2f',
                   delimiter=',', header=header, comments='')

                if os.path.isfile(filepath):
                    success = True
                else:
                    success = False

        elif ftype == 'epw':

            for n in range(0, n_sample):
                filepath = os.path.join(
                        outpath, 'syn_{0}_{1}.epw'.format(stcode, n))

                np.savetxt(filepath, np.squeeze(ts[:, :, n]), '%6.2f',
                   delimiter=',',
                   header=header )
    
                if os.path.isfile(filepath):
                    success = True
                else:
                    success = False            
            
        else:

            for n in range(0, n_sample):
                filepath = os.path.join(
                        outpath, 'syn_{0}_{1}.csv'.format(stcode, n))

                np.savetxt(filepath, np.squeeze(ts[:, :, n]), '%6.2f',
                   delimiter=',',
                   header=header + ','.join(std_cols))
    
                if os.path.isfile(filepath):
                    success = True
                else:
                    success = False

    return n_files

# ----------- End give_weather function. -----------


#def write_epw(ts, header):
#
#    np.sav
#    
## ----------- End write_csv function. -----------
#
#
#def write_csv(ts, header):
#
#    np.savetxt(filepath, np.squeeze(ts[:, :, n]), '%6.2f',
#               delimiter=',',
#               header=header + ','.join(std_cols))
#
#    if os.path.isfile(filepath):
#        success = True
#    else:
#        success = False
#
#    return success

# ----------- End write_csv function. -----------


# Useful strings and constants.

# List of keywords that identify TMY and AMY files.
keywords = dict(tmy=('nrel', 'iwec', 'ishrae', 'cwec',
                     'igdg', 'tmy3', 'meteonorm'),
                amy=('ncdc', 'nsrdb', 'nrel_indiasolar',
                     'ms', 'WY2', 'nasa_saudi'))
wformats = ('epw', 'espr', 'csv')

# List of values that could be NaNs.
nanlist = ('9900', '-9900', '9999', '99', '-99', '9999.9', '999.9', ' ', '-')

# Names of the columns in EPW files. Usually ignore the last
# three columns.
epw_colnames = ['Year', 'Month', 'Day', 'Hour', 'Minute', 'QualFlags',
                'TDB', 'TDP', 'RH', 'ATMPR', 'ETRH', 'ETRN', 'HIR',
                'GHI', 'DNI', 'DHI', 'GHE', 'DNE', 'DHE', 'ZL',
                'WDR', 'WSPD', 'TSKY', 'OSKY', 'VIS', 'CHGT',
                'PWO', 'PWC', 'PWT', 'AOPT', 'SDPT',
                'SLAST', 'UnknownVar1', 'UnknownVar2', 'UnknownVar3']

# A generic ESPR header that can be used in a print or str
# command with format specifiers.
espr_generic_header = """*CLIMATE
# ascii weather file from {0},
# defined in: {1}
# col 1: Diffuse solar on the horizontal (W/m**2)
# col 2: External dry bulb temperature   (Tenths DEG.C)
# col 3: Direct normal solar intensity   (W/m**2)
# col 4: Prevailing wind speed           (Tenths m/s)
# col 5: Wind direction     (clockwise deg from north)
# col 6: Relative humidity               (Percent)
{2}               # site name
 {3},{4},{5},{6}   # year, latitude, long diff, direct normal rad flag
 {7},{8}    # period (julian days)"""

# The standard columns used by indra.
std_cols = ('month', 'day', 'hour', 'tdb', 'tdp', 'rh',
            'ghi', 'dni', 'dhi', 'wspd', 'wdr')
