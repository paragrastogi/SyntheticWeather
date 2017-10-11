import os
import numpy as np
import pandas as pd
import csv
from scipy import interpolate

import petites

"""
This file contains functions to:
    1. load weather data from "typical" and "actual" (recorded) weather
       data files.
    2. Write out synthetic weather data to EPW or ESPr weather file formats.
    3. Associated helper functions.
"""

__author__ = "Parag Rastogi"

# ISSUES TO ADDRESS
# 1. Harmonize WMO numbers - if the incoming number is 5 digits,
# add leading zero (e.g., Geneva)
# 2. Implement something to convert GHI to DNI and DHI.
# Maybe use Erbs model like before.

# %%

# Useful strings and constants.

# List of keywords that identify TMY and AMY files.
keywords = dict(tmy=("nrel", "iwec", "ishrae", "cwec",
                     "igdg", "tmy3", "meteonorm"),
                amy=("ncdc", "nsrdb", "nrel_indiasolar",
                     "ms", "WY2", "nasa_saudi"))
wformats = ("epw", "espr", "csv")

# List of values that could be NaNs.
nanlist = ("9900", "-9900", "9999", "99", "-99", "9999.9", "999.9", " ", "-")

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
std_cols = ("year", "month", "day", "hour", "tdb", "tdp", "rh",
            "ghi", "dni", "dhi", "wspd", "wdr")

# %%


def get_weather(stcode, fpath, ftype="espr"):

    # This function calls the relevant reader based on the ftype.

    # Initialise as a non-object.
    wdata = None
    locdata = None
    header = None

    if os.path.isfile(fpath):
        print("Running weather file reader for station " +
              "{0}. Expecting format {1}.\r\n".format(stcode, ftype))
    else:
        print("I cannot find file {0}.".format(fpath) +
              " Returning empty dataframe.\r\n")
        return wdata, locdata, header

    # Load data for given station.

    if ftype == "pickle":
        try:
            wdata_array = pd.read_pickle(fpath)
            return wdata_array
        except Exception as err:
            print("You asked me to read a pickle but I could not. " +
                  "Trying all other formats.\r\n")
            print("Error: " + str(err))
            wdata = None

    elif ftype == "epw":

        try:
            wdata, locdata, header = read_epw(fpath)
        except Exception as err:
            print("Error: " + str(err))
            wdata = None

    elif ftype == "espr":

        try:
            wdata, locdata, header = read_espr(fpath)
        except Exception as err:
            print("Error: " + str(err))
            wdata = None

    elif ftype == "csv":

        try:
            wdata = pd.read_csv(fpath, header=0)
            wdata.columns = ["month", "day", "hour", "tdb", "tdp", "rh",
                             "ghi", "dni", "dhi", "wspd", "wdr"]
            # Location data is nonsensical, except for station code,
            # which will be reassigned later in this function.
            locdata = dict(loc=stcode, lat="00", long="00",
                           tz="00", alt="00", wmo="000000")
            header = "# Unknown incoming file format " + \
                     "(not epw or espr)\r\n" + \
                     "# Dummy location data: " + \
                     "loc: {0}".format(locdata["loc"]) + \
                     "lat: {0}".format(locdata["lat"]) + \
                     "long: {0}".format(locdata["long"]) + \
                     "tz: {0}".format(locdata["tz"]) + \
                     "alt: {0}".format(locdata["alt"]) + \
                     "wmo: {0}".format(locdata["wmo"]) + \
                     "\r\n"

        except Exception as err:
            print("Error: " + str(err))
            wdata = None
            header = None
            locdata = None

    # End ftype if statement.

    # The first try didn"t work for some reason.

    if wdata is None:
        print("I could not read the file you gave me with the format " +
              "you specified. Trying all readers.\r\n")

        # Once more unto the breach.

        for fmt in wformats:
            if fmt == "epw":
                try:
                    wdata, locdata, header = read_epw(fpath)
                except Exception as err:
                    print("Error: " + str(err))
                    wdata = None
            elif fmt == "espr":
                try:
                    wdata, locdata, header = read_espr(fpath)
                except Exception as err:
                    print("Error: " + str(err))
                    wdata = None
            elif fmt == "csv":
                try:
                    wdata = pd.read_csv(fpath, header=0)
                    wdata.columns = ["month", "day", "hour", "tdb", "tdp",
                                     "rh", "ghi", "dni", "dhi", "wspd", "wdr"]
                    locdata = dict(loc="xxx", lat="00", long="00",
                                   tz="00", alt="00", wmo="000000")
                    header = "# Unknown incoming file format " + \
                             "(not epw or espr)\r\n" + \
                             "# Dummy location data: " + \
                             "loc: {0}".format(locdata["loc"]) + \
                             "lat: {0}".format(locdata["lat"]) + \
                             "long: {0}".format(locdata["long"]) + \
                             "tz: {0}".format(locdata["tz"]) + \
                             "alt: {0}".format(locdata["alt"]) + \
                             "wmo: {0}".format(locdata["wmo"]) + \
                             "\r\n"

                except Exception as err:
                    print("Error: " + str(err))
                    wdata = None
                    header = None
                    locdata = None

            if wdata is not None:
                break

    # End wformats for loop and if wdata is None.

    locdata["loc"] = stcode

    if wdata is None:
        print("All attempts to read your file were unsuccesful, " +
              "returning empty table.")
        return wdata, locdata, header

    else:

        if len(np.unique(wdata.year.values)) > 1:
            # Incoming file is probably a TMY or TRY file,
            # so insert a dummy year.
            wdata["year"] = 2017

        wdata_array = (np.vstack([
                wdata.year.values, wdata.month.values,
                wdata.day.values, wdata.hour.values,
                wdata.tdb.values, wdata.tdp.values,
                wdata.rh.values, wdata.ghi.values,
                wdata.dni.values, wdata.dhi.values,
                wdata.wspd.values, wdata.wdr.values
                ])).T
        # Note that the second column is day_of_month here but in the main
        # indra script it will be converted to day_of_year.

        return wdata_array, locdata, header

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
            # The iterator "m" starts at zero.

            if rem <= 0:
                # Iterator has now reached the incomplete month.

                # The previous month is the correct month.
                # Not subtracting 1 because the iterator starts at 0.
                month[d] = m

                # Add the negative remainder to the previous month"s days.
                dom[d] = rem + prev_ndays
                break

            # Subtract number of days in this month from the day.
            rem -= ndays
            # Store number of days from previous month.
            prev_ndays = ndays

    return month, dom

# End function day_of_month


# %%
epw_colnames = ["Year", "Month", "Day", "Hour", "Minute", "QualFlags",
                    "TDB", "TDP", "RH", "ATMPR", "ETRH", "ETRN", "HIR",
                    "GHI", "DNI", "DHI", "GHE", "DNE", "DHE", "ZL",
                    "WDR", "WSPD", "TSKY", "OSKY", "VIS", "CHGT",
                    "PWO", "PWC", "PWT", "AOPT", "SDPT",
                    "SLAST", "UnknownVar1", "UnknownVar2", "UnknownVar3"]


def read_epw(fpath, epw_colnames=epw_colnames):

    # Names of the columns in EPW files. Usually ignore the last
    # three columns.

    # Number of header lines expected.
    hlines = 8

    # Uniform date index for all tmy weather data tables.
    dates = pd.date_range("1/1/2017", periods=8760, freq="H")

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
            if "epw" in f \
                and any(s in f.lower() for s in keywords["tmy"]) \
                    and not any(s in f.lower() for s in keywords["amy"]):
                tmy_filelist.append(f)
    else:
        tmy_filelist = [filepaths]

    didx = 0
    typdata = pd.DataFrame()

    for f in tmy_filelist:

        # Read table, ignoring header lines.
        wdata_typ = pd.read_csv(f, delimiter=",", skiprows=hlines,
                                header=None, names=epw_colnames)

        wdata_typ.index = dates

        if len(wdata_typ.columns) == 35:
            # Some files have three extra columns
            # (usually the TMY files from USDOE).
            # Delete those columns if found.
            wdata_typ = wdata_typ.drop(["unknownvar1", "unknownvar2",
                                        "unknownvar3"], axis=1)

        # Read header and assign all metadata.
        header = list()
        hf = open(f, "r")
        for ln in range(0, hlines):
            header.append(hf.readline())

        infoline = (header[0].strip()).split(",")

        locdata = dict(loc=infoline[1], lat=infoline[6], long=infoline[7],
                       tz=infoline[8], alt=infoline[9], wmo=infoline[5])

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

        print("Could not locate a file with given station name." +
              " Returning empty table.\r\n")

    return typdata, locdata, header

# ----------- END read_epw function -----------


def read_espr(fpath):

    # Missing functionality - reject call if path points to binary file.

    # Uniform date index for all tmy weather data tables.
    dates = pd.date_range("1/1/2017", periods=8760, freq="H")

    fpath_fldr, fpath_name = os.path.split(fpath)
    sitename = fpath_name.split(sep=".")
    sitename = sitename[0]

    with open(fpath, "r") as f:
        content = f.readlines()

    hlines = 12

    # Split the contents into a header and body.
    header = content[0:hlines]

    # Find the year of the current file.
    yline = [line for line in header if "year" in line]

    if "," in yline[0]:
        yline_split = yline[0].split(",")
    else:
        yline_split = yline[0].split()
    year = yline_split[0].strip()

    locline = [line for line in header if ("latitude" in line)][0]
    siteline = [line for line in header if ("site name" in line)][0]

    if "," in locline:
        locline = locline.split(",")
    else:
        locline = locline.split()

    if "," in siteline:
        siteline = siteline.split(",")
    else:
        siteline = siteline.split()

    locdata = dict(loc=siteline[0], lat=locline[1], long=locline[2],
                   tz="00", alt="0000", wmo="000000")
    # ESP-r files do not contain timezone, altitude, or WMO number.

    body = content[hlines:]

    del content

    # Find the lines with day tags.
    daylines = [[idx, line] for [idx, line] in enumerate(body)
                if "day" in line]

    dataout = np.zeros([8760, 11])

    dcount = 0

    for idx, day in daylines:

        # Get the next 24 lines.
        daylist = np.asarray(body[idx+1:idx+25])

        # Split each line of the current daylist into separate strings.
        if "," in daylist[0]:
            splitlist = [element.split(",") for element in daylist]
        else:
            splitlist = [element.split() for element in daylist]

        # Convert each element to a integer, then convert the resulting
        # list to a numpy array.
        daydata = np.asarray([list(map(int, x)) for x in splitlist])

        # Today"s time slice.
        dayslice = range(dcount, dcount+24, 1)

        # This will split the day-month header line on the gaps.
        if "," in day:
            splitday = day.split(",")
        else:
            splitday = day.split(" ")

        # Remove blanks.
        splitday = [x for x in splitday if x != ""]
        splitday = [x for x in splitday if x != " "]

        # Month.
        dataout[dayslice, 0] = np.repeat(int(splitday[-1]), len(dayslice))

        # Day of month.
        dataout[dayslice, 1] = np.repeat(int(splitday[2]), len(dayslice))

        # Hour (of day).
        dataout[dayslice, 2] = np.arange(0, 24, 1)

        # tdb, input is in deci-degrees, convert to degrees.
        dataout[dayslice, 3] = daydata[:, 1]/10

        # tdp is calculated after this loop.

        # rh, in percent.
        dataout[dayslice, 5] = daydata[:, 5]

        # ghi is calculated after this loop.

        # dni, in W/m2.
        dataout[dayslice, 7] = daydata[:, 2]

        # dhi, in W/m2.
        dataout[dayslice, 8] = daydata[:, 0]

        # wspd, input is in deci-m/s.
        dataout[dayslice, 9] = daydata[:, 3]/10

        # wdr, clockwise deg from north.
        dataout[dayslice, 10] = daydata[:, 4]

        dcount += 24

    # tdp, calculated from tdb and rh.
    dataout[:, 4] = petites.tdb2tdp(dataout[:, 3], dataout[:, 5])

    # ghi, in W/m2.
    dataout[:, 6] = dataout[:, 7] + dataout[:, 8]

    # wspd can have bogus values (999)
    dataout[dataout[:, 10] >= 999., 10] = np.nan
    idx = np.arange(0, dataout.shape[0])
    duds = np.logical_or(np.isinf(dataout[:, 10]), np.isnan(dataout[:, 10]))
    int_func = interpolate.interp1d(
            idx[np.logical_not(duds)], dataout[np.logical_not(duds), 10],
            kind="nearest", fill_value="extrapolate")
    dataout[duds, 10] = int_func(idx[duds])

    dataout = np.concatenate((np.reshape(np.repeat(int(year), 8760),
                                         [-1, 1]), dataout), axis=1)

    clmdata = pd.DataFrame(data=dataout, index=dates,
                           columns=["year", "month", "day", "hour",
                                    "tdb", "tdp", "rh",
                                    "ghi", "dni", "dhi", "wspd", "wdr"])

    return clmdata, locdata, header

# ----------- END read_espr function -----------


def give_weather(ts, locdata, stcode, header,
                 masterfile="che_geneva.iwec.a", ftype="espr",
                 s_shift=0, fpath_out=".", std_cols=std_cols):

    if len(ts.shape) == 2:
        ts = np.atleast_3d(ts)

    n_samples = ts.shape[-1]

    success = np.zeros(n_samples, dtype="bool")

    uy = np.asarray(np.unique(ts[:, 0, 0]), dtype=int)

    for n in range(0, n_samples):

        for y, year in enumerate(uy):

            if fpath_out == ".":
                # Make a standardised name for output file.
                filepath = os.path.join(
                        fpath_out, "wf_out_{0}_{1}".format(
                                year, n + s_shift))
            else:
                # Files need to be renamed so strip out the extension.
                fpath_out = fpath_out.replace(".a", "").replace(
                        ".epw", "").replace(".csv", "")

                # Use the name that has come in.
                if n_samples == 1:
                    filepath = fpath_out
                else:
                    filepath = fpath_out + "_{0}_{1}".format(
                                year, n + s_shift)

            ts_curr = ts[ts[:, 0, n] == year, :, n]

            if ftype == "espr":

                if filepath[-2:] != ".a":
                    filepath = filepath + ".a"

                esp_master, locdata, header = read_espr(masterfile)

                # Replace the year in the header.
                yline = [line for line in header if "year" in line]
                yval = yline[0].split(",")
                yline[0] = yline[0].replace(yval[0], str(year))
                header = [yline[0] if "year" in line else line
                          for line in header]
                # Cut out the last new-line character since numpy savetxt
                # puts in a newline character after the header anyway.
                header[-1] = header[-1][:-1]

                # Replace the relevant columns with "new" data.
                esp_master.loc[:, "tdb"] = np.squeeze(
                        ts_curr[:, [c for c, name in enumerate(std_cols)
                                if name == "tdb"]])*10  # deci-Degrees.
                esp_master.loc[:, "rh"] = np.squeeze(
                        ts_curr[:, [c for c, name in enumerate(std_cols)
                                if name == "rh"]])
                esp_master.loc[:, "dhi"] = np.squeeze(
                        ts_curr[:, [c for c,  name in enumerate(std_cols)
                                if name == "dhi"]])
                esp_master.loc[:, "dni"] = np.squeeze(
                        ts_curr[:, [c for c,  name in enumerate(std_cols)
                                if name == "dni"]])
                esp_master.loc[:, "wspd"] = np.squeeze(
                        ts_curr[:, [c for c,  name in enumerate(std_cols)
                                if name == "wspd"]])*10  # deci-m/s.

                # Save month and day to write out to file as separate rows.
                monthday = (esp_master.loc[:, ["day", "month"]]).astype(int)

                # Drop those columns that will not be written out.
                esp_master = esp_master.drop(
                        ["year", "month", "day", "hour", "ghi", "tdp"],
                        axis=1)
                # Re-arrange the columns into the espr clm file order.
                esp_master = esp_master[["dhi", "tdb", "dni", "wspd",
                                         "wdr", "rh"]]
                # Convert all data to int.
                esp_master = esp_master.astype(int)

                master_aslist = esp_master.values.tolist()

                for md in range(0, monthday.shape[0], 25):
                    md_list = [str("* day  {0} month  {1}".format(
                            monthday["day"][md], monthday["month"][md]))]
                    master_aslist.insert(md, md_list)

                # Write the header to file - though the delimiter is
                # mostly meaningless in this case.
                with open(filepath, "w") as f:
                    spamwriter = csv.writer(f, delimiter="\n", quotechar="",
                                            quoting=csv.QUOTE_NONE,
                                            escapechar=" ",
                                            lineterminator="\n")
                    spamwriter.writerow(["".join(header)])

                # Now append the actual data.
                with open(filepath, "a") as f:
                    spamwriter = csv.writer(f, delimiter=",", quotechar="",
                                            quoting=csv.QUOTE_NONE,
                                            lineterminator="\n")
                    for line in master_aslist:
                        spamwriter.writerow(line)

                if os.path.isfile(filepath):
                    success[n] = True
                else:
                    success[n] = False

                # End espr writer.

            elif ftype == "epw":

                if filepath[-4:] != ".epw":
                    filepath = filepath + ".epw"

                epw_fmt = ["%4u", "%2u", "%2u", "%2u", "%2u", "%44s"] + \
                    ((np.repeat("%5.2f", len(epw_colnames)-(6+3))).tolist())

                epw_master, locdata, header = read_epw(masterfile)
                # Cut out the last new-line character since numpy savetxt
                # puts in a newline character after the header anyway.
                header[-1] = header[-1][:-1]

                epw_master.loc[:, "tdb"] = np.squeeze(
                        ts_curr[:, [c for c, name in enumerate(std_cols)
                                if name == "tdb"]])
                epw_master.loc[:, "tdp"] = np.squeeze(
                        ts_curr[:, [c for c, name in enumerate(std_cols)
                                if name == "tdp"]])
                epw_master.loc[:, "rh"] = np.squeeze(
                        ts_curr[:, [c for c, name in enumerate(std_cols)
                                if name == "rh"]])
                epw_master.loc[:, "ghi"] = np.squeeze(
                        ts_curr[:, [c for c, name in enumerate(std_cols)
                                if name == "ghi"]])
                epw_master.loc[:, "dni"] = np.squeeze(
                        ts_curr[:, [c for c, name in enumerate(std_cols)
                                if name == "dni"]])
                epw_master.loc[:, "dhi"] = np.squeeze(
                        ts_curr[:, [c for c,  name in enumerate(std_cols)
                                if name == "dhi"]])
                epw_master.loc[:, "wspd"] = np.squeeze(
                        ts_curr[:, [c for c,  name in enumerate(std_cols)
                                if name == "wspd"]])
                epw_master.loc[:, "wdr"] = np.squeeze(
                        ts_curr[:, [c for c,  name in enumerate(std_cols)
                                if name == "wdr"]])

                np.savetxt(filepath, epw_master.values, fmt=epw_fmt,
                           delimiter=",", header="".join(header),
                           comments="")

                if os.path.isfile(filepath):
                    success[n] = True
                else:
                    success[n] = False

                # End EPW writer.

            else:

                if filepath[-4:] != ".csv":
                    filepath = filepath + ".csv"

                np.savetxt(filepath, np.squeeze(ts_curr), "%5.2f",
                           delimiter=",", header=header + " ".join(std_cols),
                           comments="#")

                if os.path.isfile(filepath):
                    success[n] = True
                else:
                    success[n] = False

    print("You asked for {0} files to be written out. ".format(n_samples) +
          "I was able to write out {0} files successfully.".format(
                  np.sum(success)))
#    return 1

# ----------- End give_weather function. -----------
