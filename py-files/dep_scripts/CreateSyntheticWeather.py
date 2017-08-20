"""This is a port of the script originally written with my thesis,
 called CreateSyntheticWeather.m .
 Script originally written by Parag Rastogi
"""

import numpy as np
# import matplotlib.pyplot as plt
import os.path
import pandas as pd
import scipy as sp
# from scipy import optimize
# from scipy import io  # This is to read MAT files.
# import code
import fourier_fit_funcs as fourier
import SimArima
# import csv

# This is just to time the execution of the script.
# Doesn't serve any programmatic purpose.
from datetime import datetime
startTime = datetime.now()

__author__ = 'Parag Rastogi'

print('Welcome to the Synthetic Weather Time Series Generator, v0\r\n')
# location = input('Please enter the name of the location you will be '
#                  'working with today. The first three letters of the '
#                  'string entered will be used as the location/region label.'
#                  '\r\n')
location = 'GENEVA'
loccode = location.upper()[0:3]

# This is the letter of the drive where the files and data are/will be
# drive = 'e:\\'
pathInEPWfldr = os.path.join('..', 'GEN')
nameInEPWfile = 'GEN_IWEC.epw'
pathInEPWfile = os.path.join(pathInEPWfldr, nameInEPWfile)

# These are the names of the columns in EPW files
colnames = ['Year', 'Month', 'Day', 'Hour', 'Minute', 'QualFlags',
            'TDB', 'TDP', 'RH', 'ATMPR', 'ETRH', 'ETRN', 'HIR', 'GHI',
            'DNI', 'DHI', 'GHE', 'DNE', 'DHE', 'ZL', 'WDR', 'WSPD', 'TSKY',
            'OSKY', 'VIS', 'CHGT', 'PWO', 'PWC', 'PWT', 'AOPT', 'SDPT',
            'SLAST', 'UnknownVar1', 'UnknownVar2', 'UnknownVar3']
# Convert the names to be in lowercase
colnames = [x.lower() for x in colnames]
tmytable = pd.read_csv(pathInEPWfile, delimiter=',', skiprows=8,
                       header=None, names=colnames)

# Initialise time index variables
N = tmytable.shape[0]  # Should be 8760
t = np.linspace(1, N, endpoint=True, num=N)  # Julian Hour Index
nm = 12  # Number of months
# mt = (2015,1,1) # Month number

# These are the 'raw' time series, i.e., read from the EPW file.
raws = dict(tdb=np.array(tmytable.tdb) + 273.15, rh=np.array(tmytable.rh),
            ghi=np.array(tmytable.ghi),
            dni=np.array(tmytable.dhi), dhi=np.array(tmytable.dhi))
# RH -> Relative Humidity (%)
# GHI -> Global Horizontal Solar Radiation [Wh/m2]
# DNI -> Direct Normal Solar Radiation [Wh/m2]
# DHI -> Diffuse Horizontal Solar Radiation [Wh/m2]

# Limit below which a reading is considered to be ZERO
ghi_limit = 1  # w/m2
# This limit applies to all solar radiation quantities.
# That is to say that if the global horizontal radiation
# was effectively zero, everything else will be zero too.

ghicensor = raws['ghi'] <= ghi_limit
raws['ghi'][np.ix_(ghicensor)] = 0
raws['ghi'][np.ix_(ghicensor)] = 0
raws['ghi'][np.ix_(ghicensor)] = 0

# The curve_fit function outputs two things:
params = {'tdb': sp.optimize.curve_fit(fourier.fit_tdb, t, raws['tdb']),
          'rh': sp.optimize.curve_fit(fourier.fit_rh, t, raws['rh']),
          'tdb_low': sp.optimize.curve_fit(fourier.fit_tdb_low, t,
                                           raws['tdb']),
          'tdb_high': sp.optimize.curve_fit(fourier.fit_tdb_high, t,
                                            raws['tdb']),
          'rh_low': sp.optimize.curve_fit(fourier.fit_rh_low, t, raws['rh'])}

# Call the fourier fit function with the calculated parameters to get the
# values of the fourier fit at each time step
ffit = {'tdb': fourier.fit('tdb', t, *params['tdb'][0]),
        'rh': fourier.fit('rh', t, *params['rh'][0]),
        'tdb_low': fourier.fit('tdb_low', t, *params['tdb_low'][0]),
        'tdb_high': fourier.fit('tdb_high', t, *params['tdb_high'][0]),
        'rh_low': fourier.fit('rh_low', t, *params['rh_low'][0])}

# Now subtract the low- and high-frequency fourier fits
# (whichever is applicable) from the raw values to get the
# 'de-meaned' values (values from which the mean has been removed).
demeaned = dict(tdb=raws['tdb'] - ffit['tdb'],
                rh=raws['rh'] - ffit['rh'])

# Declare the parameters that will be fed to the R script.
# If arp or maq>1, then the R script will try to fit models
# with each combination of arp and maq.
# The script will then output the best model (based on smallest BIC).
arp = 4
maq = 4
sarp = 1
smaq = 1
npaths = 50

# Write the de-meaned data series to CSV file.
pathTSdm = os.path.join(pathInEPWfldr,
                        '{0}_dm.csv'.format(nameInEPWfile[0:-4]))

# Construct a 2D array to write to csv file
dm = np.zeros((demeaned['tdb'].shape[0], 2))
dm[:, 0] = demeaned['tdb']
dm[:, 1] = demeaned['rh']

np.savetxt(pathTSdm, dm, delimiter=',', fmt='%f')

print('Calling the R script... \r\n')

# Call R inside the following script
[bsout, modlist] = SimArima.callrcmd(pathTSdm, arp, maq, sarp, smaq, npaths)

# Add the results from simulating the SARMA model (i.e.,
# the output from R) back to the seasonal/fourier fits that
# were removed from the original data earlier.
tdbrecon = np.zeros(bsout['tdb'].shape)

# code.interact(local={**locals(), **globals()})

for path in range(0, bsout['tdb'].shape[1], 1):
    tdbrecon[:, path] = ffit['tdb'] + bsout['tdb'][:, path]

rhrecon = np.zeros(bsout['rh'].shape)
for path in range(0, bsout['rh'].shape[1], 1):
    rhrecon[:, path] = ffit['rh'] + bsout['rh'][:, path]

arma_deg = dict(tdb=modlist['tdb'][0].split(' '),
                rh=modlist['rh'][0].split(' '))
for e in range(0, len(arma_deg['tdb'])):
    arma_deg['tdb'][e] = int(arma_deg['tdb'][e])
for e in range(0, len(arma_deg['rh'])):
    arma_deg['rh'][e] = int(arma_deg['rh'][e])

arma_coeff = dict(tdb=modlist['tdb'][1].split(' '),
                  rh=modlist['rh'][1].split(' '))
for e in range(0, len(arma_coeff['tdb'])):
    arma_coeff['tdb'][e] = float(arma_coeff['tdb'][e])
for e in range(0, len(arma_coeff['rh'])):
    arma_coeff['rh'][e] = float(arma_coeff['rh'][e])

# Create a dictionary to save the ARMA model
ARmodel = dict(tdb=dict(degrees=arma_deg['tdb'], coeffs=arma_coeff['tdb'],
                        bic=float(modlist['tdb'][2])),
               rh=dict(degrees=arma_deg['rh'], coeffs=arma_coeff['rh'],
                       bic=float(modlist['rh'][2])))

# This AR model should be written to a file so it can be shown to the user.

# recon = {'tdb': tdbrecon, 'rh': rhrecon}

# Write the synthetic outputs to files
np.savetxt('tdb_syn.csv', tdbrecon)
np.savetxt('rh_syn.csv', rhrecon)

print('End of current script ....\r\n')
print('Script took {0:0.1f} seconds'.format(datetime.now() - startTime))

# Disabled to run on server.
# Go to the command line with code.interact. The stuff inside the brackets is
# to make sure that the variables are transmitted.


# except (RuntimeError,TypeError,ValueError):
#     # Go to the command line with code.interact.
#     # The stuff inside the brackets is
#     # to make sure that the variables are transmitted.
# code.interact(local={**locals(), **globals()})
