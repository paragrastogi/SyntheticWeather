# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 17:57:25 2017

@author: rasto

Check outputs of the weather generator by plotting things.

"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.graphics.tsaplots import plot_pacf

import default_colours as colours
from wfileio import get_weather

from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')

l_start = 0
l_end = 365*24
n_samples = 10
stcode = "lgw"
figpath = stcode

syn_path = os.path.join("lgw", "syn.npy")

try:
    with open(syn_path, "rb") as f:
        ts = pickle.load(f)
except Exception as err:
    ts = np.load(syn_path)

xy_train, locdata, header = get_weather(
        stcode, os.path.join("lgw", "GBR_London_Gatwick.a"), "espr")


column_names = ('year', 'month', 'day', 'hour', 'tdb', 'tdp', 'rh',
                'ghi', 'dni', 'dhi', 'wspd', 'wdr')

# %%

plotrange = range(l_start, l_end)
p = 4
#n = range(0, n_samples+1)

# Line plot.
ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                facecolor='w', edgecolor='k').add_subplot(111)

p2 = plt.plot(plotrange, ts[plotrange, p, 0:9], linewidth=0.5,
              alpha=1)
p1 = plt.plot(plotrange, xy_train[plotrange, p],
              color=colours.blackest, linewidth=1.5)

# ax.set_xticks(range(int(l_start), int(l_end)+l_step, l_step*4))
ax.grid()
plt.xlim(plotrange[0], plotrange[-1])
plt.xlabel('Hour')
plt.title(column_names[p])
plt.ylabel('[degC]')  # Change according to plotted variable.
plt.legend(['Recorded', 'Synthetic'])  # , 'q1', 'q2','q3'])

# figname = os.path.join(figpath, 'line_{0}.pdf'.format('tdb'))
# plt.savefig(figname)
plt.show()

# %%

n = 7

# ACF plot.
ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                facecolor='w', edgecolor='k').add_subplot(111)

p1 = plot_acf(xy_train[plotrange, p], ax=ax, lags=24, alpha=0.05)
p1 = plot_acf(ts[plotrange, p, n], ax=ax, lags=24, alpha=0.05)

plt.xlabel('Lags [Hours]')
plt.ylabel('ACF')
plt.title('TDB')
plt.legend(['', 'Recorded', '', 'Synthetic'])

# PACF plot.
ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                facecolor='w', edgecolor='k').add_subplot(111)

p1 = plot_pacf(xy_train[plotrange, p], ax=ax, lags=72, alpha=0.05)
p1 = plot_pacf(ts[plotrange, p, n], ax=ax, lags=72, alpha=0.05)

plt.xlabel('Lags [Hours]')
plt.ylabel('PACF')
plt.title('TDB')
plt.legend(['', 'Recorded', '', 'Synthetic'])

# Corr plot.
ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                facecolor='w', edgecolor='k').add_subplot(111)

p1 = plt.scatter(xy_train[plotrange, p], xy_train[plotrange, p+1],
                 color=colours.blackest)
p1 = plt.scatter(ts[plotrange, p, n], ts[plotrange, p+1, n],
                 color=colours.reddest)

plt.xlabel('Lags [Hours]')
plt.ylabel('Corr.')
plt.title('TDB vs TDP')
plt.legend(['Recorded', 'Synthetic'])
