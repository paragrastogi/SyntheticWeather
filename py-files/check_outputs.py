# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 17:57:25 2017

@author: rasto

Check outputs of the weather generator by plotting things.

"""

import matplotlib.pyplot as plt

l_start = 0
l_end = 31*24
n_samples = 10


plotrange = range(l_start, l_end)
p = 3
n = range(0, n_samples+1)

# Line plot.
ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
                facecolor='w', edgecolor='k').add_subplot(111)

p1 = plt.plot(plotrange, ts_curr_in[plotrange, p],
              color=colours.orange, linewidth=1.5, zorder=1)
p2 = plt.plot(plotrange, xout_un[plotrange, p, :], linewidth=0.5,
              alpha=1)
#p31 = plt.plot(plotrange, actualq1.tdb.values, linewidth=2,
#               alpha=1, color=colours.blackest, zorder = 3, marker='o')
#p32 = plt.plot(plotrange, actualmed.tdb.values, linewidth=2,
#               alpha=1, color=colours.blackest, zorder = 4, marker='o')
#p33 = plt.plot(plotrange, actualq3.tdb.values, linewidth=2,
#               alpha=1, color=colours.blackest, zorder = 5, marker='o')

# ax.set_xticks(range(int(l_start), int(l_end)+l_step, l_step*4)) color=colours.orange,
ax.grid()
plt.xlim(plotrange[0], plotrange[-1])
plt.xlabel('Hour')
plt.title(column_names[p])
plt.ylabel('[degC]') # Change according to plotted variable.
plt.legend(['Recorded', 'Synthetic']) # , 'q1', 'q2','q3'])

figname = os.path.join(figpath, 'line_{0}.pdf'.format('tdb'))
#plt.savefig(figname)
plt.show()
 