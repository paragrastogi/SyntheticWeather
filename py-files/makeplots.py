#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 16:44:06 2017

@author: rasto
"""



# =============================================================================
    #
    # # This bit is only for plotting (to see how generated values correspond
    # # to summary statistics) and has nothing to do with the generator.
    #
    # # Calculate the iqr and median to plot later.
    # actualq3 = gw.weather_stats(historicaldata, ['month'], 'q3')
    # actualmed = gw.weather_stats(historicaldata, ['month'], 'med')
    # actualq1 = gw.weather_stats(historicaldata, ['month'], 'q1')
    #
    # temp = np.repeat(actualq1.values,
    #                  int(currentdata.shape[0]/actualq1.shape[0]), axis=0)
    # temp2 = pd.DataFrame(data=temp, columns=actualq1.columns)
    # actualq1 = temp2
    # del (temp, temp2)
    #
    # temp = np.repeat(actualq3.values,
    #                  int(currentdata.shape[0]/actualq3.shape[0]), axis=0)
    # temp2 = pd.DataFrame(data=temp, columns=actualq3.columns)
    # actualq3 = temp2
    # del (temp, temp2)
    #
    # temp = np.repeat(actualmed.values,
    #                  int(currentdata.shape[0]/actualmed.shape[0]), axis=0)
    # temp2 = pd.DataFrame(data=temp, columns=actualmed.columns)
    # actualmed = temp2
    # del (temp, temp2)

    # ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
    #                facecolor='w', edgecolor='k').add_subplot(111)
    # plt.plot(actualq1.tdb.values)
    # plt.plot(actualq3.tdb.values)
    # plt.show

    # =============================================================================

    # %%

#    plotrange = range(l_start, l_end)
#    p = 3
#    n = range(0, n_samples+1)
#
#    # Line plot.
#    ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
#                    facecolor='w', edgecolor='k').add_subplot(111)
#
#    p1 = plt.plot(plotrange, ts_curr_in[plotrange, p],
#                  color=colours.orange, linewidth=1.5, zorder=1)
#    p2 = plt.plot(plotrange, xout_un[plotrange, p, :], linewidth=0.5,
#                  alpha=1)
#    #p31 = plt.plot(plotrange, actualq1.tdb.values, linewidth=2,
#    #               alpha=1, color=colours.blackest, zorder = 3, marker='o')
#    #p32 = plt.plot(plotrange, actualmed.tdb.values, linewidth=2,
#    #               alpha=1, color=colours.blackest, zorder = 4, marker='o')
#    #p33 = plt.plot(plotrange, actualq3.tdb.values, linewidth=2,
#    #               alpha=1, color=colours.blackest, zorder = 5, marker='o')
#
#    # ax.set_xticks(range(int(l_start), int(l_end)+l_step, l_step*4)) color=colours.orange,
#    ax.grid()
#    plt.xlim(plotrange[0], plotrange[-1])
#    plt.xlabel('Hour')
#    plt.title(column_names[p])
#    plt.ylabel('[degC]') # Change according to plotted variable.
#    plt.legend(['Recorded', 'Synthetic']) # , 'q1', 'q2','q3'])
#
#    figname = os.path.join(figpath, 'line_{0}.pdf'.format('tdb'))
#    #plt.savefig(figname)
#    plt.show()
#
#    # %%
#
#    p = 2
#    n = 0  # random.randint(0, n_samples)
#    plotrange = range(l_start, l_end)
#
#    # ACF plot.
#    ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
#                    facecolor='w', edgecolor='k').add_subplot(111)
#
#    p1 = plot_acf(ts_curr_in[plotrange, p], ax=ax, lags=72, alpha=0.05)
#    p1 = plot_acf(xout_un[plotrange, p, n], ax=ax, lags=72, alpha=0.05)
#
#    plt.xlabel('Lags [Hours]')
#    plt.ylabel('ACF')
#    plt.title('TDB')
#    plt.legend(['', 'Recorded', '', 'Synthetic'])
#
#    # PACF plot.
#    ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
#                    facecolor='w', edgecolor='k').add_subplot(111)
#
#    p1 = plot_pacf(ts_curr_in[plotrange, p], ax=ax, lags=72, alpha=0.05)
#    p1 = plot_pacf(xout_un[plotrange, p, n], ax=ax, lags=72, alpha=0.05)
#
#    plt.xlabel('Lags [Hours]')
#    plt.ylabel('PACF')
#    plt.title('TDB')
#    plt.legend(['', 'Recorded', '', 'Synthetic'])
#
#    ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
#                    facecolor='w', edgecolor='k').add_subplot(111)
#
#    p1 = plt.scatter(ts_curr_in[plotrange, 2], ts_curr_in[plotrange, 3], color = colours.blackest)
#    p1 = plt.scatter(xout_un[plotrange, 2, n], xout_un[plotrange, 3, n],  color=colours.reddest)
#
#    plt.xlabel('Lags [Hours]')
#    plt.ylabel('Corr.')
#    plt.title('TDB vs TDP')
#    plt.legend(['Recorded', 'Synthetic'])
#    # %%
#
#    # Repeat the RGB spec for grey n_sample times to plot n_sample synthetic series.
#    grey_reps = np.repeat(np.atleast_2d(colours.grey).T,n_samples, axis=1).T
#
#    # Histogram.
#    ax = plt.figure(num=None, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k').add_subplot(111)
#
#
#    p1 = plt.hist(ts_curr_in[plotrange, p], normed=True, cumulative=False,
#                  bins=50, histtype='step', color=colours.blackest, linewidth=2,
#                  zorder=1)
#    p2 = plt.hist(xout_un[plotrange, p, :], normed=True, cumulative=False,
#                  bins=50, histtype='step', color=grey_reps)
#    # =============================================================================
#    # p3 = plt.hist(historicaldata.tdb.values[np.logical_not(np.isnan(
#    #         historicaldata.tdb.values))], normed=True, cumulative=False,
#    #               bins=50, histtype='step', zorder=2, linewidth=2, color=colours.orange)
#    # =============================================================================
#
#    plt.xlabel('Temp.')
#    plt.ylabel('CDF')
#    plt.legend(['Recorded', 'Synthetic'], loc='upper right')
#
#    figname = os.path.join(figpath, 'hist_{0}.pdf'.format('tdb'))
#    plt.savefig(figname)
#    plt.show()
#    # %% TDP
#
#    plotrange = range(l_start, l_end)
#    p = 5
#
#    ax = plt.figure(num=None, figsize=(8, 6), dpi=80,
#                    facecolor='w', edgecolor='k').add_subplot(111)
#
#    p1 = plt.plot(plotrange, ts_curr_in[plotrange, p],
#                  color=colours.orange, linewidth=0.5, zorder=n_samples+1)
#    p2 = plt.plot(plotrange, xout_un[plotrange, p, :], linewidth=0.75,
#                  alpha=1, color=colours.grey, zorder=1)
#    #p2 = plt.plot(plotrange, actualstats.tdb, linewidth=0.75,
#    #              alpha=1, color=colours.lgrey)
#
#    # ax.set_xticks(range(int(l_start), int(l_end)+l_step, l_step*4))
#    # labels = [item.get_text() for item in ax.get_xticklabels()]
#    # labels = range(int(l_start/24), int(l_end/24)+1, 4)
#    # ax.set_xticklabels(labels)
#    ax.grid()
#    # plt.xlim(l_start, l_end)
#    plt.xlabel('Hour')
#    plt.ylabel('[degC]')
#    plt.title(column_names[p])
#    plt.legend(['Recorded', 'Synthetic'])
#
#    figname = os.path.join(figpath, 'line_{0}.pdf'.format('tdp'))
#    plt.savefig(figname)
#    plt.show()
#
#    # %%
#
#    # Histogram.
#    plt.figure(num=None, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k')
#
#    p1 = plt.hist(ts_curr_in[plotrange, p], normed=True, cumulative=False,
#                  bins=50, histtype='step', color=colours.orange, linewidth=2,
#                  zorder=1)
#    p2 = plt.hist(xout_un[plotrange, p, :], normed=True, cumulative=False,
#                  bins=50, histtype='step', zorder=n_samples+1)
#
#    plt.xlabel('Temp.')
#    plt.ylabel('CDF')
#    plt.legend(['Recorded', 'Synthetic'], loc='upper right')
#
#    figname = os.path.join(figpath, 'hist_{0}.pdf'.format('tdp'))
#    plt.savefig(figname)
#    plt.show()
