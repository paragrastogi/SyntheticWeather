# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 20:29:16 2017

@author: rasto
"""
import numpy as np
from statsmodels.tsa.statespace.sarimax import SARIMAX
from sys import stdout
# from tqdm import tqdm


def select_models(arp, maq, sarp, smaq, s, ts_in):

    aic_curr = 0
    selaic = np.infty

    counter = 0

    # Loop through all possible combinations of ar, ma, sar, and sma lags.

    for p in arp:
        for q in maq:
            for pp in sarp:
                for qq in smaq:

                    if p == 0 and q == 0:
                        continue

                    model = SARIMAX(
                        ts_in, order=(p, 0, q),
                        seasonal_order=(pp, 0, qq, s),
                        trend=None)
                    # model_type.append('s')

                    try:
                        mod_fit_curr = model.fit(
                                disp=0, cov_type="robust",
                                full_output=True)
                        aic_curr = mod_fit_curr.aic

                        if np.isnan(aic_curr):
                            continue

                        else:

                            if counter > 0:
                                if aic_curr < selaic:
                                    selaic = aic_curr
                                    selmdl = mod_fit_curr
                                # else:
                                #    print("aic_curr > selaic")

                            elif counter == 0:
                                selaic = aic_curr
                                selmdl = mod_fit_curr

                            counter += 1
                            # print("counter {0}".format(counter))
                            # print("aic_curr {0}".format(aic_curr))
                            # print("selaic {0}".format(selaic))

                    except Exception as err:
                        # print('fit threw an error')
                        continue

                    # Print out a heartbeat.
                    stdout.write("...{0}".format(counter))
                    stdout.flush()

                # End qq loop.
            # End pp loop.
        # End q loop.
    # End p loop.

    resid = selmdl.resid

    return selmdl, resid
