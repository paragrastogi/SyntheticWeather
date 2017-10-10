# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 20:29:16 2017

@author: rasto
"""
import numpy as np
from statsmodels.tsa.statespace.sarimax import SARIMAX
# from statsmodels.tsa.arima_model import ARIMA


def select_models(arp, maq, sarp, smaq, s, ts_in):

    model_fit = list()
    model_type = list()
    aic = np.zeros(len(arp) * len(maq) * len(sarp) * len(smaq))

    pp = 0
    qq = 0

    counter = 0
    # Model for TDB.
    for p in arp:
        for q in maq:
            for pp in sarp:
                for qq in smaq:

                    if p == 0 and q == 0:
                        aic[counter] = None
                        counter += 1
                        model_fit.append(None)
                        model_type.append(None)
                        continue

#                    if pp == 0 and qq == 0:
#                        print('fitting arima.')
#                        model = ARIMA(ts_in, order=(p, 0, q), trend=None)
#                        model_type.append('a')
#                    else:
                    print('fitting sarima.')
                    model = SARIMAX(
                            ts_in, order=(p, 0, q),
                            seasonal_order=(pp, 0, qq, s),
                            trend=None)
                    model_type.append('s')

                    try:
                        mod_temp = model.fit(disp=0)
                        aic[counter] = mod_temp.aic

                    except Exception as err:
                        print('fit threw an error')
                        mod_temp = None
                        aic[counter] = None

                    model_fit.append(mod_temp)

                    counter += 1

    aicfilter = aic == np.nanmin(aic)
    selmdl = [model for (model, idx) in zip(model_fit, aicfilter)
              if idx]
    selmdl_type = [mod_type for (mod_type, idx) in
                   zip(model_type, aicfilter) if idx]

    if isinstance(selmdl, (list, tuple)):
        selmdl = selmdl[0]
        selmdl_type = selmdl_type[0]

    resid = selmdl.resid

    return selmdl, selmdl_type, resid
