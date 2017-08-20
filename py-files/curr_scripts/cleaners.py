import numpy as np
import pandas as pd


def solar_col_clean(wdata2):

    copy = wdata2.copy()

    metstat = copy.GHI_metstat
    measure = copy.GHI_measure

    idx = 0
    temp = copy.GHI_suny

    copy = copy.drop(['GHI_suny', 'GHI_metstat', 'GHI_measure'], axis=1)

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

    copy = pd.concat([temp, copy], axis=1, join='inner')
    copy = copy.rename(columns={'GHI_suny': 'ghi'})

    # Do the same for the DNI columns.

    metstat = copy.DNI_metstat
    measure = copy.DNI_measure

    idx = 0
    temp = copy.DNI_suny

    copy = copy.drop(['DNI_suny', 'DNI_metstat', 'DNI_measure'], axis=1)

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

    copy = pd.concat([temp, copy], axis=1, join='inner')
    copy = copy.rename(columns={'DNI_suny': 'dni'})

    # And for DHI.

    metstat = copy.DHI_metstat
    measure = copy.DHI_measure

    idx = 0
    temp = copy.DHI_suny

    copy = copy.drop(['DHI_suny', 'DHI_metstat', 'DHI_measure'], axis=1)

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

    copy = pd.concat([temp, copy], axis=1, join='inner')
    copy = copy.rename(columns={'DHI_suny': 'dhi'})

    return copy
