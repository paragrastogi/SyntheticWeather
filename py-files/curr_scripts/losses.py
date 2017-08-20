import numpy as np
"""
Calculate loss functions for incoming vectors of error.
Add any error functions below - just remember to explicitly
import them in the script you are calling them from.
"""


def rmseloss(e):
    return np.sqrt((e**2).mean())


def maeloss(e):
    return np.abs(e.mean())
