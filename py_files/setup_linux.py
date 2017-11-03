"""
@author: Parag Rastogi

Run this script to compile indra into an executable for linux.
The 'outermost' script is listener.py, which parses inputs and calls the
main script 'indra'. Indra, in turn, calls a bunch of scripts for various
tasks.

This script was created from the cxfreeze-quickstart script and modified
by Parag.
"""

# cxfreeze is a python module used to compile python script(s) into an
# executable on linux or windows.
from cx_Freeze import setup, Executable

# os is useful module for various operating system-related functions.
import os
# scipy isn't actually used in this script, but it is used in indra.
# There are some issues comipiling all of scipy, so it needs some files
# to be explicitly included.
import scipy

# This csv is just a list of stations with their metadata. To be used
# in a future release to identify files/sites by latitude, longitude, etc.
includefiles_list = ['CityData.csv']

# Explicitly add scipy file.
scipy_path = os.path.dirname(scipy.__file__)
includefiles_list.append(scipy_path)

# Dependencies are automatically detected, but it might need fine tuning.
# In this case, some numpy libraries have to be mentioned explicitly.
buildOptions = dict(
        packages=['numpy.lib.format', 'numpy.core._methods', 'numpy.matlib'],
        excludes=[],
        include_files=includefiles_list)

# Script will work from console.
base = 'Console'

# Compile file listener.py to work in console with name indra.
executables = [
    Executable('call_indra.py', base=base, targetName='indra')
]

# Call cxfreeze.
setup(name='indra',
      version='1.0',
      description='synweather',
      options=dict(build_exe=buildOptions),
      executables=executables)
