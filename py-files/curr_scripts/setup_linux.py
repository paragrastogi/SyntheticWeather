from cx_Freeze import setup, Executable

import os
import scipy

includefiles_list = ['CityData.csv']

scipy_path = os.path.dirname(scipy.__file__)
includefiles_list.append(scipy_path)

# Dependencies are automatically detected, but it might need
# fine tuning.
buildOptions = dict(
        packages=['numpy.lib.format', 'numpy.core._methods', 'numpy.matlib'],
        excludes=[],
        include_files=includefiles_list)

base = 'Console'

executables = [
    Executable('listener.py', base=base, targetName='indra')
]

setup(name='indra',
      version='1.0',
      description='synweather',
      options= dict(build_exe = buildOptions),
      executables=executables)
