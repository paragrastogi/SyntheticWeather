"""
This file can be used to install all of the libraries that would be
needed to run the CreateSyntheticWeather.py script. This doesn't need
to be run everytime - just once. This script assumes you do not have admin
rights on the machine. If you do have admin rights, and wish to install the
python packages for everybody, then remove the options
'--user' from the pip.install command.
"""

import subprocess

__author__ = 'Parag Rastogi'

subprocess.run('pip.install --user numpy scipy pandas matplotlib', shell=True)
# Remove '--user' for an admin install

# Write an R file to install the forecast package.
# The repo used is a global repo redirect.
# If the link doesn't work, then this script will fail.
# In that case, change the repo and try again.

f = open('install-packages.r', 'w', encoding='utf8')

# The original of the code written below is at
# https://orfe.princeton.edu/help/r-packages
f.write('## Create the personal library if it doesn''t exist. ' +
        'Ignore a warning if the directory already exists.\n' +
        'dir.create(Sys.getenv("R_LIBS_USER"), ' +
        'showWarnings = FALSE, recursive = TRUE)\n' +
        '## Install one package.\n' +
        'install.packages("forecast", Sys.getenv("R_LIBS_USER"), ' +
        'repos ="http://cloud.r-project.org")\n' +
        # The rest are other types of installations.
        # Only included here for reference.
        '## Install a package that you have copied to the remote system.\n' +
        '## install.packages("file_name.tar.gz", Sys.getenv("R_LIBS_USER")\n' +
        '## Install multiple packages.\n' +
        '## install.packages(c("timeDate","robustbase"),' +
        ' Sys.getenv("R_LIBS_USER"), ' +
        'repos = "http://cloud.r-project.org")")')

f.close()

# # Compose the R command
rcmd = 'R CMD BATCH install-packages.r '
subprocess.run(rcmd)
