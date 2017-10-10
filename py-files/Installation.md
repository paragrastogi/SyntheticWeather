% Install and Run indra
% Parag Rastogi
% 29 September, 2017

[This readme was originally written in markdown.]

Hello!

Thank you for downloading **indra**, a synthetic weather generator created by Parag Rastogi.

**indra** is a collection of python scripts. It is written for python 3.5+ and not backwards with older versions of python for now (October 2017)<sup>(1)</sup>. To begin, you need the following modules installed:

1. `numpy`
2. `scipy`
3. `sklearn`
4. `pandas`
5. `statsmodels`

The following modules are also used but should be part of your standard python distro:

1. `os`
2. `csv`
3. `sys`
4. `random`
5. `pickle`
6. `argparse`

If these are not installed for some reason, please do install them.

Sometimes, there are issues with the mkl-optimised versions of numpy, scipy, and sklearn distributed by Anaconda. In that case, install the package `nomkl` and these modules will be 'downgraded' to their `*_nomkl` versions. That shouldn't affect the functioning of these packages at all.

In general, if you are in the directory where all the scripts for **indra** are, i.e., the directory where this file is, then you do not need to write the full path to any of the files in the commands described below. __NB:__ When using this tool from the python prompt or the Windows/Linux command line, you use `call_indra.py`. If you compile it to a command-line program, then you use `indra.exe` or `./indra`. If you look at the script `call_indra.py`, you will understand why this is the case.

--------------------

## Run with python

If you just want to work in python, then get started with the script listener.py by typing the following into your command line:

### Linux/Ubuntu

    $ python <PATH_TO_LISTENER>\call_indra.py -h

### Windows

    > python <PATH_TO_LISTENER>/call_indra.py -h

--------------------

## Compile and Run as a command-line executable.

If you want to compile **indra** into an executable, use something like cx_freeze. An executable compiled on Windows will work on windows only, and vice-versa. I tried nuitka but couldn't quite get it to work on Linux.

`cx_freeze` is a set of scripts you can use to compile **indra** into an executable. You don't need to do this if you want to work in python itself. The module `cx_freeze` has to be installed using pip, or download it manually from: https://sourceforge.net/projects/cx-freeze/

To install using the `setup_linux.py` or `setup_windows.py` script in this folder, run the following command:

### Linux/Ubuntu

    $ python setup_linux.py build

### Windows

    > python setup_windows.py build

Upon compilation you should end up with a bunch of files inside a folder call `build/<something>`. The name of the subfolder `<something>` depends on the platform you are on. In this folder, there will be an executable called **indra** (on Linux) or **indra**.exe (on Windows). Don't delete any of these apparently 'extra' files - these are all the python libraries compiled in C code!

Congratulations, you now have an executable (program) called `indra` on your hands! Chang directory (`cd`) into `build\<something>` to start running **indra**. Or just rearrange things so you don't have horribly long paths. Just do it, you're smart.

### Linux/Ubuntu

You could run the command

    chmod +x indra

and make the file executable or just use the `.\indra` syntax. If you add the folder containing the executable and its associate files to your `PATH` variable, you don't need to type the full path at the command line.

    $ <PATH_TO_INDRA>\indra -h

### Windows

Like Linux, if you add the folder containing the executable and its associate files to your `PATH` variable, you don't need to type the full path at the command line.

    > <PATH_TO_INDRA>/indra.exe -h

--------------------

### Footnotes

(1) For now, I cannot guarantee compatibility with Python 2.x. On most Linux machines, Python 2.x is installed by default and the command python points to that version. To install and run this program, therefore, please use python3 and pip3. For example:

    $ python3 call_indra.py -h

    $ pip3 install pandas
