#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# CAN'T GET THIS TO FUNCTION JUST YET - FOR NOW THE USER CAN INSTALL ALL
# PACKAGES MANUALLY.

"""
Created on Fri Sep 29 16:37:26 2017

@author: parag rastogi
"""
import argparse


def depends(installer):

    try:
        import sys
    except Exception as err:
        print("I need you to install module 'sys' before anything else." +
              "Do that and run this script again.")
        return 0

    try:
        import subprocess
    except Exception as err:
        print("I need you to install module 'subprocess' as well." +
              "Do that and run this script again.")
        return 0

    print("Trying to install everything with installer {0}...\r\n".format(
            installer))

    if 'nomkl' in sys.modules:
        print("nomkl already installed, hooray! \r\n")
    else:
        try:
            p = subprocess.Popen('{0} install nomkl'.format(installer),
                                 shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
            for line in p.stdout.readlines():
                print(line)
            retval = p.wait()
        except Exception as err:
            print("I could not install numpy. Try installing it manually.")

    return 0

    if 'numpy' in sys.modules:
        print("numpy already installed, hooray! \r\n")
    else:
        try:
            p = subprocess.Popen('{0} install numpy'.format(installer),
                                 shell=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
            for line in p.stdout.readlines():
                print(line)
            retval = p.wait()
        except Exception as err:
            print("I could not install numpy. Try installing it manually.")

    print(retval)

    if 'scipy' not in sys.modules:
        try:
            subprocess.Popen('{0} install scipy'.format(installer),
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        except Exception as err:
            print("I could not install scipy. Try installing it manually.")
    else:
        print('scipy already installed, hooray! \r\n')

    if 'sklearn' not in sys.modules:
        try:
            subprocess.Popen('{0} install sklearn'.format(installer),
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        except Exception as err:
            print("I could not install sklearn. Try installing it manually.")
    else:
        print('sklearn already installed, hooray! \r\n')

    if 'cx_freeze' not in sys.modules:
        try:
            print("cx_freeze is a set of scripts you can use to compile " +
                  "indra into an executable. You don't need to do this " +
                  "if you want to work in python itself. cx_freeze has to " +
                  "be installed using pip, or download it manually from:" +
                  "https://sourceforge.net/projects/cx-freeze/ \r\n")
            subprocess.Popen('pip install cx_freeze'.format(installer),
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        except Exception as err:
            print("I could not install cxfreeze. Try installing it manually.")

# END depends function.


# This bit makes it callable from the command line.
parser = argparse.ArgumentParser(
        description="This script checks dependencies for Indra, a " +
        "generator of synthetic weather time series. This script depends " +
        "on three modules already being installed: sys, argparse, and " +
        "subprocess. These should be part of your standard python " +
        "distro. Indra is written for python 3.5+ and not backwards " +
        "with older versions of python at all.\r\n")

parser.add_argument('--installer', type=str, choices=['pip', 'conda'],
                    default='pip', help='Enter the name of your python ' +
                    'package installer. Options are [pip] and [conda].')
args = parser.parse_args()

installer = args.installer

if __name__ == '__main__':
    depends(installer)
