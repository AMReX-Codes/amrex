#!/usr/bin/env python

# Copyright 2016-2020 Andrew Myers, David Grote, Maxence Thevenet
# Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


"""
setup.py file for WarpX
"""

import sys
import argparse

from setuptools import setup

argparser = argparse.ArgumentParser(add_help=False)
argparser.add_argument('--with-libwarpx', type=str, default=None, help='Install libwarpx with the given value as DIM. This option is only used by the makefile.')
args, unknown = argparser.parse_known_args()
sys.argv = [sys.argv[0]] + unknown

if args.with_libwarpx:
    package_data = {'pywarpx' : ['libwarpx%s.so'%args.with_libwarpx]}
else:
    package_data = {}

setup (name = 'pywarpx',
       version = '20.03',
       packages = ['pywarpx'],
       package_dir = {'pywarpx':'pywarpx'},
       description = """Wrapper of WarpX""",
       package_data = package_data,
       install_requires=['picmistandard', 'periodictable']
       )
