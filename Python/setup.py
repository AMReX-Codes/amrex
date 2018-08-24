#!/usr/bin/env python

"""
setup.py file for WarpX
"""

import sys
import argparse

from distutils.core import setup

try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    from distutils.command.build_py import build_py

argparser = argparse.ArgumentParser(add_help=False)
argparser.add_argument('--with-libwarpx', action='store_true', help='Install libwarpx. This option is only used by the makefile.')
args, unknown = argparser.parse_known_args()
sys.argv = [sys.argv[0]] + unknown

if args.with_libwarpx:
    package_data = {'pywarpx' : ['libwarpx.so']}
else:
    package_data = {}

setup (name = 'pywarpx',
       packages = ['pywarpx'],
       package_dir = {'pywarpx':'pywarpx'},
       description = """Wrapper of WarpX""",
       package_data = package_data,
       cmdclass={'build_py': build_py}
       )
