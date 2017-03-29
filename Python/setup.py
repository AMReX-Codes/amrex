#!/usr/bin/env python

"""
setup.py file for WarpX
"""

from distutils.core import setup

setup (name = 'pywarpx',
       packages = ['pywarpx'],
       package_dir = {'pywarpx':'pywarpx'},
       description = """Wrapper of WarpX""",
       package_data = {'pywarpx' : ['libwarpx.so']},
       )
