#!/usr/bin/env python

"""
setup.py file for WarpX
"""

from distutils.core import setup

try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    from distutils.command.build_py import build_py


setup (name = 'pywarpx',
       packages = ['pywarpx'],
       package_dir = {'pywarpx':'pywarpx'},
       description = """Wrapper of WarpX""",
       package_data = {'pywarpx' : ['libwarpx.so']},
       cmdclass={'build_py': build_py}
       )
