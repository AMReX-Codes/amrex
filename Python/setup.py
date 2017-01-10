#!/usr/bin/env python

"""
setup.py file for WarpX
"""

from distutils.core import setup, Extension

import numpy

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

include_dirs = [numpy_include, '../Source']

example_module = Extension('pywarpx._warpxC',
                           swig_opts=['-outdir','pywarpx'],
                           sources=['warpxC.i'],
                           library_dirs=['.'],
                           libraries=['warpx'],
                           include_dirs = include_dirs,
                           )

setup (name = 'pywarpx',
       packages = ['pywarpx'],
       package_dir = {'pywarpx':'pywarpx'},
       description = """Wrapper of WarpX""",
       ext_modules = [example_module],
       )
