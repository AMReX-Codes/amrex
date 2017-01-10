#!/usr/bin/env python

"""
setup.py file for WarpX
"""

from distutils.core import setup, Extension


example_module = Extension('pywarpx._warpxC',
                           swig_opts=['-outdir','pywarpx'],
                           sources=['warpxC.i'],
                           library_dirs=['.'],
                           libraries=['warpx'],
                           )

setup (name = 'pywarpx',
       packages = ['pywarpx'],
       package_dir = {'pywarpx':'pywarpx'},
       description = """Wrapper of WarpX""",
       ext_modules = [example_module],
       )
