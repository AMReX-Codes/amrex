#!/usr/bin/env python

"""
setup.py file for WarpX
"""

import os
from distutils.core import setup, Extension

import numpy

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

boxlib_home = os.environ.get('BOXLIB_HOME', '../../BoxLib')
boxlib_includes = ['Src/C_BaseLib', 
                   'Src/C_ParticleLib',
                   'Src/C_BoundaryLib',
                   'Src/C_AmrCoreLib',
                   'Tools/C_scripts']
boxlib_includes = [os.path.join(boxlib_home, ii) for ii in boxlib_includes]

include_dirs = [numpy_include, '../Source'] + boxlib_includes

example_module = Extension('pywarpx._warpxC',
                           swig_opts=['-c++', '-outdir', 'pywarpx'],
                           sources=['warpxC.i'],
                           library_dirs=['.'],
                           libraries=['warpx'],
                           include_dirs = include_dirs,
                           define_macros = [('BL_USE_MPI','1'), ('BL_SPACEDIM','3'), ('BL_FORT_USE_UNDERSCORE','1'), ('USE_PARTICLES', None)],
                           extra_compile_args = ['-std=c++11'],
                           )

setup (name = 'pywarpx',
       packages = ['pywarpx'],
       package_dir = {'pywarpx':'pywarpx'},
       description = """Wrapper of WarpX""",
       ext_modules = [example_module],
       )
