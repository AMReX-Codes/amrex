#!/usr/bin/env python

"""
setup.py file for WarpX
"""

import os
from distutils.core import setup, Extension
import platform
import numpy

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

amrex_home = os.environ.get('AMREX_HOME', '../../amrex')
amrex_includes = ['Src/Base', 
                   'Src/Particle',
                   'Src/Boundary',
                   'Src/AmrCore',
                   'Tools/scripts']
amrex_includes = [os.path.join(amrex_home, ii) for ii in amrex_includes]

include_dirs = [numpy_include, '../Source'] + amrex_includes

definesstring = os.environ.get('DEFINES','')
defines = definesstring.split(' ')

#cpp11_flags = [] #['-std=c++11']
#if platform.system() == "Darwin":
#    macosx_deployment_target = platform.mac_ver()[0]
#    os.environ['MACOSX_DEPLOYMENT_TARGET'] = macosx_deployment_target
#    cpp11_flags.append("-stdlib=libc++")

example_module = Extension('pywarpx._warpxC',
                           swig_opts=['-c++', '-outdir', 'pywarpx'] + defines,
                           sources=['warpxC.i'],
                           library_dirs=['.'],
                           libraries=['warpx'],
                           include_dirs = include_dirs,
                           #define_macros = [('BL_USE_MPI','1'), ('BL_SPACEDIM','3'), ('BL_FORT_USE_UNDERSCORE','1'), ('USE_PARTICLES', None)],
                           #extra_compile_args = cpp11_flags,
                           )

setup (name = 'pywarpx',
       packages = ['pywarpx'],
       package_dir = {'pywarpx':'pywarpx'},
       description = """Wrapper of WarpX""",
       ext_modules = [example_module],
       )
