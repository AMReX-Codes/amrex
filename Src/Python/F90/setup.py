'''PyBoxLib setup script.'''

import glob
import os
import re

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install import install
from distutils.command.build import build
from distutils.spawn import find_executable

from subprocess import call

import numpy as np


class build_boxlib(build):

    user_options = build.user_options + [ ('disable-mpi', None, "disable MPI support") ]
    boolean_options = build.boolean_options + [ 'disable-mpi' ]

    def initialize_options(self):
        build.initialize_options(self)
        self.disable_mpi = 0

    def run(self):
        build.run(self)

        use_mpi = self.disable_mpi == 0

        print 'running build_boxlib'
        self.mkpath(self.build_temp)
        def compile():
            print '*' * 80

            fc  = os.environ.get('FC', 'mpif90')
            cc  = os.environ.get('CC', 'mpicc')
            if use_mpi:
                mpihome = os.environ.get('MPIHOME', None)
                if mpihome is None:
                    mpicc = find_executable('mpicc')
                    if mpicc is None:
                        raise ValueError("'mpicc' not found.  Please install MPI so that 'mpicc' and 'mpicxx' are in your PATH, or set MPIHOME appropriately.")
                    mpihome = os.path.dirname(os.path.dirname(mpicc))

            cmd = [ 'make', 'OUT=' + self.build_temp, 'FC=' + fc, 'CC=' + cc ]
            if use_mpi:
                cmd += [ 'MPIHOME=' + mpihome ]
            else:
                cmd += [ 'USE_MPI=FALSE' ]

            call(cmd)
            print '*' * 80

        self.execute(compile, [], 'compiling boxlib')

        self.mkpath(self.build_lib)
        target_files = [ 'libboxlib.so' ]
        if not self.dry_run:
            for target in target_files:
                self.copy_file(os.path.join(self.build_temp, target),
                               os.path.join(self.build_lib, 'fboxlib'))

blnpy = Extension('blnpy',
                  sources = ['src/boxlib_numpy_c.c'],
                  include_dirs = [np.get_include()],
                  library_dirs = ['/Users/mwemmett/opt/lib'],
                  libraries = ['boxlib'])

setup(
    name         = "PyBoxLib",
    packages     = ['fboxlib'],
    author       = "Matthew Emmett and Marc Day",
    author_email = "mwemmett@lbl.gov",
    description  = "Python wrappers for BoxLib.",
    license      = "XXX",
    keywords     = "BoxLib",
    url          = "https://ccse.lbl.gov/BoxLib/",
    ext_modules  = [ blnpy ],

    cmdclass     = {
        'build': build_boxlib,
        }
    )
