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
            if use_mpi:
                cc = os.environ.get('CC', 'mpicc')
                cxx = os.environ.get('CXX', 'mpic++')
                mpihome = os.environ.get('MPIHOME', None)
                if mpihome is None:
                    mpicc   = find_executable('mpicc')
                    if mpicc is None:
                        raise ValueError("'mpicc' not found.  Please install MPI so that 'mpicc' and 'mpicxx' are in your PATH, or set MPIHOME appropriately.")
                    mpihome = os.path.dirname(os.path.dirname(mpicc))
                call([ 'make', 'MPI_HOME=' + mpihome, 'CC=' + cc, 'CXX=' + cxx, 'OUT=' + self.build_temp ])
            else:
                call([ 'make', 'USE_MPI=FALSE', 'OUT=' + self.build_temp ])
            print '*' * 80

        self.execute(compile, [], 'compiling boxlib')

        self.mkpath(self.build_lib)
        target_files = [ '_bl1.so', '_bl2.so', '_bl3.so' ]
        if not self.dry_run:
            for target in target_files:
                self.copy_file(os.path.join(self.build_temp, target),
                               os.path.join(self.build_lib, 'boxlib'))


setup(

    name         = "PyBoxLib",
    packages     = ['boxlib'],
    author       = "Matthew Emmett and Marc Day",
    author_email = "mwemmett@lbl.gov",
    description  = "Python wrappers for BoxLib.",
    license      = "XXX",
    keywords     = "BoxLib",
    url          = "https://ccse.lbl.gov/BoxLib/",

    cmdclass     = {
        'build': build_boxlib,
        }
    )
