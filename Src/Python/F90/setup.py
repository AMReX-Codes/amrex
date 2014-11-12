'''PyBoxLib setup script.'''

import os
import sys

from distutils.core import setup
from distutils.util import get_platform
from distutils.extension import Extension
from distutils.command.install import install
from distutils.command.build import build
from distutils.spawn import find_executable

from subprocess import call

import numpy as np


class build_boxlib(build):
    def run(self):
        build.run(self)
        print 'running build_boxlib'
        self.mkpath(self.build_lib)
        target_files = [ 'libboxlib.so' ]
        if not self.dry_run:
            for target in target_files:
                self.copy_file('build/lib/libboxlib.so',
                               os.path.join(self.build_lib, 'fboxlib'))

# build_temp = 'build/temp.%s-%s' % (get_platform(), sys.version[0:3])
blnpy = Extension('fboxlib.blnpy',
                  sources = ['src/boxlib_numpy_c.c'],
                  include_dirs = [np.get_include()],
                  library_dirs = ['build/lib'],
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
