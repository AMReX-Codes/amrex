'''PyBoxLib setup script.'''

import os

from distutils.core import setup
from distutils.command.install import install
from distutils.command.build import build

class build_boxlib(build):
    def run(self):
        build.run(self)
        print 'running build_boxlib'
        self.mkpath(self.build_lib)
        target_files = [ 'fcboxlib.so' ]
        if not self.dry_run:
            self.copy_file('build/fcboxlib.so',
                           os.path.join(self.build_lib, 'fboxlib'))

setup(
    name         = "PyBoxLib",
    packages     = ['fboxlib'],
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
