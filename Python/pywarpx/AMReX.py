from .Bucket import Bucket

from .WarpX import warpx
from .Amr import amr
from .Geometry import geometry
from .Algo import algo
from .Langmuirwave import langmuirwave
from .Interpolation import interpolation
from .Laser import laser
from . import Particles
from .Particles import particles, particles_list

import ctypes
from ._libwarpx import libwarpx
from ._libwarpx import amrex_init

class AMReX(object):

    def create_argv_list(self):
        argv = []
        argv += warpx.attrlist()
        argv += amr.attrlist()
        argv += geometry.attrlist()
        argv += algo.attrlist()
        argv += langmuirwave.attrlist()
        argv += interpolation.attrlist()
        argv += particles.attrlist()
        argv += laser.attrlist()

        if not particles_list:
            # --- This is needed in case only species_names has been set,
            # --- assuming that only the built in particle types are being used.
            for pstring in particles.species_names.split(' '):
                particles_list.append(getattr(Particles, pstring))

        for particle in particles_list:
            argv += particle.attrlist()

        return argv

    def init(self):
        argv = ['warpx'] + self.create_argv_list()
        amrex_init(argv)

    def finalize(self, finalize_mpi=1):
        libwarpx.amrex_finalize(finalize_mpi)

    def write_inputs(self, filename='inputs'):
        argv = self.create_argv_list()
        with open(filename, 'w') as ff:

            for arg in argv:
                ff.write('{0}\n'.format(arg))

