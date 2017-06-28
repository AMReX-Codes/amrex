from .Bucket import Bucket

from .WarpX import warpx
from .Amr import amr
from .Geometry import geometry
from .Algo import algo
from .Langmuirwave import langmuirwave
from .Interpolation import interpolation
from . import Particles
from .Particles import particles, particles_list

import ctypes
from ._libwarpx import libwarpx
from ._libwarpx import amrex_init

class AMReX(object):

    def init(self):
        argv = ['warpx']
        argv += warpx.attrlist()
        argv += amr.attrlist()
        argv += geometry.attrlist()
        argv += algo.attrlist()
        argv += langmuirwave.attrlist()
        argv += interpolation.attrlist()
        argv += particles.attrlist()

        if not particles_list:
            # --- This is needed in case only species_names has been set,
            # --- assuming that only the built in particle types are being used.
            for pstring in particles.species_names.split(' '):
                particles_list.append(getattr(Particles, pstring))

        for particle in particles_list:
            argv += particle.attrlist()

        amrex_init(argv)

    def finalize(self, finalize_mpi=1):
        libwarpx.amrex_finalize(finalize_mpi)
