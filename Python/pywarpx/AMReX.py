from .Bucket import Bucket

from .WarpX import warpx
from .Amr import amr
from .Geometry import geometry
from .Algo import algo
from .Langmuirwave import langmuirwave
from .Interpolation import interpolation
from .Particles import particles, electrons

import ctypes
from ._libwarpx import libwarpx
from ._libwarpx import amrex_init

class AMReX(object):

    def init(self):
        argv = []
        argv += warpx.attrlist()
        argv += amr.attrlist()
        argv += geometry.attrlist()
        argv += algo.attrlist()
        argv += langmuirwave.attrlist()
        argv += interpolation.attrlist()
        argv += particles.attrlist()
        argv += electrons.attrlist()

        amrex_init(argv)

    def finalize(self, finalize_mpi=1):
        libwarpx.amrex_finalize(finalize_mpi)
