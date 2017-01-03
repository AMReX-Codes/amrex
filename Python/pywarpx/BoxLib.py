from .Bucket import Bucket

from .WarpX import warpx
from .Amr import amr
from .Geometry import geometry
from .Algo import algo
from .Langmuirwave import langmuirwave
from .Interpolation import interpolation

from . import _warpxC

class BoxLib(object):

    def init(self):
        argv = []
        argv += warpx.attrlist()
        argv += amr.attrlist()
        argv += geometry.attrlist()
        argv += algo.attrlist()
        argv += langmuirwave.attrlist()
        argv += interpolation.attrlist()

        _warpxC.boxlib_init(argv)

    def finalize(self, finalize_mpi=1):
        _warpxC.boxlib_finalize(finalize_mpi)
