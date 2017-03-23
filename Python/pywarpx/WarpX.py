from .Bucket import Bucket
from ._libwarpx import libwarpx

class WarpX(Bucket):
    """
    A Python wrapper for the WarpX C++ class
    """

    def init(self):
        libwarpx.warpx_init()

    def evolve(self, nsteps=-1):
        libwarpx.warpx_evolve(nsteps)

    def finalize(self):
        libwarpx.warpx_finalize()

    def getProbLo(self, direction):
        return libwarpx.warpx_getProbLo(direction)

    def getProbHi(self, direction):
        return libwarpx.warpx_getProbHi(direction)

warpx = WarpX('warpx')
