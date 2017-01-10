from .Bucket import Bucket
from . import warpxC

class WarpX(Bucket):

    def init(self):
        warpxC.warpx_init()

    def evolve(self, nsteps=None):
        if nsteps is None:
            warpxC.warpx_evolve()
        else:
            warpxC.warpx_evolve(nsteps)

    def finalize(self):
        warpxC.warpx_finalize()

warpx = WarpX('warpx')
