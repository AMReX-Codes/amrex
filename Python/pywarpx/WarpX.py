from .Bucket import Bucket
from . import _warpxC

class WarpX(Bucket):

    def init(self):
        _warpxC.warpx_init()

    def evolve(self, nsteps=None):
        if nsteps is None:
            _warpxC.warpx_evolve()
        else:
            _warpxC.warpx_evolve(nsteps)

    def finalize(self):
        _warpxC.warpx_finalize()

warpx = WarpX('warpx')
