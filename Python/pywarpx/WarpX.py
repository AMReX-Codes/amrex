from .Bucket import Bucket
from . import warpxC

class WarpX(Bucket):
    """
    A Python wrapper for the WarpX C++ class
    """

    def init(self):
        self.warpx = warpxC.WarpX.GetInstance()
        self.warpx.InitData()

    def evolve(self, nsteps=-1):
        self.warpx.Evolve(nsteps)

    def finalize(self):
        warpxC.WarpX.ResetInstance()

    def getProbLo(self, direction):
        return self.warpx.Geom()[0].ProbLo(direction)
        #return warpxC.warpx_getProbLo(direction)

    def getProbHi(self, direction):
        return self.warpx.Geom()[0].ProbHi(direction)
        #return warpxC.warpx_getProbHi(direction)

warpx = WarpX('warpx')
