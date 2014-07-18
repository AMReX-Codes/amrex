"""PyBoxLib.

Please see the README file in BoxLib/Src/Python for a technical
overview.
"""

import numpy as np

class BOXLIB(object):
    """Dimensionally aware BoxLib import wrapper."""
    def __init__(self):
        self.mod = None
    def __getitem__(self, dim):
        if self.mod is None:
            self.dim = dim
            self.mod = __import__('bl%d' % dim, globals(), locals(), [], -1)
            self.mod.StartParallel()
        else:
            if dim > 0:
                assert(dim == self.dim)
        return self.mod

bl = BOXLIB()

def rank():
    return bl[0].rank()

def size():
    return bl[0].size()

def Box(lo=[], hi=[]):
    """Dimensionally aware Box wrapper."""
    dim = len(lo)
    lo  = bl[dim].IntVect(*lo)
    hi  = bl[dim].IntVect(*hi)
    bx  = bl[dim].Box(lo, hi)
    return bx

def RealBox(lo=[], hi=[]):
    """Dimensionally aware Box wrapper."""
    dim = len(lo)
    bx  = bl[dim].RealBox(np.asarray(lo), np.asarray(hi))
    return bx

def lo(bx):
    lo = bx.smallEnd()
    return [ lo[x] for x in range(len(lo)) ]

def hi(bx):
    hi = bx.bigEnd()
    return [ hi[x] for x in range(len(hi)) ]

def BoxArray(boxes):
    """Dimensionally aware BoxArray wrapper."""
    dim = len(boxes[0].size())
    ba  = bl[dim].BoxArray()
    ba.resize(len(boxes))
    for i, bx in enumerate(boxes):
        ba.set(i, bx)
    return ba

def MultiFab(ba, ncomp=1, nghost=0):
    """Dimensionally aware MultiFab wrapper."""
    dim = len(ba.get(0).size())
    return bl[dim].MultiFab(ba, ncomp, nghost)

