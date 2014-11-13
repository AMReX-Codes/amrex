"""PyBoxLib boxarray class."""

import base
from ctypes import *


from pybl import bl


class boxarray(base.BLObject):
  """BoxArray."""

  def create(self, boxes):
    """Create a boxarray from a list of boxes."""

    nboxes = len(boxes)
    dim    = len(boxes[0][0])

    cboxes = (2*dim*nboxes*c_int)()
    for b, box in enumerate(boxes):
      for t in (0, 1):
        for c, comp in enumerate(boxes[b][t]):
          cboxes[c*2*nboxes+t*nboxes+b] = boxes[b][t][c]

    bl.pybl_create_boxarray_from_boxes(cboxes, nboxes, dim, byref(self.cptr))

    self.dim = dim

  def maxsize(self, size):
    dim = len(size)
    csizes = (dim*c_int)()
    for i, s in enumerate(size):
      csizes[i] = s
    bl.pybl_boxarray_maxsize(csizes, dim, self.cptr)

  @property
  def nboxes(self):
    nb = c_int()
    bl.pybl_boxarray_nboxes(self.cptr, byref(nb))
    return nb.value


  # def from_regrid(self, lmultifab, buffer_width=0):
  #   """Create a new boxarray from a logical multifab that tags cells
  #   that should be refined."""

  #   self.oid = bl.regrid(lmultifab.oid, buffer_width)
