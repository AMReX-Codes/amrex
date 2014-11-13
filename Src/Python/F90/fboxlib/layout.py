"""PyBoxLib layout class."""

from ctypes import *

import base

from boxarray import boxarray
from pybl import bl, box


class layout(base.BLObject):
  """BoxLib layout."""


  def create(self, boxarray=None, boxes=None, is_periodic=True):
    """Create a layout from a list of boxes."""

    # XXX: pmask etc

    if boxarray is not None and boxes is not None:
      raise ValueError('both boxarray and boxes specified')

    if boxes is not None:
      nboxes = len(boxes)
      dim    = len(boxes[0][0])

      # create C boxes
      cboxes = (2*dim*nboxes*c_int)()
      for b, box in enumerate(boxes):
        for t in (0, 1):
          for c, comp in enumerate(boxes[b][t]):
            cboxes[c*2*nboxes+t*nboxes+b] = boxes[b][t][c]

      pmask = (dim*c_int)()
      if not isinstance(is_periodic, list):
        is_periodic = dim * [ is_periodic ]
      for i, periodic in enumerate(is_periodic):
        pmask[i] = 1 if periodic else 0

      # create the layout
      bl.pybl_create_layout_from_boxes(cboxes, nboxes, dim, pmask, byref(self.cptr))

    if boxarray is not None:
      pmask = (boxarray.dim*c_int)()
      for i in range(boxarray.dim):
        pmask[i] = 1
      bl.pybl_create_layout_from_boxarray(boxarray.cptr, boxarray.dim, pmask, byref(self.cptr))

  @property
  def nboxes(self):
    """Return number of boxes."""
    nboxes = c_int()
    bl.pybl_layout_nboxes(self.cptr, byref(nboxes))
    return nboxes.value

  @property
  def local_boxes(self):
    """Return list of local boxes."""
    boxes = []
    for k in range(self.nboxes):
      if self.local(k+1):
        boxes.append(k+1)
    return boxes


  def local(self, box):
    """Return True if box is local."""
    local = c_int()
    bl.pybl_layout_local(self.cptr, box, byref(local))
    return local.value != 0


  def get_box(self, i):
    bx = box()
    bl.pybl_layout_get_box(self.cptr, byref(bx), i)
    return bx


  def from_regrid(self, lmfab):
    """Create a layout from re-gridding a tagged multifab."""

    ba = boxarray()
    ba.from_regrid(lmfab)

    self.create(boxarray=ba)


class mllayout(base.BLObject):
  """BoxLib multi-level layout."""

  def create(self, layouts):
    """Create a multi-level layout from a list of layouts."""

    nlayouts = len(layouts)

    cptrs = (nlayouts*c_void_p)()
    for i, l in enumerate(layouts):
      cptrs[i] = l.cptr

    bl.pybl_create_ml_layout_from_layouts(cptrs, c_int(nlayouts), byref(self.cptr))
