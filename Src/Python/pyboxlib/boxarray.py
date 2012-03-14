"""PyBoxLib boxarray class."""

import base

from pybl import bl


class boxarray(base.BLObject):
  """BoxArray."""

  def create(self, boxes=[]):
    """Create a boxarray from a list of boxes."""

    self.oid = fboxlib.create_boxarray_from_boxes(boxes)

  def from_regrid(self, lmultifab, buffer_width=0):
    """Create a new boxarray from a logical multifab that tags cells
    that should be refined."""

    self.oid = fboxlib.regrid(lmultifab.oid, buffer_width)
