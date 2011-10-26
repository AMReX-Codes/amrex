"""PyBoxLib layout class."""

import numpy as np

from pyfboxlib import fboxlib
from boxarray import boxarray
import base


class layout(base.BLObject):
  """BoxLib layout."""

  def create(self, boxarray=None, boxes=None):
    """Create a layout from a list of boxes."""

    if boxarray is not None and boxes is not None:
      raise ValueError('both boxarray and boxes specified')

    if boxes is not None:
      self.oid = fboxlib.create_layout_from_boxes(boxes)

    if boxarray is not None:
      self.oid = fboxlib.create_layout_from_boxarray(boxarray.oid)

  def from_regrid(self, lmfab):
    """Create a layout from re-gridding a tagged multifab."""

    ba = boxarray()
    ba.from_regrid(lmfab)

    self.create(boxarray=ba)

class mllayout(base.BLObject):
  """BoxLib multi-level layout."""

  def create(self, layouts):
    """Create a multi-level layout from a list of layouts."""

    self.oid = fboxlib.create_ml_layout_from_layouts([ x.oid for x in layouts ])
