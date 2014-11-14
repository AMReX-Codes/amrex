"""PyBoxLib layout class."""

import numpy as np
import fboxlib.fcboxlib as fcboxlib

class layout():
  """BoxLib layout."""

  def __init__(self, boxarray=None, cptr=None):
    if cptr:
      self.cptr = cptr
    else:
      pmask = np.ones(boxarray.dim, np.int32)
      self.cptr = fcboxlib.layout_create_from_boxarray(boxarray.cptr, pmask)

  @property
  def nboxes(self):
    """Return number of boxes."""
    return fcboxlib.layout_nboxes(self.cptr)

  @property
  def local_boxes(self):
    """List of local boxes."""
    return [ n for n in range(1, self.nboxes+1) if self.is_local(n) ]

  def is_local(self, n):
    return fcboxlib.layout_local(self.cptr, n)

  def get_box(self, n):
    if 1 <= n <= self.nboxes:
      return fcboxlib.layout_get_box(self.cptr, n)
    return None

  def echo(self):
    fcboxlib.layout_print(self.cptr)
