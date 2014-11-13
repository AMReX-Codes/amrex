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

  def echo(self):
    fcboxlib.layout_print(self.cptr)
