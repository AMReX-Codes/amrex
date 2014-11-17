"""PyBoxLib boxarray class."""

import numpy as np
import fboxlib.fcboxlib as fcboxlib

class boxarray():
  """BoxArray."""

  def __init__(self, boxes):
    """Create a boxarray from a list of boxes."""

    nboxes = len(boxes)
    dim    = len(boxes[0][0])
    cboxes = np.zeros((nboxes,2,dim), np.int32)
    for b, box in enumerate(boxes):
      for t in (0, 1):
        for c, comp in enumerate(boxes[b][t]):
          cboxes[b,t,c] = boxes[b][t][c]
    self.cptr = fcboxlib.boxarray_create_from_boxes(cboxes, nboxes, dim)

  @property
  def nboxes(self):
    return fcboxlib.boxarray_nboxes(self.cptr)

  @property
  def dim(self):
    return fcboxlib.boxarray_dim(self.cptr)

  def maxsize(self, size):
    size = np.asarray(size, np.int32)
    fcboxlib.boxarray_maxsize(self.cptr, size)

  def echo(self):
    fcboxlib.boxarray_print(self.cptr)
