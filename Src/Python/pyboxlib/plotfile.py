"""PyBoxLib plotfile class."""

from ctypes import *

import numpy as np
import base

from pybl import bl
from fab import fab

class plotfile(base.BLObject):
  """PlotFile."""

  def __init__(self, *args, **kwargs):

    super(plotfile, self).__init__(*args, **kwargs)

    self.c_int_attrs = [ 'dim', 'nvars', 'flevel' ]
    self.ini_c_int_attrs()


  def create(self, dname):
    """Create a plotfile given a name."""

    create = getattr(bl, 'pybl_create_plotfile')
    create(dname, len(dname), byref(self.cptr))

    self.get_info()


  def get_info(self):

    if self.associated:
      get_info = getattr(bl, 'pybl_get_plotfile_info')
      get_info(self.cptr,
               byref(self._dim),
               byref(self._nvars),
               byref(self._flevel))
