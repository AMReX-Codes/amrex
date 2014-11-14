"""PyBoxLib multifab class."""

import numpy as np
import fboxlib.fcboxlib as fcboxlib

from fboxlib.fab import fab
from fboxlib.layout import layout

class multifab():
  """MultiFAB."""

  def __init__(self, layout=None, nc=1, ng=0, cptr=None):
    if cptr:
      self.cptr = cptr
    else:
      self.cptr = fcboxlib.multifab_create_from_layout(layout.cptr, nc, ng)

  def echo(self):
    fcboxlib.multifab_print(self.cptr)

  # def copy(self, dest):
  #   """Copy self into *dest* (multifab)."""
  #   bl.pybl_multifab_copy(dest.cptr, self.cptr)

  def fill_boundary(self):
    """Fill ghost cells between processors."""
    fcboxlib.multifab_fill_boundary(self.cptr)

  def fab(self, i):
    """Return fab corresponding to box *i*."""
    return fab(self, i)

  @property
  def layout(self):
    return layout(cptr=fcboxlib.multifab_layout(self.cptr))

  @property
  def dim(self):
    d, _, _, _ = fcboxlib.multifab_info(self.cptr)
    return d

  @property
  def nboxes(self):
    _, d, _, _ = fcboxlib.multifab_info(self.cptr)
    return d

  @property
  def ncomp(self):
    _, _, d, _ = fcboxlib.multifab_info(self.cptr)
    return d

  @property
  def nghost(self):
    _, _, _, d = fcboxlib.multifab_info(self.cptr)
    return d

  def write(self, dirname, header):
    fcboxlib.multifab_write(self.cptr,
                            dirname, len(dirname),
                            header, len(header))

  # def read(self, dirname, header):
  #   bl.pybl_multifab_read(dirname, len(dirname), header, len(header),
  #                         byref(self.cptr))
  #   self.get_info()


class lmultifab():
  """Logical MultiFAB."""

  def __init__(self, layout=None, cptr=None):
    if cptr:
      self.cptr = cptr
    else:
      self.cptr = fcboxlib.lmultifab_create_from_layout(layout.cptr)

  def echo(self):
    fcboxlib.lmultifab_print(self.cptr)

  def fab(self, i):
    """Return fab corresponding to box *i*."""
    return fab(self, i, logical=True, squeeze=True)

  @property
  def dim(self):
    d, _, _, _ = fcboxlib.lmultifab_info(self.cptr)
    return d

  @property
  def nboxes(self):
    _, d, _, _ = fcboxlib.lmultifab_info(self.cptr)
    return d

  @property
  def ncomp(self):
    _, _, d, _ = fcboxlib.lmultifab_info(self.cptr)
    return d

  @property
  def nghost(self):
    _, _, _, d = fcboxlib.lmultifab_info(self.cptr)
    return d
