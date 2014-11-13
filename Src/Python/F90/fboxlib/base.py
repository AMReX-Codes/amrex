"""Base class for interacting with the BoxLib/Python interface."""

from ctypes import *
from pybl import bl


class BLObject(object):


  def __init__(self, cptr=None):
    if cptr is None:
      self.cptr   = c_void_p(None)
    else:
      self.ctpr = cptr
    self.c_int_attrs = []


  def __str__(self):
    return 'pyboxlib type: %s, cptr: %s' % (self.__class__.__name__, self.cptr)


  def echo(self):
    pr = 'pybl_print_' + self.__class__.__name__
    pr = getattr(bl, pr, None)
    if pr is not None:
      pr(self.cptr)


  @property
  def associated(self):
    return self.cptr.value is not None


  def init_c_int_attrs(self):

    for name in self.c_int_attrs:
      setattr(self, '_' + name, c_int(0))


  def __getattr__(self, name):

    if name in self.c_int_attrs:
      attr = getattr(self, '_' + name)
      return attr.value

    # raise AttributeError()
