"""Base class for interacting with the BoxLib/Python interface."""

from ctypes import *
from pybl import bl


class BLObject(object):


  def __init__(self):
    self.cptr = c_void_p(None)


  def __str__(self):
    return 'pyboxlib type: %s, cptr: %s' % (self.__class__.__name__, self.cptr)


  def echo(self):
    pr = 'pybl_print_' + self.__class__.__name__
    pr = getattr(bl, pr, None)
    pr(self.cptr)


  @property
  def associated(self):
    return self.cptr.value is not None
