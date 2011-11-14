"""Base class for interacting with the BoxLib/Python object store.

See also src/blobjects.py.
"""

from pyfboxlib import fboxlib

class BLObject(object):

  def __init__(self):
    self._oid = None

  def __str__(self):
    return 'pyboxlib type: %s, oid: %d' % (self.__class__.__name__, self.oid)

  # XXX: delete method

  def echo(self):
    pr = 'print_' + self.__class__.__name__
    pr = getattr(fboxlib, pr, None)

    if pr:
      pr(self.oid)

  @property
  def oid(self):
    if self._oid is None:
      raise ValueError('PyBoxLib: OID not associated')
    return self._oid

  @oid.setter
  def oid(self, value):
    assert self._oid is None
    assert isinstance(value, int)
    self._oid = value
