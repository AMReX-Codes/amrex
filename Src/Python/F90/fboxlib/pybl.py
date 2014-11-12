"""PyBoxLib prototypes and misc wrappers."""

__all__ = [ 'bl', 'open', 'mpi_size', 'mpi_rank', 'close' ]

###############################################################################
# load libpyfboxlib, set prototypes etc

import os

from ctypes import *

#bl = CDLL(os.path.dirname(__file__) + '/libboxlib.so')
bl = CDLL('/Users/mwemmett/opt/lib/libboxlib.so')
c_int_p = POINTER(c_int)
c_double_p = POINTER(c_double)

bl.pybl_create_multifab_from_layout.restype = None
bl.pybl_create_multifab_from_layout.argtypes = [
    c_void_p, c_int, c_int, c_void_p ]

bl.pybl_get_multifab_info.restype = None
bl.pybl_get_multifab_info.argtypes = [
    c_void_p, c_int_p, c_int_p, c_int_p, c_int_p ]

bl.pybl_create_layout_from_boxes.restype = None
bl.pybl_create_layout_from_boxes.argtypes = [
    c_void_p, c_int, c_int, c_void_p, c_void_p ]

bl.pybl_multifab_write.argtypes = [
    c_void_p, c_char_p, c_int, c_char_p, c_int ]

bl.pybl_multifab_read.argtypes = [
    c_char_p, c_int, c_char_p, c_int, c_void_p ]

bl.pybl_mpi_reduce_max.argtypes = [
    c_double_p ]

bl.pybl_boxarray_maxsize.argtypes = [
    c_int_p, c_int, c_void_p ]

# bl.multifab_as_numpy.argtypes = [
#     c_void_p, c_int ]

#bl.multifab_as_numpy.restype = c_double_p

c_int3 = 3*c_int
class box(Structure):
    _fields_ = [ ('dim', c_int),
                 ('lo', c_int3),
                 ('hi', c_int3) ]

###############################################################################
# misc routines

def open():
    bl.pybl_open()

def mpi_size():
    size = c_int()
    bl.pybl_mpi_size(byref(size))
    return size.value

def mpi_rank():
    rank = c_int()
    bl.pybl_mpi_rank(byref(rank))
    return rank.value

def mpi_reduce_max(l):
    g = c_double(float(l))
    bl.pybl_mpi_reduce_max(byref(g))
    return g.value

def close():
    bl.pybl_close()
