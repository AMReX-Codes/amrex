"""PyBoxLib prototypes and misc wrappers."""

__all__ = [ 'bl', 'open', 'mpi_size', 'mpi_rank', 'close' ]


###############################################################################
# load libpyfboxlib, set prototypes etc

from ctypes import *

bl = CDLL('./libpyfboxlib.so')

c_int_p = POINTER(c_int)

bl.pybl_create_multifab_from_layout.restype = None
bl.pybl_create_multifab_from_layout.argtypes = [
    c_void_p, c_int, c_int, c_int, c_void_p ]

bl.pybl_get_multifab_info.restype = None
bl.pybl_get_multifab_info.argtypes = [
    c_void_p, c_int_p, c_int_p, c_int_p, c_int_p ]

bl.pybl_create_layout_from_boxes.restype = None
bl.pybl_create_layout_from_boxes.argtypes = [
    c_void_p, c_int, c_int, c_void_p ]



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


def close():
    bl.pybl_close()
