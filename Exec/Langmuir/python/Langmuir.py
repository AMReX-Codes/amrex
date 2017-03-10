#!/usr/bin/env python

import sys
import ctypes

libwarpx = ctypes.CDLL("libwarpx.so")

LP_c_char = ctypes.POINTER(ctypes.c_char)
LP_LP_c_char = ctypes.POINTER(LP_c_char)

libwarpx.amrex_init.argtypes = (ctypes.c_int, LP_LP_c_char)

argc = len(sys.argv)
argv = (LP_c_char * (argc+1))()
for i, arg in enumerate(sys.argv):
    enc_arg = arg.encode('utf-8')
    argv[i] = ctypes.create_string_buffer(enc_arg)

libwarpx.amrex_init(argc, argv)

libwarpx.warpx_init()

libwarpx.warpx_evolve(10)

libwarpx.warpx_finalize()

libwarpx.amrex_finalize()
