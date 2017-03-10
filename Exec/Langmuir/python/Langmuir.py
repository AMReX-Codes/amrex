#!/usr/bin/env python

import sys
import ctypes
import numpy as np
import matplotlib.pyplot as plt

libwarpx = ctypes.CDLL("libwarpx.so")

def get_positions(species_number):
    '''

    Return a numpy array containing the positions of all the
    particles in the given species.

    '''
    num_particles = libwarpx.warpx_getNumParticles(species_number)
    f = libwarpx.warpx_getParticlePositions
    f.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=2, 
                                   shape=(num_particles, 3), flags="OWNDATA")
    return libwarpx.warpx_getParticlePositions(species_number)


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

positions = get_positions(0)
plt.plot(positions[:,0], positions[:,1], '.')
plt.savefig('particles.png')

libwarpx.warpx_evolve(10)

libwarpx.warpx_finalize()

libwarpx.amrex_finalize()
