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


def get_particle_data(species_number, start_comp, num_comp):
    '''

    Return a numpy array containing the particle data
    for the given species. The attributes returned will start
    at start comp, and num_comp attributes per particle will
    be returned.

    '''
    num_particles = libwarpx.warpx_getNumParticles(species_number)
    f = libwarpx.warpx_getParticleData
    f.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=2, 
                                       shape=(num_particles, num_comp), flags="OWNDATA")
    return libwarpx.warpx_getParticleData(species_number, start_comp, num_comp)


def get_particle_ids(species_number):
    '''

    Return a numpy array containing the particle id numbers

    '''
    num_particles = libwarpx.warpx_getNumParticles(species_number)
    f = libwarpx.warpx_getParticleIDs
    f.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=1, 
                                       shape=num_particles, flags="OWNDATA")
    return libwarpx.warpx_getParticleIDs(species_number)


def get_particle_cpu(species_number):
    '''

    Return a numpy array containing the particle cpu numbers

    '''
    num_particles = libwarpx.warpx_getNumParticles(species_number)
    f = libwarpx.warpx_getParticleCPU
    f.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=1, 
                                       shape=num_particles, flags="OWNDATA")
    return libwarpx.warpx_getParticleCPU(species_number)

def get_electric_field(level, direction):
    '''

    Returns a list of numpy arrays containing the electric field data
    on each grid on the specified level.

    '''

    f = libwarpx.warpx_getEfield
    f.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double) )
    shapes = ctypes.POINTER(ctypes.c_int)()
    size = ctypes.c_int(0)
    data = libwarpx.warpx_getEfield(0, 0, ctypes.byref(size), ctypes.byref(shapes))

    grid_data = []
    for i in range(size.value):
        shape=(shapes[3*i+0], shapes[3*i+1], shapes[3*i+2])
        grid_data.append(np.ctypeslib.as_array(data[i], shape))

    del shapes
    del data
    return grid_data


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

data = get_particle_data(0, 0, 2)
print(data)
print(data.shape)

print(get_particle_ids(0))
print(get_particle_cpu(0))

grid_data = get_electric_field(0, 0)

plt.clf()
plt.pcolormesh(grid_data[1][8,:,:])
plt.savefig("field.png")

libwarpx.warpx_finalize()

libwarpx.amrex_finalize()
