#!/usr/bin/env python

import matplotlib.pyplot as plt
import warpx

with warpx.script():

    # run for ten time steps
    warpx.evolve(10)

    x = warpx.get_particle_x(0)
    y = warpx.get_particle_y(0)
    plt.plot(x[0], y[0], '.')
    plt.savefig("particles.png")

    # this returns a list of numpy arrays that hold the magnetic field 
    # data in the x-direction on each grid for level 0
    grid_data = warpx.get_mesh_magnetic_field(0, 0, False)

    # plot a slice through the second grid 
    plt.clf()
    plt.pcolormesh(grid_data[1][9,:,:])
    plt.savefig("field.png")

