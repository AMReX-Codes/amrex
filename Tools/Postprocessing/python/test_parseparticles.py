#!/usr/bin/env python

# simple script showing how to make plots of particles using the parseparticles
# module
#

import parseparticles
import pylab

def main(fileList):
    """
    This script shows how to make plots using the parseparticles module.

    usage: ./test.py particleFile1 [particleFile2 particleFile3 ... ]
    """

    # this returns a dict whose keys are a unique identifier (based on 
    # id and CPU) and values are the actual particle objects
    particlesDict = parseparticles.parseParticleFile(fileList)

    # get just the particle objects
    particles = particlesDict.values()

    print "there are %d unique particles" % len(particles)

    # plots - this is based on the plotting done in the original 
    # parseparticles.py script, which has since become a module for importing
    pylab.clf()

    # global domain extrema
    xmin = 1.e33
    xmax = -1e33
    ymin = 1.e33
    ymax = -1.e33

    # loop over the particle objects and create the particle paths
    for particle in particles:
        coords, time = particle.flatten()

        # find the extrema of the particles over time to make all frames the
        # same size
        xmin = min(xmin,min(coords[0,:]))
        xmax = max(xmax,max(coords[0,:]))
        ymin = min(ymin,min(coords[1,:]))
        ymax = max(ymax,max(coords[1,:]))

        pylab.scatter([coords[0,0]], [coords[1,0]], marker="x")
        pylab.plot(coords[0,:], coords[1,:])

    pylab.xlabel("x")
    pylab.ylabel("y")

    pylab.savefig("particle_paths.png")

    # make an animation -- note: this assumes that all particles exist at all
    # timesteps


    # grab the total number of timesteps from the first particle's history
    totalSteps = len(particles[0].history)

    # nstep = 0
    # while (nstep < totalSteps):
    #     pylab.clf()

    #     for particle in particles:
    #         pylab.scatter([particle.history[nstep].xyz[0]], 
    #                       [particle.history[nstep].xyz[1]],
    #                       marker="o", s=1.0, edgecolor="None")

    #     pylab.xlabel("x")
    #     pylab.ylabel("y")

    #     pylab.axis([xmin,xmax,ymin,ymax])
    #     a = pylab.gca()
    #     a.set_aspect("equal")

    #     pylab.savefig("particles_%04d.png" % nstep)

    #     nstep += 1

if __name__ == "__main__":
    import sys
#    import cProfile

    if len(sys.argv) == 1: 
        print main.__doc__
        sys.exit()

    # pass in the particle files as arguments to this function
    main(sys.argv[1:])

# this is for profiling
#    cProfile.run("main(sys.argv[1:])","profile.tmp2")

