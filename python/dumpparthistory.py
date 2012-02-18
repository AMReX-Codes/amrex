#!/usr/bin/env python

# a simple routine to parse particle files and dump out the particle
# histories into separate files (1 file per particle) so that they can
# be plotted elsewhere

#import math
import numpy
import sys
import parseparticles


#-----------------------------------------------------------------------------
def main(files):

    # this returns a dict whose keys are a unique identifier (based on 
    # id and CPU) and values are the actual particle objects
    particlesDict = parseparticles.parseParticleFile(files)

    # get just the particle objects
    particles = particlesDict.values()


    # print out some info and test the routines
    print "number of unique particles = ", len(particles)


    # output particle instance info (debugging)
    #n = 0
    #while (n < len(particles)):
    #    print n
    #    print particles[n]
    #    print " "
    #    n += 1


    # dump out the data, particle-by-particle
    n = 0
    while (n < len(particles)):

        # get numpy arrays containing the time and coordinate
        # information for particle 0
        coords, time = particles[n].flatten()
        dim = particles[n].dim 

        # output to a flie
        of = open("particle_history.%03d" % (n), 'w')

        # header
        #   first the particle id information
        idstuff = str(particles[n])
        for line in idstuff.split("\n"):
            of.write("# %s\n" % (line))
        
        #   next the column labels
        if (dim == 1):
            of.write("# %20s %20s\n" % ("time", "x"))            
        elif (dim == 2):
            of.write("# %20s %20s %20s\n" % ("time", "x", "y"))            
        elif (dim == 3):
            of.write("# %20s %20s %20s %20s\n" % ("time", "x", "y", "z"))
            
        
        # t, x, [y, [z]] in columns
        i = 0
        while (i < len(particles[n].history)):
            
            if (dim == 1):
                of.write("  %20.10f %20.10f\n" % 
                         (time[i], coords[0,i]))
            elif (dim == 2):
                of.write("  %20.10f %20.10f %20.10f\n" % 
                         (time[i], coords[0,i], coords[1,i]))
            elif (dim == 3):
                of.write("  %20.10f %20.10f %20.10f %20.10f\n" % 
                         (time[i], coords[0,i], coords[1,i], coords[2,i]))

            i += 1
               

        of.close()


        n += 1



        


#-----------------------------------------------------------------------------
if __name__== "__main__":

    if (len(sys.argv) == 1):
        print "ERROR: no particle data files specified\n"
        sys.exit(2)

    main(sys.argv[1:])




        
