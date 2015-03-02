#!/usr/bin/env python
#
# a script showing how to use the helmeos module
# it reads T, rho, X data from a sample data file, calculates abar and zbar
# then calls the eos to get the various thermodynamic quantities

import sys
import numpy
import helmeos

# particular to this data file
lineSkip = 7


###############################################################################
#  NETWORK STUFF:
# eventually this should be automated with a selection of network and f2py
nspec = 9

aion = numpy.zeros(nspec,numpy.float64)
zion = numpy.zeros(nspec,numpy.float64)

aion[0] = 1.0
aion[1] = 3.0
aion[2] = 4.0
aion[3] = 12.0
aion[4] = 14.0
aion[5] = 16.0
aion[6] = 20.0
aion[7] = 24.0
aion[8] = 56.0

zion[0] = 1.0
zion[1] = 2.0
zion[2] = 2.0
zion[3] = 6.0
zion[4] = 7.0
zion[5] = 8.0
zion[6] = 10.0
zion[7] = 12.0
zion[8] = 28.0
###############################################################################

def main(inputFile):

    # make sure we have the python wrappers built around the fortran code
    #
    # this will build the shared object file, fhelmEOS.so, which contains
    # links to the fortran helmeos.f90 code.  
    # the helmeos python module (imported at the top of this script) uses
    # the fhelmEOS.so information
    helmeos.fwrap()

    # read in some data
    fh = open(inputFile,'r')

    X = numpy.zeros(nspec,numpy.float64)
    
    lineno = 0
    for line in fh:
        lineno += 1
        if lineno < lineSkip: continue
        line = line.split()
        
        line.pop(0)
        line.pop(0)

        temp = float(line.pop(0))
        dens = float(line.pop(0))
        X[:] = line

        break

    fh.close()

    # calculate abar and zbar
    #
    # eventually this should be method'ified
    abar = 1.0/numpy.sum(X/aion)
    zbar = numpy.sum(X*zion/aion)*abar

    # dump some output
    print 'temp: %s' % temp
    print 'dens: %s' % dens
    print 'X:', X, type(X)
    print 'sum X: %s' % numpy.sum(X)
    print 'abar: %s' % abar
    print 'zbar: %s' % zbar

    # call the eos
    #
    # the data is stored in an helmeos.EOS_dataType class object that
    # contains the data as attributes (e.g pressure is stored in
    # helmeos.EOS_dataType.p)
    data = helmeos.eos(helmeos.input_rt,abar,zbar,dens,temp)

    # print all of the attributes of the helmeos.EOS_dataType
    #
    # this uses the overwritten helmeos.EOS_dataType.__str__() method
    #    essentially this prints things like "p = <data.p>"
    print data

if __name__ == "__main__":

    if len(sys.argv) == 1:
        print "Usage: test.py <inputFile>"
        print "     a sample <inputFile> is provided in eos_data.txt"
        sys.exit(1)

    main(sys.argv[1])
