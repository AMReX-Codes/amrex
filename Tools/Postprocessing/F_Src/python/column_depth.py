#!/usr/bin/env python
import sys
import numpy

def countValidLines(iterableObject):
    """ Here we count the number of lines in the iteratableObject that don't
    start with a comment (#) character.  The iterableObject can be a file,
    in which case we rewind the file to the beginning."""
    import types

    numGoodLines = 0
    for line in iterableObject:
        if not line.lstrip().startswith("#"): numGoodLines +=1

    # if this is a file, rewind it
    if type(iterableObject) is types.FileType: iterableObject.seek(0)
    return numGoodLines

def columnDepth(inputFile,columnDepthStart):
    """ This method takes in a file specifying the density as a function
    of radius and calculates the column depth as a function of radius,
    starting with some specified initial column depth.  That is, we
    integrate some profile, and output the integrated profile to some
    output file.

    The input format is given as :

    # some headers stuff
    # more headers
    # any number of header lines
    <radius_values> <density_values>

    The output profile will be given at the same coordinates as the
    input file."""

    if inputFile:
        try: 
            fh = open(inputFile,'r')
            nLines = countValidLines(fh)
            data = fh.read()
            fh.close()
        except IOError:
            print 'Could not open file %s for reading.' % inputFile
            sys.exit(1)
    else:
        data = sys.stdin.readlines()
        nLines = countValidLines(data)
        
    print nLines
    print data[0]
        
    r = numpy.zeros(nLines,numpy.float64)
    density = numpy.zeros(nLines,numpy.float64)
    columnDepth = numpy.zeros(nLines,numpy.float64)

    columnDepth[-1] = columnDepthStart

    pt = 0
    for line in data:
        if line.lstrip().startswith("#"): continue
        line = line.split()
        r[pt] = line[0]
        density[pt] = line[1]

        pt += 1

    # trapezoidal integration
    for pt in range(nLines-2,-1,-1):
        dr = r[pt]-r[pt+1]

        columnDepth[pt] = -0.5e0*dr*(density[pt]+density[pt+1]) + \
            columnDepth[pt+1]

    for pt in range(nLines):
        print "%s %s" % (r[pt],columnDepth[pt])

if __name__ == "__main__":
    from optparse import OptionParser
    # this should be updated to something more recent, like argparse
    parser = OptionParser()
    parser.add_option("-f","--file", dest="file",
                      help="input FILE containing density profile; if not specified, use stdin",
                      default=None,
                      metavar="FILE")
    parser.add_option("-s","--start", dest="start",
                      default=0.0e0,
                      help="starting value of column depth; defaults to %default")


    (options, args) = parser.parse_args()

    columnDepth(options.file,options.start)
