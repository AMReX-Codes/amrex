#!/usr/bin/env python

# a simple script to plot some select data from a 1-d flame plotfile
#
# 2012-08-08 M. Zingale

import fsnapshot
import numpy
import pylab
import matplotlib
import os
import sys
import getopt
import math
import string
import mpl_toolkits.axes_grid1


speciesPlotThresh = 0.05

class var(object):

    def __init__(self, name, index, vmin, vmax):
        self.name = name
        self.min  = vmin
        self.max  = vmax
        self.index = index

   


#==============================================================================
# do_plot
#==============================================================================
def do_plot(plotfile, outFile, eps, dpi):


    #--------------------------------------------------------------------------
    # construct the output file name
    #--------------------------------------------------------------------------
    if (outFile == ""):
        outFile = os.path.normpath(plotfile) 

        if (not eps):
            outFile += ".png"

        else:
            outFile += ".eps"

    else:
        # make sure the proper extension is used
        if (not eps):
            if (not string.rfind(outFile, ".png") > 0):
                outFile = outFile + ".png"

        else:
            if (not string.rfind(outFile, ".eps") > 0):
                outFile = outFile + ".eps"


    #--------------------------------------------------------------------------
    # read in the meta-data from the plotfile
    #--------------------------------------------------------------------------
    (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile)

    time = fsnapshot.fplotfile_get_time(plotfile)

    (xmin, xmax, ymin, ymax, zmin, zmax) = \
        fsnapshot.fplotfile_get_limits(plotfile)


    #--------------------------------------------------------------------------
    # get the variable information
    #--------------------------------------------------------------------------
    nvar = fsnapshot.fplotfile_get_nvar(plotfile)

    varinfo = []
    
    # note, the variable indices are 1-based, since they live in
    # Fortran land.
    n = 1
    while (n <= nvar):
        (varname, varmin, varmax, ierr) = \
            fsnapshot.fplotfile_get_varinfo(plotfile, n)        

        if (ierr == 0):
            varinfo.append(var(varname.rstrip(), n, varmin, varmax))
        else:
            print "ERROR: invalid variable name"
            sys.exit(2)
        
        n += 1
    

    # n = 0
    # while (n < len(varinfo)):
    #     print varinfo[n].name, varinfo[n].index
    #     n += 1



    #--------------------------------------------------------------------------
    # plot
    #--------------------------------------------------------------------------

    # <<< density and temperature >>>
    pylab.clf()

    sp = pylab.subplot(311)

    sp.set_yscale('log')


    # get the density
    x = numpy.zeros( (nx), dtype=numpy.float64)
    rho = numpy.zeros( (nx), dtype=numpy.float64)

    (rho, x, npts, err) = \
        fsnapshot.fplotfile_get_data_1d(plotfile, "density", rho, x)

    if (not err == 0):
        sys.exit(2)


    # get the temperature
    T = numpy.zeros( (nx), dtype=numpy.float64)

    (T, x, npts, err) = \
        fsnapshot.fplotfile_get_data_1d(plotfile, "tfromp", T, x)

    if (not err == 0):
        sys.exit(2)

    
    # density plot
    pylab.plot(x, T, color="blue")
    pylab.ylabel(r"temperature (K)", color="blue")

    ax = pylab.gca()

    # the offset text is the 1.e8 that appears under the axis labels
    # when doing scientific notation
    ax.tick_params(labeltop='off')
    ax.tick_params(labelbottom='off')
    ax.xaxis.offsetText.set_visible(False)


    # temperature plot
    sp2 = pylab.twinx()
    sp2.set_yscale('log')

    pylab.plot(x, rho, color="red")
    pylab.ylabel(r"density (g cm$^{-3}$)", color="red")



    # <<< species >>>

    sp = pylab.subplot(312)
    
    # find where the species lie
    n = 0
    nSpeciesStart = -1
    nSpeciesEnd = -1

    while (n < len(varinfo)):
        if (varinfo[n].name.find("X(") >= 0 and nSpeciesStart == -1):
            nSpeciesStart = n
            n += 1
            continue

        if (varinfo[n].name.find("X(") == -1 and 
            nSpeciesStart >= 0 and nSpeciesEnd == -1):
            nSpeciesEnd = n-1

        n += 1


    # loop over all the species.  Only plot those that have abundances
    # that exceed the threshold speciesPlotThresh
    X = numpy.zeros( (nx), dtype=numpy.float64)

    n = nSpeciesStart
    while (n <= nSpeciesEnd):

        # is this species abundant enough?
        if (varinfo[n].max < speciesPlotThresh):
            n += 1
            continue


        (X, x, npts, err) = \
            fsnapshot.fplotfile_get_data_1d(plotfile, varinfo[n].name, X, x)

        if (not err == 0):
            sys.exit(2)

        pylab.plot(x, X, label=varinfo[n].name, lw=1)

        n += 1

    leg = pylab.legend(loc=1,labelspacing=0.1)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')

    leg.draw_frame(0)

    ax = pylab.gca()

    # the offset text is the 1.e8 that appears under the axis labels
    # when doing scientific notation
    ax.tick_params(labeltop='off')
    ax.tick_params(labelbottom='off')
    ax.xaxis.offsetText.set_visible(False)


    # <<< enuc >>>
    sp = pylab.subplot(313)

    #sp.set_yscale('log')


    # get enuc
    enuc = numpy.zeros( (nx), dtype=numpy.float64)

    (enuc, x, npts, err) = \
        fsnapshot.fplotfile_get_data_1d(plotfile, "enucdot", enuc, x)

    if (not err == 0):
        sys.exit(2)

    
    pylab.plot(x, enuc, color="blue")
    pylab.ylabel(r"nuclear energy generation rate (erg/g/s)")

    pylab.xlabel(r"x (cm)")

    ax = pylab.gca()
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))



    #--------------------------------------------------------------------------
    # save the figure
    #--------------------------------------------------------------------------
    f = pylab.gcf()
    f.set_size_inches(6.0,12.0)

    if (not eps):
        pylab.savefig(outFile, bbox_inches='tight', dpi=dpi, pad_inches=0.33)
    else:
        pylab.savefig(outFile, bbox_inches='tight', pad_inches=0.33)



#==============================================================================
# usage
#==============================================================================
def usage():
    usageStr = """
    ./plot1dflame.py [options] plotfile 

    Make a plot of select data from a 1-d flame plotfile.

    Options:

       -o outfile   save the plot to the file outfile

       --eps        make an EPS plot instead of a PNG

       --dpi value  (PNG only) make the plot with the dpi specified by
                    value

    Note: this script requires the fsnapshot.so library, compiled with
    f2py using the GNUmakefile in data_processing/python_plotfile/

    """              
    print usageStr



#==============================================================================
# main
#==============================================================================
if __name__== "__main__":

    outFile = ""
    eps = 0
    dpi = 100

    try: opts, next = getopt.getopt(sys.argv[1:], "o:", 
                                    ["eps","dpi="])

    except getopt.GetoptError:
        print "invalid calling sequence"
        usage()
        sys.exit(2) 
               

    for o, a in opts:

        if o == "-o":
            outFile = a

        if o == "--eps":
            eps = 1

        if o == "--dpi":
            try: dpi = int(a)
            except ValueError:
                print "invalid value for --dpi"
                sys.exit(2)



    try: plotfile = next[0]
    except IndexError:
        print "ERROR: plotfile not specified"
        usage()
        sys.exit(2)

    do_plot(plotfile, outFile, eps, dpi)

