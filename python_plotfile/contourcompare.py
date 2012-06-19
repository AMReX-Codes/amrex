#!/usr/bin/env python

# a simple script to plot contours of a single variable from 2 (or 3) different 
# BoxLib plotfiles on the same axes for comparison.  This uses the 
# matplotlib library.
#
# 2012-04-12 M. Zingale

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


#==============================================================================
# do_plot
#==============================================================================
def do_plot(plotfile1, plotfile2, plotfile3, component, outFile, 
            log, minval, maxval, ncontours, eps, dpi, 
            xmin, xmax, ymin, ymax,
            label1, label2, label3):


    #--------------------------------------------------------------------------
    # construct the output file name
    #--------------------------------------------------------------------------
    if (outFile == ""):
        outFile = "compare_" + component

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
    # read in the data from plotfile1
    #--------------------------------------------------------------------------
    (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile1)

    time = fsnapshot.fplotfile_get_time(plotfile1)

    (xmin1, xmax1, ymin1, ymax1, zmin1, zmax1) = \
        fsnapshot.fplotfile_get_limits(plotfile1)

    x1 = xmin1 + numpy.arange( (nx), dtype=numpy.float64 )*(xmax1 - xmin1)/nx
    y1 = ymin1 + numpy.arange( (ny), dtype=numpy.float64 )*(ymax1 - ymin1)/ny


    if (not nz == -1):
        print "ERROR: 2-d support only"
        sys.exit(2)


    # read in the main component
    data1 = numpy.zeros( (nx, ny), dtype=numpy.float64)

    (data1, err) = fsnapshot.fplotfile_get_data_2d(plotfile1, component, data1)
    if (not err == 0):
        sys.exit(2)

    data1 = numpy.transpose(data1)


    if log:
        data1 = numpy.log10(data1)


    extent1 = [xmin1, xmax1, ymin1, ymax1]


    #--------------------------------------------------------------------------
    # read in the data from plotfile2
    #--------------------------------------------------------------------------
    (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile2)

    time = fsnapshot.fplotfile_get_time(plotfile2)

    (xmin2, xmax2, ymin2, ymax2, zmin2, zmax2) = \
        fsnapshot.fplotfile_get_limits(plotfile2)

    x2 = xmin2 + numpy.arange( (nx), dtype=numpy.float64 )*(xmax2 - xmin2)/nx
    y2 = ymin2 + numpy.arange( (ny), dtype=numpy.float64 )*(ymax2 - ymin2)/ny


    if (not nz == -1):
        print "ERROR: 2-d support only"
        sys.exit(2)


    # read in the main component
    data2 = numpy.zeros( (nx, ny), dtype=numpy.float64)

    (data2, err) = fsnapshot.fplotfile_get_data_2d(plotfile2, component, data2)
    if (not err == 0):
        sys.exit(2)

    data2 = numpy.transpose(data2)


    if log:
        data2 = numpy.log10(data2)


    extent2 = [xmin2, xmax2, ymin2, ymax2]


    #--------------------------------------------------------------------------
    # read in the data from plotfile3 -- if present
    #--------------------------------------------------------------------------
    if (not plotfile3 == ""):
        (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotfile3)

        time = fsnapshot.fplotfile_get_time(plotfile3)

        (xmin3, xmax3, ymin3, ymax3, zmin3, zmax3) = \
            fsnapshot.fplotfile_get_limits(plotfile3)

        x3 = xmin3 + numpy.arange( (nx), dtype=numpy.float64 )*(xmax3 - xmin3)/nx
        y3 = ymin3 + numpy.arange( (ny), dtype=numpy.float64 )*(ymax3 - ymin3)/ny


        if (not nz == -1):
            print "ERROR: 2-d support only"
            sys.exit(2)


        # read in the main component
        data3 = numpy.zeros( (nx, ny), dtype=numpy.float64)

        (data3, err) = fsnapshot.fplotfile_get_data_2d(plotfile3, component, data3)
        if (not err == 0):
            sys.exit(2)

        data3 = numpy.transpose(data3)


        if log:
            data3 = numpy.log10(data3)


        extent3 = [xmin3, xmax3, ymin3, ymax3]



    #--------------------------------------------------------------------------
    # find data limits, etc.
    #--------------------------------------------------------------------------
    if (not extent1 == extent2):
        print "ERROR: extent of domains do not agree"
        sys.exit(2)

    if (not plotfile3 == ""):
        if (not extent1 == extent3):
            print "ERROR: extent of domains do not agree"
            sys.exit(2)


    extent = extent1

    if (not xmin == None):
        extent[0] = xmin

    if (not xmax == None):
        extent[1] = xmax

    if (not ymin == None):
        extent[2] = ymin

    if (not ymax == None):
        extent[3] = ymax


    if (minval == None):
        minval = min(numpy.min(data1), numpy.min(data2))
        if (not plotfile3 == ""):
            minval = min(minval, numpy.min(data3))

    if (maxval == None):
        maxval = max(numpy.max(data1), numpy.max(data2))
        if (not plotfile3 == ""):
            maxval = max(maxval, numpy.max(data3))


    levels = numpy.arange(ncontours)*(maxval-minval)/(ncontours-1)


    #--------------------------------------------------------------------------
    # make the figure
    #--------------------------------------------------------------------------
    cs1 = pylab.contour(x1, y1, data1, ncontours, colors='k', levels=levels)
    cs2 = pylab.contour(x2, y2, data2, ncontours, colors='r', levels=levels)
    if (not plotfile3 == ""):
        cs3 = pylab.contour(x3, y3, data3, ncontours, colors='g', levels=levels)

    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    #pylab.clabel(cs, fontsize=9, inline=1)#, fmt=formatter)

    pylab.axis(extent)

    ax = pylab.gca()
    ax.set_aspect("equal")

    fig1 = ax.get_figure()

    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    pylab.xlabel("x")
    pylab.ylabel("y")

    if (not label1 == None):
        trans=matplotlib.transforms.offset_copy(ax.transData, x=0, y=-0.5,
                                                fig=fig1, units='inches')
        
        pylab.text(xmin, ymin, "%s" % (label1), 
                   verticalalignment="bottom", transform = trans, 
                   clip_on=False, fontsize=10, color="k")            


    if (not label2 == None):
        trans=matplotlib.transforms.offset_copy(ax.transData, x=0, y=-0.65,
                                                fig=fig1, units='inches')

        pylab.text(xmin, ymin, "%s" % (label2), 
                   verticalalignment="bottom", transform = trans, 
                   clip_on=False, fontsize=10, color="r")            


    if (not label3 == None):
        trans=matplotlib.transforms.offset_copy(ax.transData, x=0, y=-0.8,
                                                fig=fig1, units='inches')

        pylab.text(xmin, ymin, "%s" % (label3), 
                   verticalalignment="bottom", transform = trans, 
                   clip_on=False, fontsize=10, color="g")            


    if (not eps):
        pylab.savefig(outFile, bbox_inches='tight', dpi=dpi, pad_inches=0.5)
    else:
        pylab.savefig(outFile, bbox_inches='tight', pad_inches=0.5)



#==============================================================================
# usage
#==============================================================================
def usage():
    usageStr = """
    ./contourcompare.py [options] component plotfile1 plotfile2 [plotfile3] 

    Plot contours of variable "component" from both plotfile1 and
    plotfile2 (and plotfile3, if present) on the same axes for
    comparison.

    Options:

       -o outfile    save the plot to the file outfile

       -m value      set the minimum data range for the plot to value
       -M value      set the maximum data range for the plot to value

       -n value      set the number of contours to use

       -x value      x minimum for plot
       -X value      x maximum for plot

       -y value      y minimum for plot
       -Y value      y maximum for plot

       --log         plot the logarithm (base-10) of the data

       --eps         make an EPS plot instead of a PNG

       --dpi value   (PNG only) make the plot with the dpi specified by
                     value

       --label1 str  label for file1
       --label2 str  label for file2
       --label3 str  label for file3

    Note: this script requires the fsnapshot.so library, compiled with
    f2py using the GNUmakefile in data_processing/python_plotfile/

    """              
    print usageStr



#==============================================================================
# main
#==============================================================================
if __name__== "__main__":

    outFile = ""
    log = 0
    eps = 0

    minvar = None
    maxvar = None

    dpi = 100

    ncontours = 10

    xmin = None
    xmax = None
    ymin = None
    ymax = None

    label1 = None
    label2 = None
    label3 = None

    try: opts, next = getopt.getopt(sys.argv[1:], "o:m:M:n:x:X:y:Y:", 
                                    ["log","eps","dpi=",
                                     "label1=","label2=", "label3="])
    except getopt.GetoptError:
        print "invalid calling sequence"
        usage()
        sys.exit(2) 
               

    for o, a in opts:

        if o == "-o":
            outFile = a

        if o == "-m":
            try: minvar = float(a)
            except ValueError:
                print "invalid value for -m"
                sys.exit(2)

        if o == "-M":
            try: maxvar = float(a)
            except ValueError:
                print "invalid value for -M"
                sys.exit(2)

        if o == "-n":
            try: ncontours = int(a)
            except ValueError:
                print "invalid value for -n"
                sys.exit(2)

        if o == "-x":
            try: xmin = float(a)
            except ValueError:
                print "invalid value for -x"
                sys.exit(2)

        if o == "-X":
            try: xmax = float(a)
            except ValueError:
                print "invalid value for -X"
                sys.exit(2)

        if o == "-y":
            try: ymin = float(a)
            except ValueError:
                print "invalid value for -y"
                sys.exit(2)

        if o == "-Y":
            try: ymax = float(a)
            except ValueError:
                print "invalid value for -Y"
                sys.exit(2)

        if o == "--log":
            log = 1

        if o == "--eps":
            eps = 1

        if o == "--dpi":
            try: dpi = int(a)
            except ValueError:
                print "invalid value for --dpi"
                sys.exit(2)

        if o == "--label1":
            label1 = a

        if o == "--label2":
            label2 = a

        if o == "--label3":
            label3 = a



    try: component = next[0]
    except IndexError:
        print "ERROR: no component specified"
        usage()
        sys.exit(2)    

    try: plotfile1 = next[1]
    except IndexError:
        print "ERROR: plotfile1 not specified"
        usage()
        sys.exit(2)

    try: plotfile2 = next[2]
    except IndexError:
        print "ERROR: plotfile2 not specified"
        usage()
        sys.exit(2)

    try: plotfile3 = next[3]
    except IndexError:
        plotfile3 = ""



    if (not minvar == None) and log == 1:
        minvar = math.log10(minvar)

    if (not maxvar == None) and log == 1:
        maxvar = math.log10(maxvar)

    do_plot(plotfile1, plotfile2, plotfile3, component, outFile, 
            log, minvar, maxvar, ncontours, eps, dpi, 
            xmin, xmax, ymin, ymax,
            label1, label2, label3)
