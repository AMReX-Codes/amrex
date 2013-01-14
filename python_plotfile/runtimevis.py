#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')   # this is important for batch mode on machines w/o a display
import numpy
import pylab
import os
import sys
import getopt
import ConfigParser
import fsnapshot

import math
import string
from mpl_toolkits.axes_grid1 import ImageGrid


# support doing runtime visualization in 2-d.
#
# a vis.in file is read in that is in the form
#
# [general]
# key = value
#
# [density]
# min = X
# max = Y
#
# ...
# 
# where each variable to be plotted gets its own block.  We then count
# the number of variables and plot them.
#
# We use the matplotlib Imagegrid to make the plot axes easy to setup.
#
# We look at the aspect ratio of the data to ensure that we use the
# best grid layout.
#
# Note: we always do a plot of 1280x720 pixels (or some multiple thereof).
# as this is 720p HD resolution (good for youtube).


#-----------------------------------------------------------------------------
class variable:

    def __init__(self, name="", minval=None, maxval=None, log=0):
        self.name = name
        self.min = minval
        self.max = maxval
        self.log = log
        self.data = None


    def __str__(self):
        if self.min == None:
            minStr = "None"
        else:
            minStr = `self.min`

        if self.max == None:
            maxStr = "None"
        else:
            maxStr = `self.max`

        str = "%s: range = [%s, %s], log = %d" % (self.name, minStr, maxStr, self.log)
        return str


class grid:

    def __init__ (self, xmin=0.0, ymin=0.0, xmax=1.0, ymax=1.0, 
                  dx=0.1, dy=0.1):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        self.dx = dx
        self.dy = dy


#-----------------------------------------------------------------------------
def parseInfile(inFile):

    vars = []

    parser=ConfigParser.SafeConfigParser()
    parser.read(inFile)

    if (parser.sections() == []):
        sys.exit("ERROR: no variables defined")

    for section in parser.sections():
        vars.append(variable(section))
        
        for option in parser.options(section):

            if option == "min":
                try: value=parser.getfloat(section,option)
                except ValueError:
                    sys.exit("invalid min for %s" % (section))

                vars[len(vars)-1].min = value

            elif option == "max":
                try: value=parser.getfloat(section,option)
                except ValueError:
                    sys.exit("invalid max for %s" % (section))

                vars[len(vars)-1].max = value

            elif option == "log":
                try: value=parser.getint(section,option)
                except ValueError:
                    sys.exit("invalid log for %s" % (section))

                vars[len(vars)-1].log = value

            else:
                sys.exit("invalid option for %s" % (section))


        #print vars[len(vars)-1]   # debugging
    return vars

    
#-----------------------------------------------------------------------------
def setupAxes(F, aspectRatio, nvar):

    if (aspectRatio == "h"):

        # for <= 3 variables, do a single column
        # for 4 <= # var <= 6, do two columns

        if (nvar <= 3):
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (nvar, 1), direction="row",
                               axes_pad = 0.5 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="15%")
            
        elif (nvar == 4):
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (2, 2), direction="row",
                               axes_pad = 0.5 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="15%")

        else:
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (3, 2), direction="row",
                               axes_pad = 0.5 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="20%")

    elif (aspectRatio == "v"):

        # always do 1 row -- just much with the spacings here

        if (nvar <= 4):

            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (1, nvar), direction="row",
                               axes_pad = 0.2 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="3%", cbar_pad="8%")

        else:

            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (1, nvar), direction="row",
                               axes_pad = 0.2 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="15%")

    else:
        
        # for <= 3 variables, do a single row
        # for 4 <= # var <= 6, do 2 rows. 
        if (nvar <= 3):
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (1, nvar), direction="row",
                               axes_pad = 0.2 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="10%")
            
        elif (nvar == 4):
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (2, 2), direction="row",
                               axes_pad = 0.5 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="15%")
        else:
            axGrid = ImageGrid(F, 111, # similar to subplot(111)
                               nrows_ncols = (2, 3), direction="row",
                               axes_pad = 0.5 ,
                               add_all=True,
                               label_mode = "L",
                               share_all = True,
                               cbar_location="top", cbar_mode="each",
                               cbar_size="5%", cbar_pad="15%")



    return axGrid


#-----------------------------------------------------------------------------
def doPlot(ax, grd, var):
    extent = [grd.xmin, grd.xmax, grd.ymin, grd.ymax]

    if var.log:
        pData = numpy.log10(var.data)
        if (not var.min == None): 
            pmin = math.log10(var.min)
        else:
            pmin = None
        if (not var.max == None):
            pmax = math.log10(var.max)
        else:
            pmax = None
    else:
        pData = var.data
        pmin = var.min
        pmax = var.max

    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3,3))

    im = ax.imshow(pData, origin="lower", interpolation="nearest",
                   vmin=pmin, vmax=pmax, extent=extent)

    ax.set_title(var.name)

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

    ax.cax.colorbar(im, format=formatter)


#-----------------------------------------------------------------------------
def main(inFile, plotFile):

    # get a list of variable objects that contains the information
    # about what to plot
    vars = parseInfile(inFile)

    nvar = len(vars)


    # get and store the grid info
    (nx, ny, nz) = fsnapshot.fplotfile_get_size(plotFile)
    if (not nz == -1):
        sys.exit("ERROR: cannot read a 3-d dataset")


    (xmin, xmax, ymin, ymax, zmin, zmax) = \
        fsnapshot.fplotfile_get_limits(plotFile)

    dx = (xmax - xmin)/nx
    x = xmin + numpy.arange( (nx), dtype=numpy.float64 )*dx

    dy = (ymax - ymin)/ny
    y = ymin + numpy.arange( (ny), dtype=numpy.float64 )*dy

    gridInfo = grid(xmin=xmin, xmax=xmax, 
                    ymin=ymin, ymax=ymax, 
                    dx=dx, dy=dy)


    time = fsnapshot.fplotfile_get_time(plotFile)


    # get the data
    for v in vars:
        data = numpy.zeros( (nx, ny), dtype=numpy.float64 )
        (data, err) = fsnapshot.fplotfile_get_data_2d(plotFile, v.name, data)
        if (not err == 0):
            sys.exit("ERROR: unable to read %s" % (v.name) )

        v.data = numpy.transpose(data)


    # find the aspect ratio:
    #
    # aspectRatio = "h" means horizontal
    #               "v" means vertical
    #               "s" means square (to some degree...)

    if (nx >= 2*ny):
        aspectRatio = "h"
    elif (ny >= 1.5*nx):
        aspectRatio = "v"
    else:
        aspectRatio = "s"



    # setup the figure
    F = pylab.figure(1, (12.8, 7.2)) 
    F.clf()

    pylab.rc("font", size=9)


    # setup the axes
    axGrid = setupAxes(F, aspectRatio, nvar)


    # plot the data
    n = 0
    while (n < nvar):
        doPlot(axGrid[n], gridInfo, vars[n])
        n += 1


    # 5 variables is a tricky case
    if (nvar == 5 and (aspectRatio == "h" or aspectRatio == "s")):
        # turn off the last axes
        axGrid[5].axis('off')
        axGrid[5].cax.axis('off')


    # write the time   
    print "writing time"
    F.text(0.1, 0.1, "t = %g s" % (time), transform = F.transFigure, color="k")



    pylab.savefig("%s.png" % (plotFile) )



if __name__ == "__main__":

    # parse the commandline options
    inFile = "vis.in"

    try: opts, next = getopt.getopt(sys.argv[1:], "i:")
    except getopt.GetoptError:
        sys.exit("ERROR: invalid calling sequence")

    for o, a in opts:
        if o == "-i":
            inFile = a


    try: plotFile = next[0]
    except IndexError:
        sys.exit("ERROR: plotfile not specified")

    main(inFile, plotFile)

