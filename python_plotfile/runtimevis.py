import numpy
import pylab
import matplotlib
import os
import sys
import getopt
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


def doPlot(ax, name, data):
    im = ax.imshow(data, origin="lower")
    ax.set_title(name)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.cax.colorbar(im)


def setupAxes(aspectRatio, nvar):

    if (aspectRatio == "h"):

        # for <= 3 variables, do a single column
        # for 4 <= # var <= 6, do two columns

        if (nvar <= 3):
            grid = ImageGrid(F, 111, # similar to subplot(111)
                             nrows_ncols = (nvar, 1),
                             direction="row",
                             axes_pad = 0.5 ,
                             add_all=True,
                             label_mode = "L",
                             share_all = True,
                             cbar_location="top",
                             cbar_mode="each",
                             cbar_size="5%",
                             cbar_pad="15%")
            
        elif (nvar == 4):
            grid = ImageGrid(F, 111, # similar to subplot(111)
                             nrows_ncols = (2, 2),
                             direction="row",
                             axes_pad = 0.5 ,
                             add_all=True,
                             label_mode = "L",
                             share_all = True,
                             cbar_location="top",
                             cbar_mode="each",
                             cbar_size="5%",
                             cbar_pad="15%")

        else:
            grid = ImageGrid(F, 111, # similar to subplot(111)
                             nrows_ncols = (3, 2),
                             direction="row",
                             axes_pad = 0.5 ,
                             add_all=True,
                             label_mode = "L",
                             share_all = True,
                             cbar_location="top",
                             cbar_mode="each",
                             cbar_size="5%",
                             cbar_pad="20%")

    elif (aspectRatio == "v"):
        
        # always do 1 row
        grid = ImageGrid(F, 111, # similar to subplot(111)
                         nrows_ncols = (1, nvar),
                         direction="row",
                         axes_pad = 0.2 ,
                         add_all=True,
                         label_mode = "L",
                         share_all = True,
                         cbar_location="top",
                         cbar_mode="each",
                         cbar_size="3%",
                         cbar_pad="8%")

    else:
        
        # for <= 3 variables, do a single row
        # for 4 <= # var <= 6, do 2 rows. 
        if (nvar <= 3):
            grid = ImageGrid(F, 111, # similar to subplot(111)
                             nrows_ncols = (1, nvar),
                             direction="row",
                             axes_pad = 0.2 ,
                             add_all=True,
                             label_mode = "L",
                             share_all = True,
                             cbar_location="top",
                             cbar_mode="each",
                             cbar_size="5%",
                             cbar_pad="15%")

        elif (nvar == 4):
            grid = ImageGrid(F, 111, # similar to subplot(111)
                             nrows_ncols = (2, 2),
                             direction="row",
                             axes_pad = 0.5 ,
                             add_all=True,
                             label_mode = "L",
                             share_all = True,
                             cbar_location="top",
                             cbar_mode="each",
                             cbar_size="5%",
                             cbar_pad="15%")
        else:
            grid = ImageGrid(F, 111, # similar to subplot(111)
                             nrows_ncols = (2, 3),
                             direction="row",
                             axes_pad = 0.5 ,
                             add_all=True,
                             label_mode = "L",
                             share_all = True,
                             cbar_location="top",
                             cbar_mode="each",
                             cbar_size="5%",
                             cbar_pad="15%")



        return grid



# main

# get the grid info
nx = 10
ny = 10


# get the data

nvar = 3

# find the aspect ratio:
#
# aspectRatio = "h" means horizontal
#               "v" means vertical
#               "s" means square (to some degree...)

if (nx >= 2*ny):
    aspectRatio = "h"
elif (ny >= 2*nx):
    aspectRatio = "v"
else:
    aspectRatio = "s"



# setup the figure
F = pylab.figure(1, (12.8, 7.2)) 
F.clf()


# setup the axes
grid = setupAxes(aspectRatio, nvar)


# plot the data
dummy = numpy.arange(nx*ny)
dummy.shape = ny, nx

n = 0
while (n < nvar):
    doPlot(grid[n], "var %s" % (n), dummy)
    n += 1

# 5 variables is a tricky case
if (nvar == 5 and (aspectRatio == "h" or aspectRatio == "s")):
    # turn off the last axes
    grid[5].axis('off')
    grid[5].cax.axis('off')


pylab.savefig("test.png")



