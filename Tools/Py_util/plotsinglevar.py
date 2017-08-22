#!/usr/bin/env python

# a simple script to plot 2-d or 3-d BoxLib data using the matplotlib
# library
#

from __future__ import print_function

import matplotlib                                                               
matplotlib.use('agg')   

import numpy
import pylab
import os
import sys
import getopt
import math
import string
import mpl_toolkits.axes_grid1

# Give an informative error if we can't get fsnapsho
try:
    import fsnapshot
except ImportError as e:
    print("*** ERROR: {} ***".format(e))
    print("\nNote: this script requires the fsnapshot.so library, compiled", \
          "\nwith f2py using the GNUmakefile in AmrPostprocessing/python.")
    sys.exit()


#==============================================================================
# do_plot
#==============================================================================
def do_plot(plotfile, component, component2, outFile, log, 
            minval, maxval, minval2, maxval2, eps, dpi, origin, annotation, 
            xmin_pass, ymin_pass, zmin_pass, xmax_pass, ymax_pass, zmax_pass):

    pylab.rc("font", size=9)


    #--------------------------------------------------------------------------
    # construct the output file name
    #--------------------------------------------------------------------------
    if (outFile == ""):
        outFile = os.path.normpath(plotfile) + "_" + component
        if (not component2 == ""): outFile += "_" + component2

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

    dx = (xmax - xmin)/nx
    x = xmin + numpy.arange( (nx), dtype=numpy.float64 )*dx

    dy = (ymax - ymin)/ny
    y = ymin + numpy.arange( (ny), dtype=numpy.float64 )*dy

    if (nz > 0):
        dz = (zmax - zmin)/nz
        z = zmin + numpy.arange( (nz), dtype=numpy.float64 )*dz


    if (nz == -1):

        #----------------------------------------------------------------------
        # 2-d plots
        #----------------------------------------------------------------------
        extent = [xmin, xmax, ymin, ymax]

        ix0 = 0
        ix = nx
        iy0 = 0
        iy = ny

        if (not xmin_pass == None):
            extent[0] = xmin_pass
            ix0 = int((xmin_pass - xmin)/dx)

        if (not ymin_pass == None):
            extent[2] = ymin_pass
            iy0 = int((ymin_pass - ymin)/dy)

        if (not xmax_pass == None):
            extent[1] = xmax_pass
            ix = int((xmax_pass - xmin)/dx)

        if (not ymax_pass == None):
            extent[3] = ymax_pass
            iy = int((ymax_pass - ymin)/dy)


        sparseX = 0
        if (ny >= 3*nx):
            sparseX = 1

        # read in the main component
        data = numpy.zeros( (nx, ny), dtype=numpy.float64)

        (data, err) = fsnapshot.fplotfile_get_data_2d(plotfile, component, data)
        if (not err == 0):
            sys.exit(2)

        data = numpy.transpose(data)


        # read in the component #2, if present
        if (not component2 == ""):
            data2 = numpy.zeros( (nx, ny), dtype=numpy.float64)

            (data2, err) = fsnapshot.fplotfile_get_data_2d(plotfile, component2, data2)
            if (not err == 0):
                sys.exit(2)

            data2 = numpy.transpose(data2)

        # log?
        if log:
            if (numpy.min(data) < 0):
                data = numpy.log10(numpy.abs(data))
            else:
                data = numpy.log10(data)

            if (not component2 == ""): 
                if (numpy.min(data2) < 0.0):
                    data2 = numpy.log10(numpy.abs(data2))
                else:
                    data2 = numpy.log10(data2)
                
            if (not minval == None): minval = math.log10(minval)
            if (not maxval == None): maxval = math.log10(maxval)

            if (not minval2 == None): minval2 = math.log10(minval2)
            if (not maxval2 == None): maxval2 = math.log10(maxval2)


        #----------------------------------------------------------------------
        # plot main component
        #----------------------------------------------------------------------
        if (not component2 == ""):
            ax = pylab.subplot(1,2,1)
            pylab.subplots_adjust(wspace=0.5)

        else:
            ax = pylab.subplot(1,1,1)
    

        divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)    

        im=pylab.imshow(data[iy0:iy,ix0:ix],origin='lower', extent=extent, vmin=minval, vmax=maxval)

        pylab.title(component)
        pylab.xlabel("x")
        pylab.ylabel("y")

        # axis labels in scientific notation with LaTeX
        ax = pylab.gca()
        ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

        if (sparseX):
            ax.xaxis.set_major_locator(pylab.MaxNLocator(3))


        # make space for a colorbar -- getting it the same size as the
        # vertical extent of the plot is surprisingly tricky.  See
        # http://matplotlib.sourceforge.net/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
        ax_cb = divider.new_horizontal(size="10%", pad=0.1)

        fig1 = ax.get_figure()
        fig1.add_axes(ax_cb)

        formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
        cb = pylab.colorbar(im, format=formatter, cax=ax_cb)

        # make the font size for the colorbar axis labels small.  Note
        # the offsetText is the 10^N that appears at the top of the
        # y-axis.
        cl = pylab.getp(cb.ax, 'ymajorticklabels')
        pylab.setp(cl, fontsize=10)
        
        cb.ax.yaxis.offsetText.set_fontsize("small")


        # do a fixed offset in pixels from the (xmin,ymin) data point
        trans=matplotlib.transforms.offset_copy(ax.transData, x=0, y=-0.5, 
                                                fig=fig1, units='inches')

        pylab.text(xmin, ymin, "time = %7.3g s" % (time), 
                   verticalalignment="bottom", transform = trans, 
                   clip_on=False, fontsize=10)

        if (not annotation == ""):
            trans=matplotlib.transforms.offset_copy(ax.transData, x=0, y=-0.65,
                                                    fig=fig1, units='inches')

            pylab.text(xmin, ymin, "%s" % (annotation), 
                       verticalalignment="bottom", transform = trans, 
                       clip_on=False, fontsize=10)            


        #----------------------------------------------------------------------
        # plot component #2
        #----------------------------------------------------------------------
        if (not component2 == ""):
            ax = pylab.subplot(1,2,2)

            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)

            im = pylab.imshow(data2[iy0:iy,ix0:ix], origin='lower', extent=extent, 
                              vmin=minval2, vmax=maxval2)

            pylab.title(component2)
            pylab.xlabel("x")
            pylab.ylabel("y")

            # axis labels in scientific notation with LaTeX
            ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
            ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

            if (sparseX):
                ax.xaxis.set_major_locator(pylab.MaxNLocator(3))

            # make space for a colorbar -- getting it the same size as
            # the vertical extent of the plot is surprisingly tricky.
            # See
            # http://matplotlib.sourceforge.net/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
            ax_cb = divider.new_horizontal(size="10%", pad=0.1)

            fig1 = ax.get_figure()
            fig1.add_axes(ax_cb)

            formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
            cb = pylab.colorbar(im, format=formatter, cax=ax_cb)

            # make the font size for the colorbar axis labels small.  Note the
            # offsetText is the 10^N that appears at the top of the y-axis.
            cl = pylab.getp(cb.ax, 'ymajorticklabels')
            pylab.setp(cl, fontsize=10)
        
            cb.ax.yaxis.offsetText.set_fontsize("small")

            #ax_cb.yaxis.tick_right()


    else:

        #----------------------------------------------------------------------
        # 3-d plot
        #----------------------------------------------------------------------

        # starting points for the figure positions

        # assume that the width of the plotting area is 0.05 to 0.95,
        # leaving 0.9 to be split amongst the 3 plots.  So each has a
        # width of 0.3

        # for the height, we will assume that the colorbar at the
        # bottom gets 0.15, and that we go until 0.95, leaving 0.8 of
        # height for the plots.
            
        pos1 = [0.05, 0.15, 0.3, 0.8]
        pos2 = [0.35, 0.15, 0.3, 0.8]
        pos3 = [0.65, 0.15, 0.3, 0.8]

        fig = pylab.figure()


        # read in the slices
        # x-y
        data_xy = numpy.zeros( (nx, ny), dtype=numpy.float64)

        indir = 3
        (data_xy, err) = \
            fsnapshot.fplotfile_get_data_3d(plotfile, component, indir, 
                                            origin, data_xy)
        if (not err == 0):
            sys.exit(2)

        data_xy = numpy.transpose(data_xy)

        if log:
            if (numpy.min(data_xy) < 0):
                data_xy = numpy.log10(numpy.abs(data_xy))
            else:
                data_xy = numpy.log10(data_xy)


        # x-z
        data_xz = numpy.zeros( (nx, nz), dtype=numpy.float64)
        (data_xz, err) = \
            fsnapshot.fplotfile_get_data_3d(plotfile, component, 2, 
                                            origin, data_xz)
        if (not err == 0):
            sys.exit(2)

        data_xz = numpy.transpose(data_xz)

        if log:
            if (numpy.min(data_xz) < 0):
                data_xz = numpy.log10(numpy.abs(data_xz))
            else:
                data_xz = numpy.log10(data_xz)
                

        # y-z
        data_yz = numpy.zeros( (ny, nz), dtype=numpy.float64)
        (data_yz, err) = \
            fsnapshot.fplotfile_get_data_3d(plotfile, component, 1, 
                                            origin, data_yz)
        if (not err == 0):
            sys.exit(2)

        data_yz = numpy.transpose(data_yz)

        if log:
            if (numpy.min(data_yz) < 0):
                data_yz = numpy.log10(numpy.abs(data_yz))
            else:
                data_yz = numpy.log10(data_yz)
                


                
        if (not minval == None): 
            if (log):
                minval = math.log10(minval)
        else:
            minval = numpy.min(data_xy)
            minval = min(minval,numpy.min(data_xz))
            minval = min(minval,numpy.min(data_yz))


        if (not maxval == None): 
            if (log):
                maxval = math.log10(maxval)
        else:
            maxval = numpy.max(data_xy)
            maxval = max(maxval,numpy.max(data_xz))
            maxval = max(maxval,numpy.max(data_yz))
            



        # x-y
        extent = [xmin, xmax, ymin, ymax]

        ix0 = 0
        ix = nx
        iy0 = 0
        iy = ny

        if (not xmin_pass == None):
            extent[0] = xmin_pass
            ix0 = int((xmin_pass - xmin)/dx)

        if (not ymin_pass == None):
            extent[2] = ymin_pass
            iy0 = int((ymin_pass - ymin)/dy)

        if (not xmax_pass == None):
            extent[1] = xmax_pass
            ix = int((xmax_pass - xmin)/dx)

        if (not ymax_pass == None):
            extent[3] = ymax_pass
            iy = int((ymax_pass - ymin)/dy)



        ax = pylab.subplot(1,3,1)
        pylab.subplots_adjust(wspace=0.4)
        #fig.add_axes(pos1)

        im=pylab.imshow(data_xy[iy0:iy,ix0:ix],origin='lower', extent=extent, 
                        vmin=minval, vmax=maxval)#, axes=pos1)

        pylab.xlabel("x")
        pylab.ylabel("y")

        # axis labels in scientific notation with LaTeX
        ax = pylab.gca()
        ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

        ax.xaxis.offsetText.set_fontsize("small")
        ax.yaxis.offsetText.set_fontsize("small")

        cl = pylab.getp(ax, 'ymajorticklabels')
        pylab.setp(cl, fontsize=10)
        cl = pylab.getp(ax, 'xmajorticklabels')
        pylab.setp(cl, fontsize=10)

        # do a fixed offset in pixels from the (xmin,ymin) data point
        fig1 = ax.get_figure()
        trans=matplotlib.transforms.offset_copy(ax.transData, x=0, y=-0.5, 
                                                fig=fig1, units='inches')

        # pylab.text(xmin_pass, ymin_pass, "time = %7.3g s" % (time), 
        #            verticalalignment="bottom", transform = trans, 
        #            clip_on=False, fontsize=10)


        # x-z
        extent = [xmin, xmax, zmin, zmax]

        ix0 = 0
        ix = nx
        iz0 = 0
        iz = nz

        if (not xmin_pass == None):
            extent[0] = xmin_pass
            ix0 = int((xmin_pass - xmin)/dx)

        if (not zmin_pass == None):
            extent[2] = zmin_pass
            iz0 = int((zmin_pass - zmin)/dz)

        if (not xmax_pass == None):
            extent[1] = xmax_pass
            ix = int((xmax_pass - xmin)/dx)

        if (not zmax_pass == None):
            extent[3] = zmax_pass
            iz = int((zmax_pass - zmin)/dz)



        ax = pylab.subplot(1,3,2)
        #fig.add_axes(pos2)

        im=pylab.imshow(data_xz[iz0:iz,ix0:ix],origin='lower', extent=extent, 
                        vmin=minval, vmax=maxval) #, axes=pos2)

        pylab.xlabel("x")
        pylab.ylabel("z")

        # axis labels in scientific notation with LaTeX
        ax = pylab.gca()
        ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

        ax.xaxis.offsetText.set_fontsize("small")
        ax.yaxis.offsetText.set_fontsize("small")

        cl = pylab.getp(ax, 'ymajorticklabels')
        pylab.setp(cl, fontsize=10)
        cl = pylab.getp(ax, 'xmajorticklabels')
        pylab.setp(cl, fontsize=10)


        # y-z
        extent = [ymin, ymax, zmin, zmax]

        iy0 = 0
        iy = ny
        iz0 = 0
        iz = nz

        if (not ymin_pass == None):
            extent[0] = ymin_pass
            iy0 = int((ymin_pass - ymin)/dy)

        if (not zmin_pass == None):
            extent[2] = zmin_pass
            iz0 = int((zmin_pass - zmin)/dz)

        if (not ymax_pass == None):
            extent[1] = ymax_pass
            iy = int((ymax_pass - ymin)/dy)

        if (not zmax_pass == None):
            extent[3] = zmax_pass
            iz = int((zmax_pass - zmin)/dz)



        ax = pylab.subplot(1,3,3)
        #fig.add_axes(pos3)

        im=pylab.imshow(data_yz[iz0:iz,iy0:iy],origin='lower', extent=extent, 
                        vmin=minval, vmax=maxval) #, axes=pos3)

        pylab.xlabel("y")
        pylab.ylabel("z")

        # axis labels in scientific notation with LaTeX
        ax = pylab.gca()
        ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

        ax.xaxis.offsetText.set_fontsize("small")
        ax.yaxis.offsetText.set_fontsize("small")

        cl = pylab.getp(ax, 'ymajorticklabels')
        pylab.setp(cl, fontsize=10)
        cl = pylab.getp(ax, 'xmajorticklabels')
        pylab.setp(cl, fontsize=10)


        # colorbar
        pylab.subplots_adjust(bottom=0.1, left=0.05, right=0.95)

        formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
        
        cax = pylab.axes([0.05, 0.06, 0.9, 0.04])
        pylab.colorbar(orientation="horizontal", cax=cax, format=formatter)

        pylab.title(component)


    #--------------------------------------------------------------------------
    # save the figure
    #--------------------------------------------------------------------------
    # try: pylab.tight_layout()  # requires matplotlib >= 1.1
    # except:
    #     pass

    if (not eps):
        pylab.savefig(outFile, bbox_inches='tight', dpi=dpi, pad_inches=0.33)
    else:
        pylab.savefig(outFile)#, bbox_inches='tight', pad_inches=0.33)



#==============================================================================
# usage
#==============================================================================
def usage():
    usageStr = """
    ./plotsinglevar.py [options] plotfile component [component2]

    Make a simple colormap plot of variable "component" from the
    BoxLib plotfile "plotfile".  Support for 2-d and 3-d.

    If two component names are specified after the plotfile name then
    two side-by-side plots will be made, one for each specified
    component.

    Options:

       -o outfile   save the plot to the file outfile

       -m value     set the minimum data range for the plot to value

       -M value     set the maximum data range for the plot to value

       --log        plot the logarithm (base-10) of the data
                    (note: if the data < 0 anywhere, then the abs
                     is taken first)

       -x num       minimum x coordinate
       -X num       maximum x coordinate

       -y num       minimum y coordinate
       -Y num       maximum y coordinate

       -z num       minimum z coordinate
       -Z num       maximum z coordinate

       --eps        make an EPS plot instead of a PNG

       --dpi value  (PNG only) make the plot with the dpi specified by
                    value

       --origin     (3-d only) slice through the origin (0,0,0) instead
                    of the center of the domain.

       --annotate string (2-d only) add an annotation under the time


    Note: this script requires the fsnapshot.so library, compiled with
    f2py using the GNUmakefile in AmrPostprocessing/python.

    """              
    print(usageStr)



#==============================================================================
# main
#==============================================================================
if __name__== "__main__":

    outFile = ""
    log = 0
    eps = 0
    minvar = None
    maxvar = None
    minvar2 = None
    maxvar2 = None
    dpi = 100
    origin = 0
    annotation = ""
    xmin = None
    ymin = None
    zmin = None
    xmax = None
    ymax = None
    zmax = None

    try: opts, next = getopt.getopt(sys.argv[1:], "o:m:M:n:N:x:y:z:X:Y:Z:", 
                                    ["log","eps","dpi=","origin","annotate="])
    except getopt.GetoptError:
        print("invalid calling sequence")
        usage()
        sys.exit(2) 
               

    for o, a in opts:

        if o == "-o":
            outFile = a

        if o == "-m":
            try: minvar = float(a)
            except ValueError:
                print("invalid value for -m")
                sys.exit(2)

        if o == "-M":
            try: maxvar = float(a)
            except ValueError:
                print("invalid value for -M")
                sys.exit(2)

        if o == "-n":
            try: minvar2 = float(a)
            except ValueError:
                print("invalid value for -n")
                sys.exit(2)

        if o == "-N":
            try: maxvar2 = float(a)
            except ValueError:
                print("invalid value for -N")
                sys.exit(2)

        if o == "-x":
            try: xmin = float(a)
            except ValueError:
                print("invalid value for -x")
                sys.exit(2)            

        if o == "-y":
            try: ymin = float(a)
            except ValueError:
                print("invalid value for -y")
                sys.exit(2)            

        if o == "-z":
            try: zmin = float(a)
            except ValueError:
                print("invalid value for -z")
                sys.exit(2)            

        if o == "-X":
            try: xmax = float(a)
            except ValueError:
                print("invalid value for -X")
                sys.exit(2)            

        if o == "-Y":
            try: ymax = float(a)
            except ValueError:
                print("invalid value for -Y")
                sys.exit(2)            

        if o == "-Z":
            try: zmax = float(a)
            except ValueError:
                print("invalid value for -Z")
                sys.exit(2)            
 
        if o == "--log":
            log = 1

        if o == "--eps":
            eps = 1

        if o == "--dpi":
            try: dpi = int(a)
            except ValueError:
                print("invalid value for --dpi")
                sys.exit(2)

        if o == "--origin":
            origin = 1

        if o == "--annotate":
            annotation = a



    try: plotfile = next[0]
    except IndexError:
        print("ERROR: plotfile not specified")
        usage()
        sys.exit(2)

    try: component = next[1]
    except IndexError:
        print("ERROR: no component specified")
        usage()
        sys.exit(2)    

    try: component2 = next[2]
    except IndexError:
        component2 = ""

    do_plot(plotfile, component, component2, outFile, 
            log, minvar, maxvar, minvar2, maxvar2, eps, dpi, origin, annotation,
            xmin, ymin, zmin, xmax, ymax, zmax)
