# (plotsinglevar.py):
# a simple script to plot 2-d or 3-d BoxLib data using the matplotlib
# library
#
# 2011-12-02 M. Zingale
#
# (plotsinglevar_parts.py):
# modify plotsinglevar.py routine to be a module for plotting 
# component data and particle data. This routine will be imported into 
# the particle plotting routine. Only intended to be used in 2-d, 
# but it should be easy to adapt it to include 3-d.
#
# 2012-3-02 R. Orvedahl
#

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
import parseparticles


#==============================================================================
# do_plot
#==============================================================================
def do_plot(plotfile, component, component2, outFile, log, 
            minval, maxval, minval2, maxval2, eps, dpi, origin,
            annotation, particles, time_ind):

    # plotfile            plotfile to read data from
    # component           first variable to plot
    # component2          optional second variable to plot
    # outFile             save plot as "outFile"
    # log (= 1,0)         plot log (base 10) of the data
    # minval (float)      specify minimum data range for plot 
    # maxval (float)      specify maximum data range for plot
    # eps (= 1,0)         make an EPS plot instead of PNG
    # dpi (int)           (PNG only) make the plot with dpi value
    # origin (= 1,0)      (3-d only) slice through origin (0,0,0)
    # annotation          (2-d only) add annotation string under time
    # particles           (2-d only) particle data
    # time_ind            array index of time cooresponding to plotfile

    #----------------------
    # check incoming values:
    #-----------------------
    if (plotfile == ""):
       print "\n---ERROR: plotfile not specified---"
       print   "          (plotsinglevar_parts.py)\n"
       usage()
       sys.exit(2)

    if (component == ""):
       print "\n---ERROR: no component specified---"
       print   "          (plotsinglevar_parts.py)\n"
       usage()
       sys.exit(2)

    if ((log != 1) and (log != 0)):
       print "\n---ERROR: invalid value for log (= 1,0)---"
       print   "          (plotsinglevar_parts.py)\n"
       usage()
       sys.exit(2)

    if ((eps != 1) and (eps != 0)):
       print "\n---ERROR: invalid value for eps (= 1,0)---"
       print   "          (plotsinglevar_parts.py)\n"
       usage()
       sys.exit(2)

    if (minval != None):
       try: minval = float(minval)
       except ValueError:
            print "\n---ERROR: invalid value for minval (= float)---"
            print   "          (plotsinglevar_parts.py)\n"
            usage()
            sys.exit(2)

    if (maxval != None):
       try: maxval = float(maxval)
       except ValueError:
            print "\n---ERROR: invalid value for maxval (= float)---"
            print   "          (plotsinglevar_parts.py)\n"
            usage()
            sys.exit(2)

    if (minval2 != None):
       try: minval2 = float(minval2)
       except ValueError:
            print "\n---ERROR: invalid value for minval2 (= float)---"
            print   "          (plotsinglevar_parts.py)\n"
            usage()
            sys.exit(2)

    if (maxval2 != None):
       try: maxval2 = float(maxval2)
       except ValueError:
            print "\n---ERROR: invalid value for maxval2 (= float)---"
            print   "          (plotsinglevar_parts.py)\n"
            usage()
            sys.exit(2)

    if ((origin != 1) and (origin != 0)):
       print "\n---ERROR: invalid value for origin (= 1,0)---"
       print   "          (plotsinglevar_parts.py)\n"
       usage()
       sys.exit(2)

    if (dpi != None):
       try: dpi = int(dpi)
       except ValueError:
            print "\n---ERROR: invalid value for dpi (= int)---"
            print   "          (plotsinglevar_parts.py)\n"
            usage()
            sys.exit(2)

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

    x = xmin + numpy.arange( (nx), dtype=numpy.float64 )*(xmax - xmin)/nx
    y = ymin + numpy.arange( (ny), dtype=numpy.float64 )*(ymax - ymin)/ny
    if (nz > 0):
        z = zmin + numpy.arange( (nz), dtype=numpy.float64 )*(zmax - zmin)/nz


    if (nz == -1):

        #----------------------------------------------------------------------
        # 2-d plots
        #----------------------------------------------------------------------
        extent = xmin, xmax, ymin, ymax

        # read in the main component
        data = numpy.zeros( (nx, ny), dtype=numpy.float64)

        (data, err) = fsnapshot.fplotfile_get_data_2d(plotfile, component, data)
        if (not err == 0):
            sys.exit(2)

        data = numpy.transpose(data)
        if (minval == None): minval = numpy.min(data)
        if (maxval == None): maxval = numpy.max(data)

        # read in the component #2, if present
        if (not component2 == ""):
            data2 = numpy.zeros( (nx, ny), dtype=numpy.float64)

            (data2, err) = fsnapshot.fplotfile_get_data_2d(plotfile, component2, data2)
            if (not err == 0):
                sys.exit(2)

            data2 = numpy.transpose(data2)
            if (minval2 == None): minval2 = numpy.min(data2)
            if (maxval2 == None): maxval2 = numpy.max(data2)

        if log:
            data = numpy.log10(data)

            if (not component2 == ""):
                data2 = numpy.log10(data2)
                minval2 = math.log10(minval2)
                maxval2 = math.log10(maxval2)
                
            minval = math.log10(minval)
            maxval = math.log10(maxval)


        #----------------------------------------------------------------------
        # plot main component
        #----------------------------------------------------------------------
        if (not component2 == ""):
            ax = pylab.subplot(1,2,1)
            pylab.subplots_adjust(wspace=0.5)

        else:
            ax = pylab.subplot(1,1,1)
    

        divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)    

        im=pylab.imshow(data,origin='lower', extent=extent, vmin=minval, vmax=maxval)

        #----------------------------------------------------------------------
        # plot the particle data
        #----------------------------------------------------------------------
        n=0
        while(n < len(particles)):

            # sometimes the length of particle history is larger than index
            if (time_ind < len(particles[n].history)):
                pylab.scatter(particles[n].history[time_ind].xyz[0], 
                           particles[n].history[time_ind].xyz[1], 
                           s=0.5,c='black',marker='o') 
            n += 1

        pylab.title(component)
        pylab.xlabel("x")
        pylab.ylabel("y")

        # axis labels in scientific notation with LaTeX
        ax = pylab.gca()
        ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

        # make space for a colorbar -- getting it the same size as the
        # vertical extent of the plot is surprisingly tricky.  See
        # http://matplotlib.sourceforge.net/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
        ax_cb = divider.new_horizontal(size="10%", pad=0.1)

        fig1 = ax.get_figure()
        fig1.add_axes(ax_cb)

        formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
        cb = pylab.colorbar(im, format=formatter, cax=ax_cb)

        # make the font size for the colorbar axis labels small.  Note


        #----------------------------------------------------------------------
        # plot component #2
        #----------------------------------------------------------------------
        if (not component2 == ""):
            ax = pylab.subplot(1,2,2)

            divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)

            im = pylab.imshow(data2, origin='lower', extent=extent, 
                              vmin=minval2, vmax=maxval2)

            #----------------------------------------------------------------
            # plot the particle data
            #----------------------------------------------------------------
            n=0
            while(n < len(particles)):

                # sometimes the length of particle history is larger than index
                if (time_ind < len(particles[n].history)):
                   pylab.scatter(particles[n].history[time_ind].xyz[0], 
                              particles[n].history[time_ind].xyz[1], 
                              s=0.5,c='black',marker='o')
                n += 1

            pylab.title(component2)
            pylab.xlabel("x")
            pylab.ylabel("y")

            # axis labels in scientific notation with LaTeX
            ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
            ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

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

        ###################################
        #  NOT SUPPORTED YET
        ###################################
        print "\n\n--- ERROR: 3-d not yet implemented ---"
        print     "          (plotsinglevar_parts.py)\n"
        sys.exit(2)

        # figure out the maximum width -- we will plot xy, xz, and yz
        w1 = xmax - xmin  # also w2
        w3 = ymax - ymin

        if (w3 > w1):
            scale = w3
        else:
            scale = w1
        

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


        # x-y
        extent = xmin, xmax, ymin, ymax

        # read in the main component
        data = numpy.zeros( (nx, ny), dtype=numpy.float64)

        indir = 3
        (data, err) = \
            fsnapshot.fplotfile_get_data_3d(plotfile, component, indir, origin, data)
        if (not err == 0):
            sys.exit(2)

        data = numpy.transpose(data)

        if log:
            data = numpy.log10(data)
                
            if (not minval == None): minval = math.log10(minval)
            if (not maxval == None): maxval = math.log10(maxval)


        ax = pylab.subplot(1,3,1)
        pylab.subplots_adjust(wspace=0.4)
        #fig.add_axes(pos1)

        im=pylab.imshow(data,origin='lower', extent=extent, 
                        vmin=minval, vmax=maxval, axes=pos1)

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

        pylab.text(xmin, ymin, "time = %7.3g s" % (time), 
                   verticalalignment="bottom", transform = trans, 
                   clip_on=False, fontsize=10)


        # x-z
        extent = xmin, xmax, zmin, zmax

        # read in the main component
        data = numpy.zeros( (nx, nz), dtype=numpy.float64)
        (data, err) = \
            fsnapshot.fplotfile_get_data_3d(plotfile, component, 2, origin, data)
        if (not err == 0):
            sys.exit(2)

        data = numpy.transpose(data)

        if log:
            data = numpy.log10(data)
                
            if (not minval == None): minval = math.log10(minval)
            if (not maxval == None): maxval = math.log10(maxval)


        ax = pylab.subplot(1,3,2)
        #fig.add_axes(pos2)

        im=pylab.imshow(data,origin='lower', extent=extent, 
                        vmin=minval, vmax=maxval, axes=pos2)

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
        extent = ymin, ymax, zmin, zmax

        # read in the main component
        data = numpy.zeros( (ny, nz), dtype=numpy.float64)
        (data, err) = \
            fsnapshot.fplotfile_get_data_3d(plotfile, component, 1, origin, data)
        if (not err == 0):
            sys.exit(2)

        data = numpy.transpose(data)

        if log:
            data = numpy.log10(data)
                
            if (not minval == None): minval = math.log10(minval)
            if (not maxval == None): maxval = math.log10(maxval)


        ax = pylab.subplot(1,3,3)
        #fig.add_axes(pos3)

        im=pylab.imshow(data,origin='lower', extent=extent, 
                        vmin=minval, vmax=maxval, axes=pos3)

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

    if (not eps):
        pylab.savefig(outFile, bbox_inches='tight', dpi=dpi, pad_inches=0.33)
    else:
        pylab.savefig(outFile, bbox_inches='tight', pad_inches=0.33)



#==============================================================================
# usage
#==============================================================================
def usage():
    usageStr = """
    do_plot(plotfile, component, component2, outFile, log, 
            minval, maxval, minval2, maxval2, eps, dpi, origin,
            annotation, particles, time_ind):

    Make a simple colormap plot of variable "component" from the
    BoxLib plotfile "plotfile". Support for 3-d not yet implemented.

    If component2 is a non-empty string then two side-by-side plots 
    will be made, one for each specified component.

    Variables:

       outFile     save the plot to the file outFile

       minval      set the minimum data range for the plot of component

       maxval      set the maximum data range for the plot of component

       minval2     set the minimum data range for the plot of component2

       maxval2     set the maximum data range for the plot of component2

       log         plot the logarithm (base-10) of the data

       eps         make an EPS plot instead of a PNG

       dpi         (PNG only) make the plot with the dpi specified by
                     value

       origin      (3-d only) slice through the origin (0,0,0) instead
                    of the center of the domain.

       annotation  (2-d only) add text "annotation" under the time

       particles   particle data

       time_ind    array index of the time cooresponding to the time in 
                    "plotfile"

    Note: this script requires the fsnapshot.so library, compiled with
    f2py using the GNUmakefile in data_processing/python_plotfile/

    """              
    print usageStr

