#!/usr/bin/env python
import math
import numpy
import pylab
import os
import argparse


def getFileParams(filename):

    # get the number of lines and fields in the data file
    # skip comment and blank lines, and also, require that
    # the time (first field) is monotonically increasing

    mf = open(filename, "r")

    numLines = 0
    numFields = -1

    oldTime = -1.0

    for line in mf:
        if (not (line.startswith("#") or line.lstrip() == "")):

            fields = line.split()
            time = float(fields[0])

            if (numFields == -1):
                numFields = len(fields)
                oldTime = time
                numLines += 1

            elif (time > oldTime):
                oldTime = time
                numLines +=1
            
    mf.close()

    return numLines, numFields


def getData(filename):

    # get the data from the data file
    # skip comment and blank lines, and also, require that
    # the time (first field) is monotonically increasing

    numLines, numFields = getFileParams(filename)

    data = numpy.zeros( (numLines, numFields), numpy.float64)

    mf = open(filename, "r")

    indexLine = 0

    oldTime = -1.0

    for line in mf:
        if (not (line.startswith("#") or line.lstrip() == "")):

            fields = line.split()
            time = float(fields[0])

            if (oldTime == -1.0):
                data[indexLine,:] = fields
                oldTime = time
                indexLine += 1

            elif (time > oldTime):
                data[indexLine,:] = fields
                oldTime = time
                indexLine +=1
            
    mf.close()

    return data


def do_diagnostics():

    # resizing
    #pylab.rcParams.update({'xtick.labelsize': 20, 
    #                       'ytick.labelsize': 20,
    #                       'text.fontsize': 20})

    #pylab.rc("axes", linewidth=2.0)
    #pylab.rc("lines", markeredgewidth=2.0)

    specFile = os.path.join(args.infile, "all_spectrum.dat")

    specData = getData(specFile)

    Ek = specData[:,1] + specData[:,3] + specData[:,5]

    # do a crude "fit"
    #n = specData.shape[0] - 1
    #k = numpy.arange(1000) + 10
    #Efit = Ek[n]*(specData[n,0]/k)**(5./3.)

    #pylab.plot(k, Efit,
    #           color="r", ls="--", label=r"$k^{-5/3}$")   #, lw=2)


    # plot the simulation data
    pylab.plot(specData[:,0], Ek,
               color="k", label=r"spectrum")   #, lw=2)


    # axis stuff
    ax = pylab.gca()
    ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
    
    ax.set_xscale('log')
    ax.set_yscale('log')

    pylab.xlim(1, 5000)
    #pylab.ylim(1.e9,1.e12)

    pylab.xlabel(r"$k$")   #, fontsize=28)

    pylab.ylabel(r"$E(k)$")   #, fontsize=28)

    leg = pylab.legend(loc=1)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize=12)
    leg.draw_frame(0)

    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.savefig(os.path.join(args.infile, "spectrum.png"))


if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=str, help='Name of input plotfile to plot spectra after processing with AmrDeriveSpectrum.')
    args = parser.parse_args()

    do_diagnostics()

