#!/usr/bin/env python
import sys
import os
import commands

fconv_slopesPath="/data_processing/MAESTRO_xrb/"
networksPath="/extern/networks/"
EOSPath="/extern/EOS/helmeos/"

def compileit(network,fParallel):
    """ 
This routine compiles the fconv_slopes.f90 routine with the specified network,
and then copies the executable to the current working directory.  This assumes 
that this script is in data_processing/python/ in order to find the location 
of the fconv_slopes.f90 routine; if it can't be found, we complain."""

    cwd = os.getcwd()

    fParallelA = os.path.abspath(os.path.expanduser(fParallel))

    # make sure the network is valid
    if os.path.exists(fParallelA + networksPath + network):
        networkPath = networksPath + network
    else:
        print "Network %s not found in %s" % (network,fParallelA+networksPath)
        sys.exit()

    print "Building fconv_slopes..."
    os.chdir(fParallelA + fconv_slopesPath)
    cmd = "make programs=fconv_slopes NETWORK=%s NDEBUG=" % networkPath
    print cmd
    (status, out) = commands.getstatusoutput(cmd)
    if status != 0: 
        print "Compilation Error:"
        print out
        sys.exit()
    cmd = "cp fconv_slopes.*.exe %s" % cwd
    print cmd
    (status, out) = commands.getstatusoutput(cmd)
    if status != 0:
        print "Couldn't copy fconv_slopes*.exe to %s" % cwd
        print out
        sys.exit()
    cmd = "ln -s %s/%s/helm_table.dat" % (fParallelA,EOSPath)
    if not os.path.exists("./helm_table.dat"):
        print cmd
        (status, out) = commands.getstatusoutput(cmd)
    if status != 0:
        print "Couldn't create softlink to %s%shelm_table.dat" % \
            (fParallelA,EOSPath)
        print out
        sys.exit()
    os.chdir(cwd)

    

def runit(network, fParallel, input, output,*args):
    
    compileit(network,fParallel)

    print "Running fconv_slopes..."
    cmd = "./fconv_slopes*debug*.exe -i %s" % input
    if output: cmd += " -o %s" % output
    (status,out) = commands.getstatusoutput(cmd)
    if status != 0:
        print out
    

if __name__ == "__main__":
    import optparse

    usage="""
This script compiles the data_processing/fconv_slopes.f90 routine with a 
specified (via the -n or --network options) reaction network, then executes 
the fconv_slopes routine on the input file (-i or --input) and dumps the result
to the output file (-o or --output; defaults to <inputFile>.out)."""
    
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-f", "--fParallel", dest="fParallel",
                      default="~/MAESTRO/fParallel/",
                      help="specify the location of the fParallel directory [default: %default]")
    parser.add_option("-n", "--network", dest="network", 
                      help="specify which NETWORK to compile with")
    parser.add_option("-i", "--input", dest="inputFile",
                      help="specify the input FILE", metavar="FILE")
    parser.add_option("-o", "--output", dest="outputFile",
                      help="specify the output FILE", metavar="FILE")

    (options, args) = parser.parse_args()

    if not options.network:
        print 'You must specify a network to use for compilation of '
        print 'fconv_slopes.f90'
        print "example: %s -n hotcno" % sys.argv[0]
        sys.exit()
    if not options.inputFile:
        print 'You must specify an input file for fconv_slopes.f90.'
        print 'example: %s -i model.hse' % sys.argv[0]
        sys.exit()

    if not os.path.exists(os.path.expanduser(options.fParallel)):
        print "Can't find the fParallel directory."
        print "Currently pointed to %s" % options.fParallel
        print "Change this with the -f option."
        sys.exit()

    runit(options.network, 
          options.fParallel,
          options.inputFile, options.outputFile)
