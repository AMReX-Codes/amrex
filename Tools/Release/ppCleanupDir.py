#!/usr/bin/env python

import os
import shutil
import sys
import getopt
import readline
import subprocess
import ppCleanup as ppc
    
def ppCleanupDir(argv):

    usage = """
    /path/to/ppCleanupDir.py CLEANWORD_FILE InDir OutDir
    InDir must exist and OutDir must not
    """
    if len(argv) != 4:
        print usage
        sys.exit(2)

    cwfile = argv[1]
    indir  = argv[2]
    outdir = argv[3]

    if not os.path.isfile(cwfile):
        print "clean word file:",cwfile,"does't exist."
        sys.exit(1)

    if not os.path.isdir(indir):
        print "In Dir", indir,"does't exist."
        sys.exit(1)

    if os.path.exists(outdir):
        print "Out Dir", outdir,"already exists."
        sys.exit(1)

    shutil.copytree(indir, outdir)

    print "\nCleaning up IFDEFs in", indir, "and store new files in", outdir
    for root, dirs, files in os.walk(outdir):
        for name in files:
            f = os.path.join(root,name)
            fext = os.path.splitext(f)[1]
            if fext in ['.f','.f90','.F','.F90','.c','.cpp','.h','.H']:
                ftmp = f+'.orig'
                os.rename(f,ftmp)
                ppc.ppCleanup(cwfile, ftmp, f)
                os.remove(ftmp)

    print "Finish cleaning up"


if __name__== "__main__":
    ppCleanupDir(sys.argv)
