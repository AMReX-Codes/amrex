#!/usr/bin/env python

# a simple script that returns (to stdout) the paths to _parameter
# files given a list of directories to search (as arguments on the
# commandline).  This is used in various makefiles to get the input
# for write_probin.py

import sys
import os

def findparams(paramFileDirs):

    params = []

    for d in paramFileDirs:
        f = os.path.normpath(d + "/_parameters")
        if (os.path.isfile(f)):
            params.append(f)

    for f in params:
        print f,


if __name__ == "__main__":

    paramFileDirs = sys.argv[1:]

    findparams(paramFileDirs)
