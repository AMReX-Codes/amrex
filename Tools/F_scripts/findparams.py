#!/usr/bin/env python

# a simple script that returns (to stdout) the paths to _parameter
# files given a list of directories to search (as arguments on the
# commandline).  This is used in various makefiles to get the input
# for write_probin.py

from __future__ import print_function

import sys
import os

def findparams(param_file_dirs):

    params = []

    for d in param_file_dirs:
        f = os.path.normpath(d + "/_parameters")
        if (os.path.isfile(f)):
            params.append(f)

    for f in params:
        print(f,)


if __name__ == "__main__":

    param_file_dirs = sys.argv[1:]

    findparams(param_file_dirs)
