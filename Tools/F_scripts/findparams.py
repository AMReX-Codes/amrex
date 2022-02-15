#!/usr/bin/env python

from __future__ import print_function

import sys
import os

def findparams(param_file_dirs):
    """
    a simple script that returns (to stdout) the paths to _parameter
    files given a list of directories to search (as arguments on the
    commandline).  This is used in various makefiles to get the input for
    write_probin.py

    """

    params = []

    for d in param_file_dirs:
        f = os.path.normpath(d + "/_parameters")
        if os.path.isfile(f):
            params.append(f)

    for f in params:
        print(f,)


if __name__ == "__main__":

    if len(sys.argv) == 1:
        sys.exit("invalid calling sequence.\n  findparams.py path-to-dir1 path-to-dir2 ...\n")

    findparams(sys.argv[1:])
