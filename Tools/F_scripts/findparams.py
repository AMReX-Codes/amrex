#!/usr/bin/env python

import sys
import os

def findparams(paramFileDirs):

    params = []

    for d in paramFileDirs:
        f = d + "/_parameters"
        if (os.path.isfile(f)):
            params.append(f)

    for f in params:
        print f,


if __name__ == "__main__":

    paramFileDirs = sys.argv[1:]

    findparams(paramFileDirs)
