#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import argparse

def find_files(args):

    files = []
    not_found = []

    if args.vpath is None:
        sys.exit("vpath is empty -- nothing to search")

    if args.files is None:
        sys.exit("files is empty -- nothing to search")

    filenames = args.files.split()
    vpath = args.vpath.split()

    for f in filenames:

        found = False

        for d in vpath:
            if os.path.isfile("{}/{}".format(d, f)):
                found = True
                files.append((f, d))
                break

        if not found:
            not_found.append(f)

    # output
    print("locations of source files:")
    for f, d in files:
        print("{} : {}".format(f, d))

    print("\n")
    print("files not found:")
    for f in not_found:
        print("{}".format(f))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--vpath", type=str, default=None,
                        metavar="'list of dirs in vpath'",
                        help="the list of directories that make will look at to find source files (e.g. the vpath)")
    parser.add_argument("--files", type=str, default=None,
                        metavar="'list of files'",
                        help="the list of source files")
    args = parser.parse_args()

    find_files(args)
