#!/usr/bin/env python

from __future__ import print_function

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for mkpkgconfig.py")

import argparse

def doit(prefix, version, cflags, libs, libpriv, fflags):
    print("# AMReX Version: "+version)
    print("")
    print("prefix="+prefix)
    print("exec_prefix=${prefix}")
    print("libdir=${prefix}/lib")
    print("includedir=${prefix}/include")
    print("")
    print("fflags="+fflags);
    print("")
    print("Name: amrex")
    print("Description: Software Framework for Block Structed AMR")
    print("Version:")
    print("URL: https://github.com/AMReX-Codes/amrex")
    print("Requires:")
    print("Cflags: -I${includedir}", cflags)
    print("Libs: -L${libdir} -lamrex", libs)
    print("Libs.private:", libpriv)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",
                        help="prefix",
                        default="")
    parser.add_argument("--version",
                        help="version",
                        default="")
    parser.add_argument("--cflags",
                        help="cflags",
                        default="")
    parser.add_argument("--libs",
                        help="libs",
                        default="")
    parser.add_argument("--libpriv",
                        help="libpriv",
                        default="")
    parser.add_argument("--fflags",
                        help="fflags",
                        default="")
    args = parser.parse_args()

    try:
        doit(prefix=args.prefix, version=args.version, cflags=args.cflags,
             libs=args.libs, libpriv=args.libpriv, fflags=args.fflags)
    except:
        # something went wrong
        print("$(error something went wrong in mkpkgconfig.py)")
