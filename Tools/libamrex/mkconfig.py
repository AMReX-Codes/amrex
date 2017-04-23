#!/usr/bin/env python

from __future__ import print_function

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for mkconfig.py")

import argparse

def doit(defines, undefines):
    print("#ifndef AMREX_CONFIG_H_")
    print("#define AMREX_CONFIG_H_")
    defs = defines.split("-D")
    for d in defs:
        dd = d.strip()
        if dd:
            v = dd.split("=")
            if len(v) == 2:
                print("#define",v[0],v[1])
            else:
                print("#define",v[0],1)
    print("#undef",undefines)
    print("#endif")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--defines",
                        help="preprocessing macros: -Dxx -Dyy",
                        default="")
    parser.add_argument("--undefines",
                        help="preprocessing macros to be undefined",
                        default="")

    args = parser.parse_args()

    try:
        doit(defines=args.defines, undefines=args.undefines)
    except:
        # something went wrong
        print("$(error something went wrong in mkconfig.py)")

