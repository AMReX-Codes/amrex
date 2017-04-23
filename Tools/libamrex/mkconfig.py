#!/usr/bin/env python

from __future__ import print_function

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for mkconfig.py")

import argparse

def doit(defines, undefines, comp, allow_diff_comp):
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

    if comp == "gnu":
        print("#ifndef __GNUC__")
        print('static_assert(false,"libamrex was compiled with GNU");')
        print("#endif")
    elif comp == "intel":
        print("#ifndef __INTEL_COMPILER")
        print('static_assert(false,"libamrex was compiled with Intel");')
        print("#endif")
    elif comp == "cray":
        print("#ifndef _CRAYC")
        print('static_assert(false,"libamrex was compiled with Cray");')
        print("#endif")
    elif comp == "pgi":
        print("#ifndef __PGI")
        print('static_assert(false,"libamrex was compiled with PGI");')
        print("#endif")
    elif comp == "llvm":
        print("#ifndef __llvm__")
        print('static_assert(false,"libamrex was compiled with Clang/LLVM");')
        print("#endif")
    else:
        sys.exit("ERROR: unknown compiler "+comp+" to mkconfig.py")

    print("#endif")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--defines",
                        help="preprocessing macros: -Dxx -Dyy",
                        default="")
    parser.add_argument("--undefines",
                        help="preprocessing macros to be undefined",
                        default="")
    parser.add_argument("--comp",
                        help="compiler",
                        choices=["gnu","intel","cray","pgi","llvm"])
    parser.add_argument("--allow-different-compiler",
                        help="allow an application to use a different compiler than the one used to build libamrex",
                        choices=["TRUE","FALSE"])
    args = parser.parse_args()

    try:
        doit(defines=args.defines, undefines=args.undefines, comp=args.comp,
             allow_diff_comp=args.allow_different_compiler)
    except:
        # something went wrong
        print("$(error something went wrong in mkconfig.py)")

