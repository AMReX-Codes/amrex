#!/usr/bin/env python

from __future__ import print_function

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for mkconfig.py")

import argparse

def doit(defines, undefines, comp, allow_diff_comp, use_omp):
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

    print("#ifdef __cplusplus");

    if allow_diff_comp == "FALSE":
        if comp == "gnu" or comp == "nag":
            comp_macro = "__GNUC__"
            comp_id    = "GNU"
        elif comp == "intel":
            comp_macro = "__INTEL_COMPILER"
            comp_id    = Intel
        elif comp == "cray":
            comp_macro = "_CRAYC"
            comp_id    = Cray
        elif comp == "pgi":
            comp_macro = "__PGI"
            comp_id    = PGI
        elif comp == "llvm":
            comp_macro = "__llvm__"
            comp_id    = "Clang/LLVM"
        elif comp == "nec":
            comp_macro = "__NEC__"
            comp_id    = "NEC"
        elif comp == "ibm":
            comp_macro = "__ibmxl__"
            comp_id    = "IBM"
        else:
            sys.exit("ERROR: unknown compiler "+comp+" to mkconfig.py")

        msg = "libamrex was built with " + comp_id + ". "
        msg = msg + "To avoid this error, reconfigure with --allow-different-compiler=yes"
        print("#ifndef " + comp_macro )
        print('static_assert(false,"'+msg+'");')
        print("#endif")

    if use_omp == "TRUE":
        print("#ifndef _OPENMP")
        print('static_assert(false,"libamrex was built with OpenMP");')
        print("#endif")
    elif use_omp == "FALSE":
        print("#ifdef _OPENMP")
        print('static_assert(false,"libamrex was built without OpenMP");')
        print("#endif")
    else:
        sys.exit("ERROR: unknown use_omp flag "+use_omp+" in mkconfig.py")

    print("#endif") #  ifdef __cplusplus

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
                        choices=["gnu","intel","cray","pgi","llvm","nag","nec","ibm"])
    parser.add_argument("--allow-different-compiler",
                        help="allow an application to use a different compiler than the one used to build libamrex",
                        choices=["TRUE","FALSE"])
    parser.add_argument("--use-omp",
                        help="use openmp",
                        choices=["TRUE","FALSE"])
    args = parser.parse_args()

    try:
        doit(defines=args.defines, undefines=args.undefines, comp=args.comp,
             allow_diff_comp=args.allow_different_compiler,
             use_omp=args.use_omp)
    except:
        # something went wrong
        print("$(error something went wrong in mkconfig.py)")
