#!/usr/bin/env python3

import sys, re
import argparse

def doit(defines, undefines, comp, allow_diff_comp):
    print("#ifndef AMREX_HAVE_NO_CONFIG_H")
    print("#define AMREX_HAVE_NO_CONFIG_H")

    # Remove -I from input
    defines = re.sub(r'-I.*?(-D|$)', r'\1', defines)

    defs = defines.split("-D")
    for d in defs:
        dd = d.strip()
        if dd:
            v = dd.split("=")
            print("#ifndef",v[0])
            if len(v) == 2:
                print("#define",v[0],v[1])
            else:
                print("#define",v[0],1)
            print("#endif")

    for ud in undefines:
        print("#undef",ud)

    print("#ifdef __cplusplus");

    if allow_diff_comp == "FALSE":
        if comp == "gnu" or comp == "nag":
            comp_macro = "__GNUC__"
            comp_id    = "GNU"
        elif comp == "intel" or comp == "intel-classic":
            comp_macro = "__INTEL_COMPILER"
            comp_id    = "Intel"
        elif comp == "intel-llvm":
            comp_macro = "__INTEL_LLVM_COMPILER"
            comp_id    = "Intel"
        elif comp == "cray":
            comp_macro = "_CRAYC"
            comp_id    = "Cray"
        elif comp == "pgi":
            comp_macro = "__PGI"
            comp_id    = "PGI"
        elif comp == "nvhpc":
            comp_macro = "__NVCOMPILER"
            comp_id    = "NVHPC"
        elif comp == "llvm":
            comp_macro = "__llvm__"
            comp_id    = "Clang/LLVM"
        elif comp == "nec":
            comp_macro = "__NEC__"
            comp_id    = "NEC"
        elif comp == "ibm":
            comp_macro = "__ibmxl__"
            comp_id    = "IBM"
        elif comp == "armclang":
            comp_macro = "__armclang_version__"
            comp_id    = "ArmClang"
        elif comp == "hip":
            comp_macro = "__HIP__"
            comp_id    = "HIP"
        elif comp == "sycl":
            comp_macro = "__INTEL_CLANG_COMPILER"
            comp_id    = "SYCL"
        else:
            sys.exit("ERROR: unknown compiler "+comp+" to mkconfig.py")

        msg = "#error libamrex was built with " + comp_id + ". "
        msg = msg + "To avoid this error, reconfigure with --allow-different-compiler=yes"
        print("#ifndef " + comp_macro )
        print(msg)
        print("#endif")

    print("#endif") #  ifdef __cplusplus

    # hipcc does not necessarily set the _OPENMP macro
    # https://rocmdocs.amd.com/en/latest/Programming_Guides/HIP-FAQ.html?highlight=_openmp#openmp-is-undefined-when-compiling-with-fopenmp
    print("#if defined(AMREX_USE_OMP) && !defined(_OPENMP) && !defined(AMREX_USE_HIP)")
    print('#error libamrex was built with OpenMP')
    print("#endif")

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
                        choices=["gnu","intel","intel-llvm","intel-classic","cray","pgi","nvhpc","llvm","nag","nec","ibm",
                                 "armclang","hip","sycl"])
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
