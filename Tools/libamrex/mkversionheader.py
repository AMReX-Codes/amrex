#!/usr/bin/env python3

import sys, re
import argparse

def doit(code, defines):
    print("#ifndef "+code+"_VERSION_H_")
    print("#define "+code+"_VERSION_H_")

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

    print("#endif // "+code+"_VERSION_H_")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--code",
                        help="code name",
                        default="AMREX")
    parser.add_argument("--defines",
                        help="preprocessing macros: -Dxx -Dyy",
                        default="")
    args = parser.parse_args()

    try:
        doit(code=args.code, defines=args.defines)
    except:
        # something went wrong
        print("$(error something went wrong in mkversionheader.py)")
