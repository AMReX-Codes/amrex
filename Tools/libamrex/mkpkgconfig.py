#!/usr/bin/env python3

import argparse
import re

def pkg_version(git_version):
    """Extract a clean MAJOR.MINOR[.PATCH] from git describe output.

    git describe produces strings like 25.06, 25.06.1, 25.06-3-g1a2b3c4, or
    25.06.1-3-g1a2b3c4-dirty.  We take the leading numeric tag as the
    pkg-config version so that version comparisons work.
    """
    m = re.match(r'^(\d+\.\d+(?:\.\d+)?)', git_version)
    return m.group(1) if m else git_version

def doit(prefix, version, cflags, libs, libpriv, fflags):
    pkg_ver = pkg_version(version)
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
    print("Description: Software Framework for Block Structured AMR")
    print("Version: "+pkg_ver)
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
