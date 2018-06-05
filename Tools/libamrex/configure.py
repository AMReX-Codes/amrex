#!/usr/bin/env python

from __future__ import print_function

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for configure.py")

import argparse

def configure(argv):
    argv[0] = "configure" # So the help message print it
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",
                        help="Install libarmex, headers and Fortran modules in PREFIX directory [default=tmp_install_dir]",
                        default="tmp_install_dir")
    parser.add_argument("--dim", 
                        help="Dimension [default=3]",
                        choices=['1','2','3'],
                        default="3")
    parser.add_argument("--with-mpi",
                        help="Use MPI [default=yes]",
                        choices=["yes","no"],
                        default="yes")
    parser.add_argument("--with-omp",
                        help="Use OpenMP [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--comp",
                        help="Compiler [default=gnu]",
                        choices=["gnu","intel","cray","pgi","llvm","nag"],
                        default="gnu")
    parser.add_argument("--debug",
                        help="Debug build [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--enable-particle",
                        help="Enable AMReX particle classes [default=yes]",
                        choices=["yes","no"],
                        default="yes")
    parser.add_argument("--enable-fortran-api",
                        help="Enable AMReX Fortran API [default=yes]",
                        choices=["yes","no"],
                        default="yes")
    parser.add_argument("--enable-linear-solver",
                        help="Enable AMReX linear solvers [default=yes]",
                        choices=["yes","no"],
                        default="yes")
    parser.add_argument("--enable-xsdk-defaults",
                        help="Enable XSDK mode [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--allow-different-compiler",
                        help="Allow an application to use a different compiler than the one used to build libamrex [default=no]",
                        choices=["yes","no"],
                        default="no")
    args = parser.parse_args()

    f = open("GNUmakefile","w")
    f.write("AMREX_INSTALL_DIR = " + args.prefix.strip() + "\n")
    f.write("DIM = "+args.dim.strip() + "\n")
    if args.with_mpi == "no":
        f.write("USE_MPI = FALSE\n")
    else:
        f.write("USE_MPI = TRUE\n")
    if args.with_omp == "no":
        f.write("USE_OMP = FALSE\n")
    else:
        f.write("USE_OMP = TRUE\n")
    f.write("COMP = " + args.comp.strip() + "\n")
    if args.debug == "yes":
        f.write("DEBUG = TRUE\n")
    else:
        f.write("DEBUG = FALSE\n")
    if args.enable_particle == "no":
        f.write("USE_PARTICLES = FALSE\n")
    else:
        f.write("USE_PARTICLES = TRUE\n")
    if args.enable_fortran_api == "no":
        f.write("USE_FORTRAN_INTERFACE = FALSE\n")
    else:
        f.write("USE_FORTRAN_INTERFACE = TRUE\n")
    if args.enable_linear_solver == "no":
        f.write("USE_LINEAR_SOLVERS = FALSE\n")
    else:
        f.write("USE_LINEAR_SOLVERS = TRUE\n")
    if args.enable_xsdk_defaults == "yes":
        f.write("AMREX_XSDK = TRUE\n")
    else:
        f.write("AMREX_XSDK = FALSE\n")
    if args.allow_different_compiler == "no":
        f.write("ALLOW_DIFFERENT_COMP = FALSE\n")
    else:
        f.write("ALLOW_DIFFERENT_COMP = TRUE\n")

    f.write("\n")

    fin = open("GNUmakefile.in","r")
    for line in fin.readlines():
        f.write(line)
    fin.close()

    f.close()

if __name__ == "__main__":
    configure(sys.argv)
