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
                        help="Install libamrex, headers and Fortran modules in PREFIX directory [default=tmp_install_dir]",
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
    parser.add_argument("--with-cuda",
                        help="Use CUDA [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--with-acc",
                        help="Use OpenACC [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--comp",
                        help="Compiler [default=gnu]",
                        choices=["gnu","intel","cray","pgi","llvm","nag","nec","ibm"],
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
    parser.add_argument("--enable-hypre",
                        help="Enable Hypre as an option for bottom solver of AMReX linear solvers [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--enable-petsc",
                        help="Enable PETSc as an option for bottom solver of AMReX linear solvers [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--enable-eb",
                        help="Enable AMReX embedded boundary capability [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--enable-xsdk-defaults",
                        help="Enable XSDK mode [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--allow-different-compiler",
                        help="Allow an application to use a different compiler than the one used to build libamrex [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--with-sensei-insitu",
                        help="Use SENSEI in situ [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--with-omp-offload",
                        help="Use OpenMP-offload [default=no]",
                        choices=["yes","no"],
                        default="no")
    parser.add_argument("--enable-tiny-profile",
                        help="Enable tiny profile [default=no]",
                        choices=["yes","no"],
                        default="no")
    args = parser.parse_args()

    f = open("GNUmakefile","w")
    f.write("AMREX_INSTALL_DIR = " + args.prefix.strip() + "\n")
    f.write("DIM = " + args.dim.strip() + "\n")
    f.write("USE_MPI = {}\n".format("FALSE" if args.with_mpi == "no" else "TRUE"))
    f.write("USE_OMP = {}\n".format("FALSE" if args.with_omp == "no" else "TRUE"))
    f.write("USE_CUDA = {}\n".format("FALSE" if args.with_cuda == "no" else "TRUE"))
    f.write("USE_ACC = {}\n".format("FALSE" if args.with_acc == "no" else "TRUE"))
    f.write("COMP = " + args.comp.strip() + "\n")
    f.write("DEBUG = {}\n".format("TRUE" if args.debug == "yes" else "FALSE"))
    f.write("USE_PARTICLES = {}\n".format("FALSE" if args.enable_particle == "no" else "TRUE"))
    f.write("USE_FORTRAN_INTERFACE = {}\n".format("FALSE" if args.enable_fortran_api == "no" else "TRUE"))
    f.write("USE_LINEAR_SOLVERS = {}\n".format("FALSE" if args.enable_linear_solver == "no" else "TRUE"))
    f.write("USE_HYPRE = {}\n".format("TRUE" if args.enable_hypre == "yes" else "FALSE"))
    f.write("USE_PETSC = {}\n".format("TRUE" if args.enable_petsc == "yes" else "FALSE"))
    f.write("USE_EB = {}\n".format("TRUE" if args.enable_eb == "yes" else "FALSE"))
    f.write("AMREX_XSDK = {}\n".format("TRUE" if args.enable_xsdk_defaults == "yes" else "FALSE"))
    f.write("ALLOW_DIFFERENT_COMP = {}\n".format("FALSE" if args.allow_different_compiler == "no" else "TRUE"))
    f.write("USE_SENSEI_INSITU = {}\n".format("FALSE" if args.with_sensei_insitu == "no" else "TRUE"))
    f.write("USE_OMP_OFFLOAD = {}\n".format("FALSE" if args.with_omp_offload == "no" else "TRUE"))
    f.write("TINY_PROFILE = {}\n".format("FALSE" if args.enable_tiny_profile == "no" else "TRUE"))
    f.write("\n")

    fin = open("GNUmakefile.in","r")
    for line in fin.readlines():
        f.write(line)
    fin.close()

    f.close()

if __name__ == "__main__":
    configure(sys.argv)
