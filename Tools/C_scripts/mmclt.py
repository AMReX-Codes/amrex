#!/usr/bin/env python3

"""
A script processing ccache log file to make Make file for Clang-Tidy

This generates a makefile for clang-tidying ccache-cache-missing files. This
could be used to speed up clang-tidy in CIs.
"""

import argparse, re, sys

def mmclt(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",
                        help="Ccache log file",
                        default="ccache.log.txt")
    parser.add_argument("--identifier",
                        help="Unique identifier for finding compilation line in the log file",
                        default="Src/Base")
    # We assume Src/Base can be used as an identifier to distinguish amrex code
    # from cmake's temporary files like build/CMakeFiles/CMakeScratch/TryCompile-hw3x4m/test_mpi.cpp
    parser.add_argument("--output",
                        help="Make file for clang-tidy",
                        default="clang-tidy-ccache-misses.mak")
    args = parser.parse_args()

    fin = open(args.input, "r")
    fout = open(args.output, "w")

    fout.write("CLANG_TIDY ?= clang-tidy\n")
    fout.write("override CLANG_TIDY_ARGS += --extra-arg=-Wno-unknown-warning-option --extra-arg-before=--driver-mode=g++\n")
    fout.write("\n")

    fout.write(".SECONDEXPANSION:\n")
    fout.write("clang-tidy: $$(all_targets)\n")
    fout.write("\t@echo SUCCESS\n\n")

    exe_re = re.compile(r" Executing .*? (-.*{}.*) -c .* -o .* (\S*)".format(args.identifier))

    count = 0
    for line in fin.readlines():
        ret_exe_re = exe_re.search(line)
        if (ret_exe_re):
            fout.write("target_{}: {}\n".format(count, ret_exe_re.group(2)))
            fout.write("\t$(CLANG_TIDY) $(CLANG_TIDY_ARGS) $< -- {}\n".format
                       (ret_exe_re.group(1)))
            fout.write("\ttouch target_{}\n\n".format(count))
            count = count + 1

    fout.write("all_targets =")
    for i in range(count):
        fout.write(" target_{}".format(i))
    fout.write("\n\n")

    fout.write("clean:\n\t$(RM) $(all_targets)\n\n")

    fout.close()
    fin.close()

if __name__ == "__main__":
    mmclt(sys.argv)
