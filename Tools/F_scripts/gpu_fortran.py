#!/usr/bin/env python

"""
Search the Fortran source code for subroutines marked with !$gpu in their specification part:

  subroutine sub(a)

    integer :: a

    !$gpu

    ...

  end subroutine

For each one, prepend AMREX_DEVICE prior to the subroutine statement:

  AMREX_DEVICE subroutine sub(a)

    integer :: a

    !$gpu

    ...

  end subroutine
"""

from __future__ import print_function

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for gpu_fortran.py")

if sys.version[0] == "2":
    reload(sys)
    sys.setdefaultencoding('utf8')

import os
import argparse



def update_fortran_subroutines(ffile):
    """For subroutines marked up with !$gpu, replace with
    AMREX_DEVICE."""

    # open the input Fortran file
    try:
        fin = open(ffile, "r")
    except IOError:
        sys.exit("Cannot open Fortran file {} for reading".format(ffile))

    # read in the file, close it, and then reopen it for writing
    lines = fin.readlines()

    fin.close()

    try:
        fout = open(ffile, "w")
    except IOError:
        sys.exit("Cannot open Fortran file {} for writing".format(ffile))

    # loop through the file and look for all subroutines

    found_subroutine = False
    subroutine = ""

    for line in lines:

        if "subroutine" in line and not found_subroutine:
            found_subroutine = True

        if "end subroutine" in line:

            # Now if the subroutine contains !$gpu, prepend
            # AMREX_DEVICE to the subroutine, if it does not
            # already contain it.

            if "!$gpu" in subroutine and "AMREX_DEVICE" not in subroutine:
                subroutine = subroutine.replace("subroutine", "AMREX_DEVICE subroutine", 1)

            fout.write(subroutine)
            subroutine = ""
            found_subroutine = False

        if found_subroutine:
            subroutine += line
        else:
            fout.write(line)

    fout.close()



if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--fortran",
                        help="the name of the Fortran file to update")

    args = parser.parse_args()

    update_fortran_subroutines(args.fortran)
