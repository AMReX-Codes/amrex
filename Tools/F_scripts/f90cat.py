#!/usr/bin/env python3

# concatenate a bunch of Fortran files, sorting them so that the modules are
# defined before they are used

from __future__ import print_function

import sys
import toposort

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for dep.py")

if sys.version[0] == "2":
    reload(sys)
    sys.setdefaultencoding('utf8')


import re
import os
import argparse

# modules to ignore in the dependencies
IGNORES = ["iso_c_binding", "iso_fortran_env", "omp_lib", "mpi"]

# regular expression for "{}module{}name", where {} can be any number
# of spaces.  We use 4 groups here, denoted by (), so the name of the
# module is the 4th group
module_re = re.compile(r"^(\s*)([Mm][Oo][Dd][Uu][Ll][Ee])(\s+)((?:[a-z][a-z_0-9]+))",
                       re.IGNORECASE|re.DOTALL)

# regular expression for "{}module{}procedure{}name"
module_proc_re = re.compile(r"(\s*)(module)(\s+)(procedure)(\s+)((?:[a-z][a-z_0-9]+))",
                            re.IGNORECASE|re.DOTALL)

# regular expression for "{}use{}modulename...".  Note this will work for
# use modulename, only: stuff, other stuff'
# see (txt2re.com)
use_re = re.compile(r"^(\s*)([Uu][Ss][Ee])(\s+)((?:[a-z_][a-z_0-9]+))",
                    re.IGNORECASE|re.DOTALL)


def is_sorted(files):
    """ this should check whether the files are properly sorted"""
    have_sorted = True

    for n, currfile in enumerate(files):
        # nothing before the current file should depend on it
        for q in range(n):
            prevfile = files[q]
            for d in prevfile.needs:
                if d in currfile.defined:
                    have_sorted = False
                    break

    return have_sorted


class SourceFile(object):
    """ hold information about one of the .f90/.F90 files """

    def __init__(self, filename):

        self.name = filename
        self.ext = os.path.splitext(filename)[1]

        self.defined = self.defined_modules()
        self.needs = self.needed_modules()

    def __str__(self):
        return self.name

    def defined_modules(self):
        """determine what modules this file provides -- we work off of the
        preprocessed file if it exists."""

        defines = []

        with open(self.name, "r") as f:

            for line in f:

                # strip off the comments
                idx = line.find("!")
                line = line[:idx]

                # we want a module definition itself, not a 'module procedure'
                # also, Fortran is case-insensitive
                rebreak = module_re.search(line)
                rebreak2 = module_proc_re.search(line)
                if rebreak and not rebreak2:
                    defines.append(rebreak.group(4).lower())

        return defines

    def needed_modules(self):
        """determine what modules this file needs.  Assume only one use
           statement per line.  Ignore any only clauses.  We work off
           of the preprocessed file if it exists."""

        depends = []

        with open(self.name, "r") as f:

            for line in f:

                # strip off the comments
                idx = line.find("!")
                line = line[:idx]

                rebreak = use_re.search(line)
                if rebreak:
                    used_module = rebreak.group(4)
                    if used_module.lower() not in IGNORES:
                        depends.append(used_module.lower())

        # remove duplicates
        depends = list(set(depends))

        return depends


def doit(files):
    """ main routine that processes the files"""


    # all_files is a list of SourceFile objects
    all_files = []

    for cf in files:
        all_files.append(SourceFile(cf))


    # create a dictionary where the keys are the file and the values are a
    # set of the files it depends on
    deps = {}
    for f in all_files:
        dep_files = []
        for need in f.needs:
            dep_files += [q.name for q in all_files if need in q.defined]
        deps[f.name] = set(dep_files)

    sorted_names = toposort.toposort_flatten(deps)

    # now make it back into a list
    sorted_files = []
    for name in sorted_names:
        sorted_files += [q for q in all_files if q.name == name]
    
    print("sort check: {}".format(is_sorted(sorted_files)))

    # now concatenate
    with open("mega_f.F90", "w") as mf:
        for f in sorted_files:
            with open(f.name, "r") as sf:
                lines = sf.readlines()
                mf.writelines(lines)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("files", metavar="source files", type=str, nargs="*",
                        help="F90 source files to find dependencies amongst")

    args = parser.parse_args()

    #try:
    doit(args.files)
    #except:
    #    # something went wrong
    #    print("$(error something went wrong in f90cat.py.)")

