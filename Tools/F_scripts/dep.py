#!/usr/bin/env python

# automatically generate Makefile dependencies for Fortran 90 source.
#
# this will output all the dependency pairs amongst the source files.


import sys
import re
import os
import argparse


IGNORES = ["iso_c_binding", "iso_fortran_env"]


def doit(prefix, search_path, files):

    # regular expression for ' use modulename, only: stuff, other stuff'
    # see (txt2re.com)
    use_re = re.compile("( )(use)(\s+)((?:[a-z_][a-z_0-9]+))", 
                        re.IGNORECASE|re.DOTALL)

    module_re = re.compile("( )(module)(\s+)((?:[a-z][a-z_0-9]+))",
                           re.IGNORECASE|re.DOTALL)

    module_proc_re = re.compile("( )(module)(\s+)(procedure)(\s+)((?:[a-z][a-z_0-9]+))",
                                re.IGNORECASE|re.DOTALL)

    # first parse the files and find all the module statements.  Keep a
    # dictionary of 'module name':filename.
    modulefiles = {}

    all_files = []

    for cf in files:
        
        # find the file in the first part of the search path it exists
        if len(search_path) > 0:
            for p in search_path:
                full_file = "{}/{}".format(p, cf)
                if os.path.isfile(full_file): 
                    break
        else:
            full_file = cf

        all_files.append(full_file)

        f = open(full_file, "r")
        
        for line in f:

            # strip off the comments
            idx = line.find("!")
            line = line[:idx]

            rebreak = module_re.search(line)
            rebreak2 = module_proc_re.search(line)
            if rebreak and not rebreak2:
                modulefiles[rebreak.group(4)] = cf

        f.close()


    # go back through the files now and look for the use statements.
    # Assume only one use statement per line.  Ignore any only clauses.
    # Build a list of dependencies for the current file and output it.
    for cf in all_files:

        f = open(cf, "r")

        for line in f:

            # strip off the comments
            idx = line.find("!")
            line = line[:idx]

            rebreak = use_re.search(line)
            if rebreak:
                used_module = rebreak.group(4)
                if used_module in IGNORES:
                    continue
                print prefix+os.path.basename(cf).replace(".f90", ".o"), ':', \
                    prefix+os.path.basename(modulefiles[used_module]).replace(".f90", ".o")

        f.close()
        print " "

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",
                        help="prefix to prepend to each dependency pair, e.g., for a build directory",
                        default="")
    parser.add_argument("--search_path",
                        help="ordered path to search for the files",
                        default="")
    parser.add_argument("files", metavar="source files", type=str, nargs="*",
                        help="F90 source files to find dependencies amongst")

    args = parser.parse_args()

    prefix = "{}/".format(os.path.normpath(args.prefix))

    doit(prefix, args.search_path.split(), args.files)



