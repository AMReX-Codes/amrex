#!/usr/bin/env python

# automatically generate Makefile dependencies for Fortran 90 source.
#
# this will output all the dependency pairs amongst the source files.

from __future__ import print_function

import sys
import re
import os
import argparse


# modules to ignore in the dependencies
IGNORES = ["iso_c_binding", "iso_fortran_env"]

module_re = re.compile("( )(module)(\s+)((?:[a-z][a-z_0-9]+))",
                       re.IGNORECASE|re.DOTALL)

module_proc_re = re.compile("( )(module)(\s+)(procedure)(\s+)((?:[a-z][a-z_0-9]+))",
                            re.IGNORECASE|re.DOTALL)

# regular expression for ' use modulename, only: stuff, other stuff'
# see (txt2re.com)
use_re = re.compile("( )(use)(\s+)((?:[a-z_][a-z_0-9]+))", 
                    re.IGNORECASE|re.DOTALL)


class Preprocessor(object):
    """ hold the information about preprocessing """

    def __init__(self, temp_dir, cpp_cmd, includes, defines):
        self.temp_dir = temp_dir
        self.cpp_cmd = cpp_cmd
        self.includes = includes
        self.defines = defines

    def preprocess(self, source_file):
        """ preprocess the file described by a SourceFile object source_file """
        pass
        

class SourceFile(object):
    """ hold information about one of the .f90/.F90 files """

    def __init__(self, filename):

        self.name = filename

        # do we need to be preprocessed?  We'll use the convention
        # that .F90 = yes, .f90 = no
        ext = os.path.splitext(filename)

        if ext in [".F90", ".F95", ".F03"]:
            self.preprocess = True
        else:
            self.preprocess = False

        # when we preprocess, the output file has a different name
        self.cpp_name = None


    def defined_modules(self):
        """ determine what modules this file provides """

        defines = []

        with open(self.name, "r") as f:

            for line in f:

                # strip off the comments
                idx = line.find("!")
                line = line[:idx]

                # we want a module definition itself, not a 'module procedure'
                rebreak = module_re.search(line)
                rebreak2 = module_proc_re.search(line)
                if rebreak and not rebreak2:
                    defines.append(rebreak.group(4))

        return defines


    def needed_modules(self):
        """determine what modules this file needs.  Assume only one use
           statement per line.  Ignore any only clauses. """

        depends = []

        with open(self.name, "r") as f:

            for line in f:

                # strip off the comments
                idx = line.find("!")
                line = line[:idx]

                rebreak = use_re.search(line)
                if rebreak:
                    used_module = rebreak.group(4)
                    if used_module not in IGNORES:
                        depends.append(used_module)

        # remove duplicates
        depends = list(set(depends))

        return depends


def doit(prefix, search_path, files):

    # first parse the files and find all the module statements.  Keep a
    # dictionary of 'module name':filename.
    module_files = {}

    all_files = []

    # find the locations of all the files, given an (optional) search
    # path
    for cf in files:
        
        # find the file in the first part of the search path it exists
        if len(search_path) > 0:
            for p in search_path:
                full_file = "{}/{}".format(p, cf)
                if os.path.isfile(full_file): 
                    break
        else:
            full_file = cf

        all_files.append(SourceFile(full_file))

    # for each file, figure out what modules they define and add those to
    # the module_files dictionary -- the module provided is the "key"
    for sf in all_files:
        provides = sf.defined_modules()

        for p in provides:
            module_files[p] = sf.name


    # go back through the files now and look for the use statements
    # to figure out what modules each file requires
    for sf in all_files:
        depends = sf.needed_modules()
        
        for d in depends:
            print("{}: {}".format(
                prefix+os.path.basename(sf.name).replace(".f90", ".o"), 
                prefix+os.path.basename(module_files[d]).replace(".f90", ".o")))

        print(" ")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",
                        help="prefix to prepend to each dependency pair, e.g., for a build directory",
                        default="")
    parser.add_argument("--search_path",
                        help="ordered path to search for the files",
                        default="")
    parser.add_argument("--cpp",
                        help="apply the C preprocessor first?",
                        action="store_true")
    parser.add_argument("--temp_dir",
                        help="directory to put temporary preprocessed files",
                        default="")
    parser.add_argument("--includes",
                        help="include files needed to preprocess the source",
                        default="")
    parser.add_argument("--defines",
                        help="defines to send to preprocess the files",
                        default="")
    parser.add_argument("files", metavar="source files", type=str, nargs="*",
                        help="F90 source files to find dependencies amongst")

    args = parser.parse_args()

    prefix = "{}/".format(os.path.normpath(args.prefix))

    doit(prefix, args.search_path.split(), args.files)
