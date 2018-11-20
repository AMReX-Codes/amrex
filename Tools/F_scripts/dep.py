#!/usr/bin/env python

# automatically generate Makefile dependencies for Fortran 90 source.
#
# this will output all the dependency pairs amongst the source files.
#
# tips for debugging:
#
#   -- do `make > out` when building your code.  Any error messages
#      about modules not found will be sent to stderr, so you will
#      easily see them.  Usually this is a case of the regex not
#      capturing how the module is defined
#
#   -- add `--debug` to the `dep.py` commandline (e.g., in
#      `C_mk/Make.rules`).  This will output a file called
#      `dependencies.out` that shows details about what files are
#      parsed, what modules they define, and what modules they require.
#
#   -- a system-provided module won't be found in your file.
#      (e.g. iso_c_binding).  Add any system-provided modules to the
#      `IGNORES` list below

from __future__ import print_function

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for dep.py")

if sys.version[0] == "2":
    reload(sys)
    sys.setdefaultencoding('latin-1')


import io
import re
import os
import argparse
import preprocess

# modules to ignore in the dependencies
IGNORES = ["iso_c_binding", "iso_fortran_env", "omp_lib", "mpi", "cudafor", "openacc"]

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


class SourceFile(object):
    """ hold information about one of the .f90/.F90 files """

    def __init__(self, filename):

        self.name = filename

        # do we need to be preprocessed?  We'll use the convention
        # that .F90 = yes, .f90 = no
        self.ext = os.path.splitext(filename)[1]

        self.preprocess = bool(self.ext in [".F90", ".F95", ".F03"])

        # when we preprocess, the output file has a different name
        self.cpp_name = None


    def search_name(self):
        """return the file name we use for searching -- this is the
        preprocessed file if it exists"""

        if self.cpp_name is not None:
            search_file = self.cpp_name
        else:
            search_file = self.name

        return search_file


    def obj(self):
        """ the name of the object file we expect to be produced -- this
        will always be based on the original name -- we do not compile the
        preprocessed files """
        return self.name.replace(self.ext, ".o")


    def defined_modules(self):
        """determine what modules this file provides -- we work off of the
        preprocessed file if it exists."""

        defines = []

        with io.open(self.search_name(), "r", encoding="latin-1") as f:

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

        with io.open(self.search_name(), "r", encoding="latin-1") as f:

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


def doit(prefix, search_path, files, cpp, debug=False):
    """ main routine that processes the files"""

    if debug:
        df = io.open("dependencies.out", "w", encoding="latin-1")


    # module_files is a dictionary where the keys are the name of the
    # module (as it appears in Fortran code) and the value associated
    # with the key is the name of the object file that provides the
    # module.
    module_files = {}

    # all_files is a list of SourceFile objects
    all_files = []

    # find the locations of all the files, given an (optional) search
    # path and preprocess them if necessary
    for cf in files:

        # find the file in the first part of the search path it exists
        if len(search_path) > 0:
            for p in search_path:
                full_file = "{}/{}".format(p, cf)
                if os.path.isfile(full_file):
                    break
        else:
            full_file = cf

        sf = SourceFile(full_file)

        # preprocess, if necessary
        if sf.preprocess and cpp is not None:
            command = cpp.preprocess(sf)

        if debug:
            df.write("source file: {}\n".format(sf.name))
            if sf.preprocess:
                df.write("preprocessed: {}\n".format(sf.cpp_name))
            df.write("\n")

        all_files.append(sf)

    # for each file, figure out what modules they define and add those to
    # the module_files dictionary -- the module provided is the "key" and
    # the value is the object file corresponding to the source file that
    # holds the module definition
    for sf in all_files:
        provides = sf.defined_modules()

        for p in provides:
            module_files[p] = sf.obj()

        if debug:
            df.write("source file: {}\n".format(sf.name))
            df.write("provides:\n")
            for p in provides:
                df.write("   {}\n".format(p))
            df.write("\n")

    # go back through the files now and look for the use statements
    # to figure out what modules each file requires
    for sf in all_files:
        depends = sf.needed_modules()

        for d in depends:
            try: provides_obj = module_files[d]
            except KeyError:
                print("warning: module {} required by {} not found".format(d, sf.name),
                      file=sys.stderr)
                print("$(warning module {} required by {} not found)".format(d, sf.name))
                continue

            # skip the case where a file provides the module it needs
            # on its own; otherwise output the dependency line
            if provides_obj != sf.obj():
                print("{}: {}".format(
                    prefix+os.path.basename(sf.obj()),
                    prefix+os.path.basename(provides_obj)))

        print(" ")

        if debug:
            df.write("source file: {}\n".format(sf.name))
            df.write("depends on:\n")
            for d in depends:
                df.write("   {}\n".format(d))
            df.write("\n")

    if debug:
        df.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",
                        help="prefix to prepend to each dependency pair, e.g., for a build directory",
                        default="")
    parser.add_argument("--search_path",
                        help="ordered path to search for the files",
                        default="")
    parser.add_argument("--cpp",
                        help="command to run C preprocessor on .F90 files first.  If omitted, then no preprocessing is done",
                        default="")
    parser.add_argument("--f90_preprocess",
                        help="command to pipe cpp output to for additional Fortran preprocessing",
                        default="")
    parser.add_argument("--temp_dir",
                        help="directory to put temporary preprocessed files",
                        default="")
    parser.add_argument("--defines",
                        help="defines to send to preprocess the files",
                        default="")
    parser.add_argument("--debug",
                        help="output a detailed log file describing each source file",
                        action="store_true")
    parser.add_argument("files", metavar="source files", type=str, nargs="*",
                        help="F90 source files to find dependencies amongst")

    args = parser.parse_args()

    if args.prefix != "":
        prefix_pass = "{}/".format(os.path.normpath(args.prefix))
    else:
        prefix_pass = "./"

    if args.temp_dir != "":
        temp_dir = args.temp_dir
    else:
        temp_dir = "./"

    # create a preprocessor object
    if args.cpp != "":
        cpp_pass = preprocess.Preprocessor(temp_dir=temp_dir, cpp_cmd=args.cpp,
                                           defines=args.defines,
                                           f90_preprocess=args.f90_preprocess)
    else:
        cpp_pass = None

    try:
        doit(prefix_pass, args.search_path.split(), args.files, cpp_pass, debug=args.debug)
    except:
        # something went wrong
        print("$(error something went wrong in dep.py.  Remake, adding the option 'DEP_CHECK_OPTS=--debug' to your make command and examine the 'dependencies.out' file)")
