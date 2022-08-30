#!/usr/bin/env python3

# a simple script that writes the build_info.f90 file that is used
# to store information for the job_info file that we store in plotfiles.

source_string = """
module build_info_module

  implicit none

  character (len=128), save :: build_date = &
"@@BUILD_DATE@@"

  character (len=128), save :: build_dir = &
"@@BUILD_DIR@@"

  character (len=128), save :: build_machine = &
"@@BUILD_MACHINE@@"

  character (len=128), save :: amrex_dir = &
"@@AMREX_DIR@@"

  character (len=128), save :: FCOMP = &
"@@FCOMP@@"

  character (len=128), save :: FCOMP_version = &
"@@FCOMP_VERSION@@"

  character (len=250), save :: f90_compile_line = &
@@F90_COMP_LINE@@

  character (len=250), save :: f_compile_line = &
@@F_COMP_LINE@@

  character (len=250), save :: C_compile_line = &
@@C_COMP_LINE@@

  character (len=250), save :: link_line = &
@@LINK_LINE@@


  character (len=128), save :: amrex_git_hash = &
"@@AMREX_HASH@@"

  character (len=128), save :: source_git_hash = &
"@@SOURCE_HASH@@"

  character (len=128), save :: extra_git_hash = &
"@@EXTRA_HASH@@"

  character (len=128), save :: extra_git_hash2 = &
"@@EXTRA_HASH2@@"

  character (len=128), save :: network_dir = &
"@@NETWORK@@"

  character (len=128), save :: integrator_dir = &
"@@INTEGRATOR@@"

  character (len=128), save :: eos_dir = &
"@@EOS@@"

  character (len=128), save :: conductivity_dir = &
"@@CONDUCTIVITY@@"

  logical, parameter :: different_build_tree = @@BUILD_TREE_LOGICAL@@
  character (len=128), save :: build_git_hash = &
"@@BUILD_HASH@@"

  @@MODULE_STUFF@@

end module build_info_module
"""

module_str = """
  integer, parameter :: NUM_MODULES=@@NUM_MODULES@@
  character (len=128), save :: modules(NUM_MODULES) = (/ &
@@MODULE_INFO@@
  /)
"""


import os
import argparse
import datetime
import subprocess

def runcommand(command):
    p = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    out = p.stdout.read().decode("utf-8")
    return out.strip()


def get_git_hash(d):
    cwd = os.getcwd()
    os.chdir(d)
    try: ghash = runcommand("git rev-parse HEAD")
    except: ghash = ""
    os.chdir(cwd)
    return ghash


usage = """
This script is intended to be used within a F90 AMReX makefile to create
a file called build_info.f90 with information about the build environment.
"""

def doit():

    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--modules", type=str, default="",
                        metavar="'list of modules'",
                        help="the list of modules that were used to build the application")
    parser.add_argument("--FCOMP", type=str, default="",
                        metavar="fortran-compiler",
                        help="the name of the Fortran compiler executable")
    parser.add_argument("--FCOMP_version", type=str, default="",
                        metavar="fortran-compiler-version",
                        help="the Fortran compiler version number, as output by the compiler itself")
    parser.add_argument("--f90_compile_line", default="",
                        metavar="f90-compiler-invocation",
                        help="the Fortran 90 full compiler invocation used to build the application, including all options")
    parser.add_argument("--f_compile_line", default="",
                        metavar="f77-compiler-invocation",
                        help="the Fortran 77 full compiler invocation used to build the application, including all options")
    parser.add_argument("--C_compile_line", default="",
                        metavar="C-compiler-invocation",
                        help="the C full compiler invocation used to build the application, including all options")
    parser.add_argument("--link_line", default="",
                        metavar="link-invocation",
                        help="the link invocation used to link the application, including all options")
    parser.add_argument("--amrex_home", default="",
                        metavar="/path/to/AMReX",
                        help="the full path to the main amrex/ directory")
    parser.add_argument("--source_home", default="",
                        metavar="/path/to/source",
                        help="the full path to the main application source directory")
    parser.add_argument("--extra_home", default="",
                        metavar="/path/to/extra-home",
                        help="the full path to an additional source directory needed to build the application")
    parser.add_argument("--extra_home2", default="",
                        metavar="/path/to/extra-home2",
                        help="the full path to an additional source directory needed to build the application")
    parser.add_argument("--network", default="", metavar="network-name",
                        help="the name of any reaction network used in the application")
    parser.add_argument("--integrator", default="", metavar="integrator-name",
                        help="the name of any integration method used in the application")
    parser.add_argument("--eos", default="", metavar="eos-name",
                        help="the name of any equation of state used in the application")
    parser.add_argument("--conductivity", default="",
                        metavar="conductivity-name",
                        help="the name of any conductivity routine used in the application")

    args = parser.parse_args()

    MAX_STR_LEN = 128
    DBL_STR_LINE_LEN = 125

    # assemble some information

    # build stuff
    build_date = str(datetime.datetime.now())
    build_dir = os.getcwd()
    build_machine = runcommand("uname -a")

    # git hashes
    running_dir = os.getcwd()

    amrex_hash = get_git_hash(args.amrex_home)
    source_hash = get_git_hash(args.source_home)

    if not args.extra_home == "":
        try: os.chdir(args.extra_home)
        except:
            extra_hash = "ERROR: directory not found"
        else:
            extra_hash = get_git_hash(args.extra_home)
            os.chdir(running_dir)

    if not args.extra_home2 == "":
        try: os.chdir(args.extra_home2)
        except:
            extra_hash2 = "ERROR: directory not found"
        else:
            extra_hash2 = get_git_hash(args.extra_home2)
            os.chdir(running_dir)

    # we may not be building in a sub-directory of the source
    # directory, in that case, store an extra hash
    source_dir = os.path.abspath(args.source_home)
    build_dir = os.path.abspath(os.getcwd())

    sdir_parts = source_dir.split("/")
    bdir_parts = build_dir.split("/")

    is_sub_dir = 1
    for n in range(len(sdir_parts)):
        if not sdir_parts[n] == bdir_parts[n]:
            is_sub_dir = 0
            break

    if not is_sub_dir:
        have_build_hash = 1
        build_hash = get_git_hash(build_dir)
    else:
        have_build_hash = 0

    # modules
    module_list = args.modules.split()

    # output
    fout = open("build_info.f90", "w")

    for line in source_string.splitlines():

        index = line.find("@@")

        if index >= 0:
            index2 = line.rfind("@@")

            keyword = line[index+len("@@"):index2]

            if keyword == "BUILD_DATE":
                newline = line.replace("@@BUILD_DATE@@", build_date[:MAX_STR_LEN])
                fout.write(newline)

            elif keyword == "BUILD_DIR":
                newline = line.replace("@@BUILD_DIR@@", build_dir[:MAX_STR_LEN])
                fout.write(newline)

            elif keyword == "BUILD_MACHINE":
                newline = line.replace("@@BUILD_MACHINE@@",
                                       build_machine[:MAX_STR_LEN])
                fout.write(newline)

            elif keyword == "AMREX_DIR":
                newline = line.replace("@@AMREX_DIR@@",
                                       args.amrex_home[:MAX_STR_LEN])
                fout.write(newline)

            elif keyword == "FCOMP":
                newline = line.replace("@@FCOMP@@", args.FCOMP)
                fout.write(newline)

            elif keyword == "FCOMP_VERSION":
                newline = line.replace("@@FCOMP_VERSION@@",
                                       args.FCOMP_version[:MAX_STR_LEN])
                fout.write(newline)

            elif keyword == "F90_COMP_LINE":
                # this can span 2 lines
                if len(args.f90_compile_line) > DBL_STR_LINE_LEN:
                    tstr = args.f90_compile_line[:DBL_STR_LINE_LEN] + "\"// &\n\"" + \
                          args.f90_compile_line[DBL_STR_LINE_LEN:2*DBL_STR_LINE_LEN]
                else:
                    tstr = args.f90_compile_line

                newline = line.replace("@@F90_COMP_LINE@@",
                                       "\"%s\"" % (tstr))
                fout.write(newline)

            elif keyword == "F_COMP_LINE":
                # this can span 2 lines
                if len(args.f_compile_line) > DBL_STR_LINE_LEN:
                    tstr = args.f_compile_line[:DBL_STR_LINE_LEN] + "\"// &\n\"" + \
                          args.f_compile_line[DBL_STR_LINE_LEN:2*DBL_STR_LINE_LEN]
                else:
                    tstr = args.f_compile_line

                newline = line.replace("@@F_COMP_LINE@@",
                                       "\"%s\"" % (tstr))
                fout.write(newline)

            elif keyword == "C_COMP_LINE":
                # this can span 2 lines
                if len(args.C_compile_line) > DBL_STR_LINE_LEN:
                    tstr = args.C_compile_line[:DBL_STR_LINE_LEN] + "\"// &\n\"" + \
                          args.C_compile_line[DBL_STR_LINE_LEN:2*DBL_STR_LINE_LEN]
                else:
                    tstr = args.C_compile_line

                newline = line.replace("@@C_COMP_LINE@@",
                                       "\"%s\"" % (tstr))
                fout.write(newline)

            elif keyword == "LINK_LINE":
                # this can span 2 lines
                if len(args.link_line) > DBL_STR_LINE_LEN:
                    tstr = args.link_line[:DBL_STR_LINE_LEN] + "\"// &\n\"" + \
                          args.link_line[DBL_STR_LINE_LEN:2*DBL_STR_LINE_LEN]
                else:
                    tstr = args.link_line

                newline = line.replace("@@LINK_LINE@@",
                                       "\"%s\"" % (tstr))
                fout.write(newline)

            elif keyword == "AMREX_HASH":
                newline = line.replace("@@AMREX_HASH@@", amrex_hash)
                fout.write(newline)

            elif keyword == "SOURCE_HASH":
                newline = line.replace("@@SOURCE_HASH@@", source_hash)
                fout.write(newline)

            elif keyword == "EXTRA_HASH":
                if not args.extra_home == "":
                    newline = line.replace("@@EXTRA_HASH@@", extra_hash)
                else:
                    newline = line.replace("@@EXTRA_HASH@@", "")

                fout.write(newline)

            elif keyword == "EXTRA_HASH2":
                if not args.extra_home2 == "":
                    newline = line.replace("@@EXTRA_HASH2@@", extra_hash2)
                else:
                    newline = line.replace("@@EXTRA_HASH2@@", "")

                fout.write(newline)

            elif keyword == "BUILD_TREE_LOGICAL":
                if have_build_hash == 1:
                    newline = line.replace("@@BUILD_TREE_LOGICAL@@",
                                           ".true.")
                else:
                    newline = line.replace("@@BUILD_TREE_LOGICAL@@",
                                           ".false.")

                fout.write(newline)

            elif keyword == "BUILD_HASH":
                if have_build_hash == 1:
                    newline = line.replace("@@BUILD_HASH@@", build_hash)
                else:
                    newline = line.replace("@@BUILD_HASH@@", "")

                fout.write(newline)

            elif keyword == "NETWORK":
                fout.write(line.replace("@@NETWORK@@", args.network))

            elif keyword == "INTEGRATOR":
                fout.write(line.replace("@@INTEGRATOR@@", args.integrator))

            elif keyword == "EOS":
                fout.write(line.replace("@@EOS@@", args.eos))

            elif keyword == "CONDUCTIVITY":
                fout.write(line.replace("@@CONDUCTIVITY@@", args.conductivity))

            elif keyword == "MODULE_STUFF":

                if len(module_list) == 0:
                    newline = line.replace("@@MODULE_STUFF@@",
                                           "integer, parameter :: NUM_MODULES=0")
                else:
                    tstr = ""
                    for n in range(len(module_list)):
                        if n < len(module_list)-1:
                            tstr += "\"%-120s\", &\n" % (module_list[n])
                        else:
                            tstr += "\"%-120s\" &" % (module_list[n])

                    newlinet = module_str.replace("@@NUM_MODULES@@", repr(len(module_list)))
                    newline = newlinet.replace("@@MODULE_INFO@@", tstr)
                fout.write(newline)

        else:
            fout.write(line)

        fout.write("\n")

    fout.close()


if __name__ == "__main__":
    doit()
