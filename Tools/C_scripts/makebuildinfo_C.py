#!/usr/bin/env python

import sys

if sys.version_info < (2, 7):
    sys.exit("ERROR: need python 2.7 or later for dep.py")

import argparse

import os
import datetime
import subprocess


source = """
namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = "@@BUILD_DATE@@";
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "@@BUILD_DIR@@";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "@@BUILD_MACHINE@@";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "@@amrex_home@@";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "@@COMP@@";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "@@COMP_VERSION@@";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "@@FCOMP@@";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "@@FCOMP_VERSION@@";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "@@CXX_comp_name@@";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "@@F_comp_name@@";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = "@@CXX_flags@@";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = "@@F_flags@@";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = "@@link_flags@@";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = "@@libraries@@";
  return libraries;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";
  @@AUX_DECLS@@
  static const char EMPT[] = "";

  switch(i)
  {
    @@AUX_CASE@@
    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int num_modules = X;
  @@NUM_MODULES@@
  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";
  @@MNAME_DECLS@@
  static const char EMPT[] = "";

  switch(i)
  {
    @@MNAME_CASE@@
    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";
  @@MVAL_DECLS@@
  static const char EMPT[] = "";

  switch(i)
  {
    @@MVAL_CASE@@
    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  @@GIT_DECLS@@
  static const char EMPT[] = "";

  switch(i)
  {
    @@GIT_CASE@@
    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  @@BUILDGIT_DECLS@@

  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  @@BUILDGIT_NAME@@

  return NAME;
}

}
"""

def runcommand(command):
    p = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    out = p.stdout.read()
    return out.strip().decode("ascii")

def get_git_hash(d):
    cwd = os.getcwd()
    os.chdir(d)
    try: hash = runcommand("git describe --always --tags --dirty")
    except: hash = ""
    os.chdir(cwd)
    return hash


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--amrex_home", help="path to the AMReX source",
                        type=str, default="")

    parser.add_argument("--COMP",
                        help="Compiler system defined in AMReX's build system",
                        type=str, default="")

    parser.add_argument("--COMP_VERSION", help="Compiler version",
                        type=str, default="")

    parser.add_argument("--CXX_comp_name",
                        help="C++ compiler command", type=str, default="")

    parser.add_argument("--CXX_flags",
                        help="C++ compiler flags", type=str, default="")

    parser.add_argument("--FCOMP",
                        help="Fortran compiler as defined by AMReX's build system (deprecated)",
                        type=str, default="")

    parser.add_argument("--FCOMP_VERSION",
                        help="Fortran compiler version (deprecated)",
                        type=str, default="")

    parser.add_argument("--F_comp_name", help="Fortran compiler command",
                        type=str, default="")

    parser.add_argument("--F_flags", help="Fortran compiler flags",
                        type=str, default="")

    parser.add_argument("--link_flags", help="linker flags", type=str, default="")

    parser.add_argument("--libraries", help="libraries linked", type=str, default="")

    parser.add_argument("--AUX",
                        help="auxillary information (EOS, network path) (deprecated)",
                        type=str, default="")

    parser.add_argument("--MODULES",
                       help="module information in the form of key=value (e.g., EOS=helmeos)",
                       type=str, default="")

    parser.add_argument("--GIT",
                        help="the directories whose git hashes we should capture",
                        type=str, default="")


    parser.add_argument("--build_git_name",
                        help="the name of the build directory if different from the main source -- we'll get the git hash of this.",
                        type=str, default="")

    parser.add_argument("--build_git_dir",
                        help="the full path to the build directory that corresponds to build_git_name",
                        type=str, default="")


    # parse and convert to a dictionary
    args = parser.parse_args()

    # build stuff
    build_date = str(datetime.datetime.now())
    build_dir = os.getcwd()
    build_machine = runcommand("uname -a")


    # git hashes
    running_dir = os.getcwd()

    if args.GIT == "":
        GIT = []
    else:
        GIT = args.GIT.split()

    ngit = len(GIT)
    git_hashes = []
    for d in GIT:
        if d and os.path.isdir(d):
            git_hashes.append(get_git_hash(d))
        else:
            git_hashes.append("")

    if args.build_git_dir != "":
        try: os.chdir(args.build_git_dir)
        except:
            build_git_hash = "directory not valid"
        else:
            build_git_hash = get_git_hash(args.build_git_dir)
            os.chdir(running_dir)
    else:
        build_git_hash = ""


    # modules
    if args.MODULES == "":
        MODULES = []
    else:
        MODULES = args.MODULES.split()

    if len(MODULES) > 0:
        mod_dict = {}
        for m in MODULES:
            k, v = m.split("=")
            mod_dict[k] = v


    # auxillary info
    if args.AUX == "":
        AUX = []
    else:
        AUX = args.AUX.split()


    fout = open("AMReX_buildInfo.cpp", "w")

    # dictionary view of the args
    dargs = vars(args)

    for line in source.splitlines():

        index = line.find("@@")

        if index >= 0:
            index2 = line.rfind("@@")
            keyword = line[index+len("@@"):index2]

            if keyword == "BUILD_DATE":
                newline = line.replace("@@BUILD_DATE@@", build_date)
                fout.write(newline)

            elif keyword == "BUILD_DIR":
                newline = line.replace("@@BUILD_DIR@@", build_dir)
                fout.write(newline)

            elif keyword == "BUILD_MACHINE":
                newline = line.replace("@@BUILD_MACHINE@@", build_machine)
                fout.write(newline)

            elif keyword == "AUX_DECLS":
                indent = index
                aux_str = ""
                for n, a in enumerate(AUX):
                    aux_str += '{}static const char AUX{:1d}[] = "{}";\n'.format(
                        indent*" ", n+1, a)

                fout.write(aux_str)

            elif keyword == "AUX_CASE":
                indent = index
                aux_str = ""
                for i in range(len(AUX)):
                    aux_str += '{}case {:1d}: return AUX{:1d};\n'.format(
                        indent*" ", i+1, i+1)

                fout.write(aux_str)

            elif keyword == "NUM_MODULES":
                num_modules = len(MODULES)
                indent = index
                fout.write("{}int num_modules = {};\n".format(
                    indent*" ", num_modules))

            elif keyword == "MNAME_DECLS":
                indent = index
                aux_str = ""
                if len(MODULES) > 0:
                    for i, m in enumerate(list(mod_dict.keys())):
                        aux_str += '{}static const char AUX{:1d}[] = "{}";\n'.format(
                            indent*" ", i+1, m)

                fout.write(aux_str)

            elif keyword == "MNAME_CASE":
                indent = index
                aux_str = ""
                if len(MODULES) > 0:
                    for i, m in enumerate(list(mod_dict.keys())):
                        aux_str += '{}case {:1d}: return AUX{:1d};\n'.format(
                            indent*" ", i+1, i+1)

                fout.write(aux_str)

            elif keyword == "MVAL_DECLS":
                indent = index
                aux_str = ""
                if len(MODULES) > 0:
                    for i, m in enumerate(list(mod_dict.keys())):
                        aux_str += '{}static const char AUX{:1d}[] = "{}";\n'.format(
                            indent*" ", i+1, mod_dict[m])

                fout.write(aux_str)

            elif keyword == "MVAL_CASE":
                indent = index
                aux_str = ""
                if len(MODULES) > 0:
                    for i, m in enumerate(list(mod_dict.keys())):
                        aux_str += '{}case {:1d}: return AUX{:1d};\n'.format(
                            indent*" ", i+1, i+1)

                fout.write(aux_str)

            elif keyword == "GIT_DECLS":
                indent = index
                git_str = ""
                for i, gh in enumerate(git_hashes):
                    git_str += '{}static const char HASH{:1d}[] = "{}";\n'.format(
                        indent*" ", i+1, gh)

                fout.write(git_str)

            elif keyword == "GIT_CASE":
                indent = index
                git_str = ""
                for i in range(len(git_hashes)):
                    git_str += '{}case {:1d}: return HASH{:1d};\n'.format(
                        indent*" ", i+1, i+1)

                fout.write(git_str)

            elif keyword == "BUILDGIT_DECLS":
                indent = index
                git_str = '{}static const char HASH[] = "{}";\n'.format(
                    indent*" ", build_git_hash)
                fout.write(git_str)

            elif keyword == "BUILDGIT_NAME":
                index = index
                git_str = '{}static const char NAME[] = "{}";\n'.format(
                    indent*" ", args.build_git_name)
                fout.write(git_str)

            elif keyword in dargs:
                # simple replacement using the commandline arguments
                newline = line.replace("@@{}@@".format(keyword),
                                       dargs[keyword].replace('"', r'\"'))
                fout.write(newline)

        else:
            fout.write(line)

        fout.write("\n")

    fout.close()
