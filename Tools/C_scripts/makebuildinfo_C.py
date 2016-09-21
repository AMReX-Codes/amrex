#!/usr/bin/env python

import sys
import os
import getopt
import datetime
import subprocess


source = """
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

const char* buildInfoGetBoxlibDir() {

  static const char BOXLIB_DIR[] = "@@BOXLIB_DIR@@";
  return BOXLIB_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "@@COMP@@";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "@@COMP_VERSION@@";
  return COMP_VERSION;
}

const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "@@FCOMP@@";
  return FCOMP;
}

const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "@@FCOMP_VERSION@@";
  return FCOMP_VERSION;
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
"""

def runcommand(command):
    p = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    out = p.stdout.read()
    return out.strip().decode("ascii")

def get_git_hash(d):
    cwd = os.getcwd()
    os.chdir(d)
    try: hash = runcommand("git rev-parse HEAD")
    except: hash = ""
    os.chdir(cwd)
    return hash


try: opts, next = getopt.getopt(sys.argv[1:], "",
                                ["boxlib_home=",
                                 "FCOMP=",
                                 "FCOMP_VERSION=",
                                 "COMP=",
                                 "COMP_VERSION=",
                                 "AUX=",
                                 "MODULES=",
                                 "GIT=",
                                 "build_git_name=",
                                 "build_git_dir="])

except getopt.GetoptError:
    sys.exit("invalid calling sequence")

boxlib_home = ""
FCOMP = ""
FCOMP_VERSION = ""
COMP = ""
COMP_VERSION = ""
AUX = []
# modules have the form "key=value", e.g. "EOS=helmeos"
MODULES = []
GIT = []
build_git_name = ""
build_git_dir = None

for o, a in opts:

    if o == "--boxlib_home": boxlib_home = a

    if o == "--FCOMP": FCOMP = a

    if o == "--FCOMP_VERSION": FCOMP_VERSION = a    

    if o == "--COMP": COMP = a

    if o == "--COMP_VERSION": COMP_VERSION = a    

    if o == "--AUX":
        if not a == "": AUX = a.split()

    if o == "--MODULES":
        if not a == "": MODULES = a.split()

    if o == "--GIT":
        if not a == "": GIT = a.split()

    if o == "--build_git_name":
        if not a == "": build_git_name = a

    if o == "--build_git_dir":
        if not a == "": build_git_dir = a

        
# build stuff
build_date = str(datetime.datetime.now())
build_dir = os.getcwd()
build_machine = runcommand("uname -a")

# git hashes
running_dir = os.getcwd()

ngit = len(GIT)
git_hashes = []
for d in GIT:
    if d and os.path.isdir(d):
        git_hashes.append(get_git_hash(d))
    else:
        git_hashes.append("")

if not build_git_dir == None:
    try: os.chdir(build_git_dir)
    except:
        build_git_hash = "directory not valid"
    else:
        build_git_hash = get_git_hash(build_git_dir)
        os.chdir(running_dir)
else:
    build_git_hash = ""


if len(MODULES) > 0:
    mod_dict = {}
    for m in MODULES:
        k, v = m.split("=")
        mod_dict[k] = v

fout = open("buildInfo.cpp", "w")

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

        elif keyword == "BOXLIB_DIR":
            newline = line.replace("@@BOXLIB_DIR@@", boxlib_home)
            fout.write(newline)

        elif keyword == "COMP":
            newline = line.replace("@@COMP@@", COMP)
            fout.write(newline)

        elif keyword == "COMP_VERSION":
            newline = line.replace("@@COMP_VERSION@@", COMP_VERSION)
            fout.write(newline)

        elif keyword == "FCOMP":
            newline = line.replace("@@FCOMP@@", FCOMP)
            fout.write(newline)

        elif keyword == "FCOMP_VERSION":
            newline = line.replace("@@FCOMP_VERSION@@", FCOMP_VERSION)
            fout.write(newline)

        elif keyword == "AUX_DECLS":
            indent = index
            aux_str = ""
            for i in range(len(AUX)):
                aux_str += '%sstatic const char AUX%1d[] = "%s";\n' % (indent*" ", i+1, AUX[i])

            fout.write(aux_str)

        elif keyword == "AUX_CASE":
            indent = index
            aux_str = ""
            for i in range(len(AUX)):
                aux_str += '%scase %1d: return AUX%1d;\n' % (indent*" ", i+1, i+1)

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
            for i in range(len(GIT)):
                git_str += '%sstatic const char HASH%1d[] = "%s";\n' % (indent*" ", i+1, git_hashes[i])

            fout.write(git_str)

        elif keyword == "GIT_CASE":
            indent = index
            git_str = ""
            for i in range(len(GIT)):
                git_str += '%scase %1d: return HASH%1d;\n' % (indent*" ", i+1, i+1)

            fout.write(git_str)

        elif keyword == "BUILDGIT_DECLS":
            indent = index
            git_str = '%sstatic const char HASH[] = "%s";\n' % (indent*" ", build_git_hash)
            fout.write(git_str)

        elif keyword == "BUILDGIT_NAME":
            index = index
            git_str = '%sstatic const char NAME[] = "%s";\n' % (indent*" ", build_git_name)
            fout.write(git_str)
            
    else:
        fout.write(line)

    fout.write("\n")

fout.close()
