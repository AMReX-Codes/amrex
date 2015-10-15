#!/usr/bin/env python

import sys
import os
import getopt
import datetime
import string
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

const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "@@FCOMP@@";
  return FCOMP;
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
    return out.strip()

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
                                 "COMP=",
                                 "AUX=",
                                 "GIT=",
                                 "build_git_name=",
                                 "build_git_dir="])

except getopt.GetoptError:
    sys.exit("invalid calling sequence")

boxlib_home = ""
FCOMP = ""
COMP = ""
AUX = []
GIT = []
build_git_name = ""
build_git_dir = None

for o, a in opts:

    if o == "--boxlib_home": boxlib_home = a

    if o == "--FCOMP": FCOMP = a

    if o == "--COMP": COMP = a

    if o == "--AUX":
        if not a == "": AUX = a.split()

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
        build_git_hash = runcommand(get_git_hash(build_git_dir))
        os.chdir(running_dir)
else:
    build_git_hash = ""

fout = open("buildInfo.cpp", "w")

for line in source.splitlines():

    index = line.find("@@")

    if index >= 0:
        index2 = line.rfind("@@")
        keyword = line[index+len("@@"):index2]

        if keyword == "BUILD_DATE":
            newline = string.replace(line, "@@BUILD_DATE@@", build_date)
            fout.write(newline)

        elif keyword == "BUILD_DIR":
            newline = string.replace(line, "@@BUILD_DIR@@", build_dir)
            fout.write(newline)

        elif keyword == "BUILD_MACHINE":
            newline = string.replace(line, "@@BUILD_MACHINE@@", build_machine)
            fout.write(newline)

        elif keyword == "BOXLIB_DIR":
            newline = string.replace(line, "@@BOXLIB_DIR@@", boxlib_home)
            fout.write(newline)

        elif keyword == "COMP":
            newline = string.replace(line, "@@COMP@@", COMP)
            fout.write(newline)

        elif keyword == "FCOMP":
            newline = string.replace(line, "@@FCOMP@@", FCOMP)
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
