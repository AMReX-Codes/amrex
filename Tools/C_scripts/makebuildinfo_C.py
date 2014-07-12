#!/usr/bin/env python

import sys
import os
import argparse
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
"""

def runcommand(command):
    p = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    out = p.stdout.read()
    return out.strip()


parser = argparse.ArgumentParser()
parser.add_argument("--boxlib_home", help="location of the BoxLib/ directory", 
                   type=str, default=None)
parser.add_argument("--FCOMP", help="the name of the Fortran compiler", 
                   type=str, default=None)
parser.add_argument("--COMP", help="the name of the C++ compiler", 
                   type=str, default=None)
parser.add_argument("--AUX", help="auxillary information (EOS, network path)", 
                   type=str, default=None, nargs="*")
parser.add_argument("--GIT", help="the directories whose git hashes we should capture", 
                   type=str, default=None, nargs="*")

args = parser.parse_args()

# build stuff                                                                                           
build_date = str(datetime.datetime.now())
build_dir = os.getcwd()
build_machine = runcommand("uname -a")

# git hashes
running_dir = os.getcwd()

ngit = len(args.GIT)
git_hashes = []
for i in range(ngit):
    os.chdir(args.GIT[i])
    git_hashes.append(runcommand("git rev-parse HEAD"))
    os.chdir(running_dir)


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
            newline = string.replace(line, "@@BOXLIB_DIR@@", args.boxlib_home)
            fout.write(newline)        

        elif keyword == "COMP":
            newline = string.replace(line, "@@COMP@@", args.COMP)
            fout.write(newline)        

        elif keyword == "FCOMP":
            newline = string.replace(line, "@@FCOMP@@", args.FCOMP)
            fout.write(newline)        

        elif keyword == "AUX_DECLS":
            indent = index
            aux_str = ""
            for i in range(len(args.AUX)):
                aux_str += '{}static const char AUX{}[] = "{}";\n'.format(indent*" ", i+1, args.AUX[i])

            fout.write(aux_str)        

        elif keyword == "AUX_CASE":
            indent = index
            aux_str = ""
            for i in range(len(args.AUX)):
                aux_str += '{}case {}: return AUX{};\n'.format(indent*" ", i+1, i+1)

            fout.write(aux_str)

        elif keyword == "GIT_DECLS":
            indent = index
            git_str = ""
            for i in range(len(args.GIT)):
                git_str += '{}static const char HASH{}[] = "{}";\n'.format(indent*" ", i+1, git_hashes[i])

            fout.write(git_str)        

        elif keyword == "GIT_CASE":
            indent = index
            git_str = ""
            for i in range(len(args.GIT)):
                git_str += '{}case {}: return HASH{};\n'.format(indent*" ", i+1, i+1)

            fout.write(git_str)

    else:
        fout.write(line)

    fout.write("\n")

fout.close()





