#!/usr/bin/env python

# a simple script that writes the build_info.f90 file that is used
# to store information for the job_info file that we store in plotfiles.

sourceString="""
module build_info_module

  implicit none

  character (len=128), save :: build_date = &
"@@BUILD_DATE@@"

  character (len=128), save :: build_dir = &
"@@BUILD_DIR@@"

  character (len=128), save :: build_machine = &
"@@BUILD_MACHINE@@"

  character (len=128), save :: boxlib_dir = &
"@@BOXLIB_DIR@@"

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


  character (len=128), save :: boxlib_git_hash = &
"@@BOXLIB_HASH@@"

  character (len=128), save :: source_git_hash = &
"@@SOURCE_HASH@@"

  character (len=128), save :: extra_git_hash = &
"@@EXTRA_HASH@@"

  character (len=128), save :: extra_git_hash2 = &
"@@EXTRA_HASH2@@"

  logical, parameter :: different_build_tree = @@BUILD_TREE_LOGICAL@@
  character (len=128), save :: build_git_hash = &
"@@BUILD_HASH@@"


  integer, parameter :: NUM_MODULES=@@NUM_MODULES@@
  character (len=128), save :: modules(NUM_MODULES) = (/ &
@@MODULE_INFO@@
  /)

end module build_info_module
"""

import sys
import os
import getopt
import datetime
import subprocess

def runcommand(command):
    p = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    out = p.stdout.read().decode("utf-8")
    return out.strip()


def get_git_hash(d):
    cwd = os.getcwd()
    os.chdir(d)
    try: hash = runcommand("git rev-parse HEAD")
    except: hash = ""
    os.chdir(cwd)
    return hash

try: opts, next = getopt.getopt(sys.argv[1:], "",
                               ["modules=",
                                "FCOMP=",
                                "FCOMP_version=",
                                "f90_compile_line=",
                                "f_compile_line=",
                                "C_compile_line=",
                                "link_line=",
                                "boxlib_home=",
                                "source_home=",
                                "extra_home=",
                                "extra_home2="])
except getopt.GetoptError:
    sys.exit("invalid calling sequence")

modules = ""
FCOMP = ""
FCOMP_version = ""
f90_compile_line = ""
f_compile_line = ""
C_compile_line = ""
link_line = ""
boxlib_home = ""
source_home = ""
extra_home = ""
extra_home2 = ""

for o, a in opts:

    if o == "--modules":
        modules = a

    if o == "--FCOMP":
        FCOMP = a

    if o == "--FCOMP_version":
        FCOMP_version = a

    if o == "--f90_compile_line":
        f90_compile_line = a

    if o == "--f_compile_line":
        f_compile_line = a

    if o == "--C_compile_line":
        C_compile_line = a
        
    if o == "--link_line":
        link_line = a

    if o == "--boxlib_home":
        boxlib_home = a

    if o == "--source_home":
        source_home = a

    if o == "--extra_home":
        extra_home = a

    if o == "--extra_home2":
        extra_home2 = a        


MAX_STRING_LENGTH=128
DBL_STRING_LINE_LENGTH=125

# assemble some information

# build stuff
build_date = str(datetime.datetime.now())
build_dir = os.getcwd()
build_machine = runcommand("uname -a")

# git hashes
runningDir = os.getcwd()

boxlib_hash = get_git_hash(boxlib_home)
source_hash = get_git_hash(source_home)

if (not extra_home == ""):
    try: os.chdir(extra_home)
    except:
        extra_hash = "ERROR: directory not found"
    else:
        extra_hash = get_git_hash(extra_home)
        os.chdir(runningDir)

if (not extra_home2 == ""):
    try: os.chdir(extra_home2)
    except:
        extra_hash2 = "ERROR: directory not found"
    else:
        extra_hash2 = get_git_hash(extra_home2)
        os.chdir(runningDir)
    
# we may not be building in a sub-directory of the source directory, in that
# case, store an extra hash
sourceDir = os.path.abspath(source_home)
buildDir = os.path.abspath(os.getcwd())

sDirParts = sourceDir.split("/")
bDirParts = buildDir.split("/")

isSubDir = 1
for n in range(len(sDirParts)):
    if not sDirParts[n] == bDirParts[n]:
        isSubDir = 0
        break

if not isSubDir:
    have_build_hash = 1
    build_hash = get_git_hash(buildDir)
else:
    have_build_hash = 0

# modules
moduleList = modules.split()


# output
fout = open("build_info.f90", "w")

for line in sourceString.splitlines():

    index = line.find("@@")

    if index >= 0:
        index2 = line.rfind("@@")

        keyword = line[index+len("@@"):index2]

        if keyword == "BUILD_DATE":
            newline = line.replace("@@BUILD_DATE@@", build_date[:MAX_STRING_LENGTH])
            fout.write(newline)

        elif keyword == "BUILD_DIR":
            newline = line.replace("@@BUILD_DIR@@", build_dir[:MAX_STRING_LENGTH])
            fout.write(newline)

        elif keyword == "BUILD_MACHINE":
            newline = line.replace("@@BUILD_MACHINE@@", 
                                   build_machine[:MAX_STRING_LENGTH])
            fout.write(newline)

        elif keyword == "BOXLIB_DIR":
            newline = line.replace("@@BOXLIB_DIR@@", 
                                   boxlib_home[:MAX_STRING_LENGTH])
            fout.write(newline)

        elif keyword == "FCOMP":
            newline = line.replace("@@FCOMP@@", 
                                   FCOMP)
            fout.write(newline)            

        elif keyword == "FCOMP_VERSION":
            newline = line.replace("@@FCOMP_VERSION@@", 
                                   FCOMP_version[:MAX_STRING_LENGTH])
            fout.write(newline)            

        elif keyword == "F90_COMP_LINE":
            # this can span 2 lines
            if (len(f90_compile_line) > DBL_STRING_LINE_LENGTH):
                str = f90_compile_line[:DBL_STRING_LINE_LENGTH] + "\"// &\n\"" + \
                      f90_compile_line[DBL_STRING_LINE_LENGTH:2*DBL_STRING_LINE_LENGTH]
            else:
                str = f90_compile_line

            newline = line.replace("@@F90_COMP_LINE@@", 
                                   "\"%s\"" % (str))
            fout.write(newline)

        elif keyword == "F_COMP_LINE":
            # this can span 2 lines
            if (len(f_compile_line) > DBL_STRING_LINE_LENGTH):
                str = f_compile_line[:DBL_STRING_LINE_LENGTH] + "\"// &\n\"" + \
                      f_compile_line[DBL_STRING_LINE_LENGTH:2*DBL_STRING_LINE_LENGTH]
            else:
                str = f_compile_line

            newline = line.replace("@@F_COMP_LINE@@", 
                                   "\"%s\"" % (str))
            fout.write(newline)

        elif keyword == "C_COMP_LINE":
            # this can span 2 lines
            if len(C_compile_line) > DBL_STRING_LINE_LENGTH:
                str = C_compile_line[:DBL_STRING_LINE_LENGTH] + "\"// &\n\"" + \
                      C_compile_line[DBL_STRING_LINE_LENGTH:2*DBL_STRING_LINE_LENGTH]
            else:
                str = C_compile_line

            newline = line.replace("@@C_COMP_LINE@@", 
                                   "\"%s\"" % (str))
            fout.write(newline)

        elif keyword == "LINK_LINE":
            # this can span 2 lines
            if len(link_line) > DBL_STRING_LINE_LENGTH:
                str = link_line[:DBL_STRING_LINE_LENGTH] + "\"// &\n\"" + \
                      link_line[DBL_STRING_LINE_LENGTH:2*DBL_STRING_LINE_LENGTH]
            else:
                str = link_line

            newline = line.replace("@@LINK_LINE@@", 
                                   "\"%s\"" % (str))
            fout.write(newline)

        elif keyword == "BOXLIB_HASH":
            newline = line.replace("@@BOXLIB_HASH@@", boxlib_hash)
            fout.write(newline)

        elif keyword == "SOURCE_HASH":
            newline = line.replace("@@SOURCE_HASH@@", source_hash)
            fout.write(newline)

        elif keyword == "EXTRA_HASH":
            if not extra_home == "":
                newline = line.replace("@@EXTRA_HASH@@", extra_hash)
            else:
                newline = line.replace("@@EXTRA_HASH@@", "")

            fout.write(newline)

        elif keyword == "EXTRA_HASH2":
            if not extra_home2 == "":
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

        elif keyword == "NUM_MODULES":
            newline = line.replace("@@NUM_MODULES@@", repr(len(moduleList)))
            fout.write(newline)

        elif keyword == "MODULE_INFO":
            str = ""
            for n in range(len(moduleList)):
                if n < len(moduleList)-1:
                    str += "\"%-120s\", &\n" % (moduleList[n])
                else:
                    str += "\"%-120s\" &" % (moduleList[n])

            newline = line.replace("@@MODULE_INFO@@", str)
            fout.write(newline)

    else:
        fout.write(line)

    fout.write("\n")

fout.close()
