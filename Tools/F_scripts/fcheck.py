#!/usr/bin/env python

# a simple routine to parse Fortran files and make sure that things are 
# declared double precision, and constants are of the form 1.0_dp_t or 
# 1.0d0.

import os
import re
import sys

# list of valid Fortran file extensions
extensions = ['.f90', '.f']


#-----------------------------------------------------------------------------
# myFiles class
#-----------------------------------------------------------------------------
class myFiles:

    # a simple class to hold a list of all the files in the search path with
    # the given extension
    def __init__(self,extension=""):
        self.extension = extension
        self.files = []


#-----------------------------------------------------------------------------
# visit function
#-----------------------------------------------------------------------------
def visit(argFiles, dirname, files):

    # invoked by os.path.walk to find all the files
    # matching our 
    for file in files:
        base, ext = os.path.splitext(file)
        if (ext == argFiles.extension):
            currentFiles.files.append(dirname + '/' + file)


#-----------------------------------------------------------------------------
# main program 
#-----------------------------------------------------------------------------

# generate the regular expressions we will use to search
# helpful site: http://re-try.appspot.com/

# look for 'real' followed by any type of space character '\s', or 
# an open parenthesis
real_re = re.compile('real', re.IGNORECASE)

# Look for some for of real with 'dp_t'.  This should match
#   real(dp_t) :: X
#   real (dp_t) :: X
#   real (kind=dp_t) :: X
#   real ( kind = dp_t) :: X
#   ...
dpt_re = re.compile('real[\s]*\([\s]*(kind\s*=\s*)*dp_t[\s]*\)', re.IGNORECASE)


# look for a floating point constant
const_fp_re = re.compile('(\d+[.]\d*|\d*[.]\d+)', re.IGNORECASE)

# Look for constants declared as
# 1.0e0_dp_t
# 1.0d0_dp_t
# 1.0d0
# i.e. that have to be double precision somehow
const_dp_re = re.compile('\d*[.]\d*([ed][+-]?\d+_dp_t|d[+-]?\d+|_dp_t)', re.IGNORECASE)

# recursively find all of the files with a given extension
for ext in extensions:

    currentFiles = myFiles(extension=ext)

    root = os.getcwd()

    os.path.walk(root, visit, currentFiles)
    
    for file in currentFiles.files:

        # open the file for parsing
        try:
            fin = open(file, 'r')
        except (OSError, IOError), error:
            print error
            sys.exit()


        # read in each line, one by one
        lineNum = 1
        line = fin.readline()
        badFile = 0
        while (line):

            # for each line where we match the pattern real_re, make sure
            # that we also match the pattern dpt_re
            if (not real_re.search(line) == None and
                dpt_re.search(line) == None):
                if (badFile == 0):
                    index = file.find(root)
                    print file[index+len(root)+1:] + ':'
                    badFile = 1

                print lineNum, ": ", line,


            # if we find a floating point constant, make sure that it is 
            # double precision
            if (not const_fp_re.search(line) == None and
                const_dp_re.search(line) == None):
                if (badFile == 0):
                    index = file.find(root)
                    print file[index+len(root)+1:] + ':'
                    badFile = 1

                print lineNum, ": ", line,



            line = fin.readline()            
            lineNum += 1

        if (badFile == 1):
            print " "



    
    
    
