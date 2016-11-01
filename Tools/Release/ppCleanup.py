#!/usr/bin/env python

import os
import shutil
import sys
import getopt
import re

def ppCleanup(cwfile, infile, oufile):

    try:
        f = open(cwfile)
    except IOError:
        print "ERROR: ", cwfile, " doesn't appear to exist."
        sys.exit(2)
    
    cleanWords = []
    for line in f.readlines():
        word = line[:-1].strip()
        if word:
            cleanWords.append(line[:-1])
    f.close()

    try:
        fin = open(infile)
    except IOError:
        print "ERROR: ", infile, " doesn't appear to exist."
        sys.exit(2)

    try:
        fou = open(oufile, 'w')
    except:
        print "ERROR: in opening ", infile
        sys.exit(2)


    ifRE = re.compile(r"#if(\s|\().*\n", re.IGNORECASE)
    ifdefRE = re.compile(r"#ifdef\s+(?P<defword>[a-zA-Z_][a-zA-Z_0-9]*).*\n", re.IGNORECASE)
    endifRE = re.compile(r"#endif.*\n", re.IGNORECASE)
    elseRE = re.compile(r"#else.*\n", re.IGNORECASE)
    ifndefRE = re.compile(r"#ifndef\s+(?P<defword>[a-zA-Z_][a-zA-Z_0-9]*).*\n", re.IGNORECASE)

    blanklineRE = re.compile(r"\s*\n")

    endofcleanedupblock = False
    defstack = []
    cleanupstack = [False]  # we start with no cleanup as default

    for line in fin.readlines():

        if blanklineRE.match(line) and endofcleanedupblock:
            # skip all blanks lines after cleaned up region
            continue

        endofcleanedupblock = False

        match_if     = ifRE.search(line)
        match_ifdef  = ifdefRE.search(line)
        mathc_else   = elseRE.search(line)
        match_ifndef = ifndefRE.search(line)
        match_endif  = endifRE.search(line)

        if match_if:

            defstack.append("I_DO_NOT_THINK_ANYONE_WOULD_USE_THIS_AS_WORD")

            if cleanupstack[-1]:
                # We are inside "CLEANUP", so we always clean up this block
                cleanupstack.append(True)
            else:
                cleanupstack.append(False)
                fou.write(line)

        elif match_ifdef:

            dw = match_ifdef.group('defword')
            defstack.append(dw)

            if cleanupstack[-1]:
                # We are inside "CLEANUP", so we always clean up this block
                cleanupstack.append(True)
            else:
                if dw in cleanWords:
                    cleanupstack.append(True)
                else:
                    cleanupstack.append(False)
                    fou.write(line)

        elif mathc_else:

            status = cleanupstack.pop()
            dw = defstack.pop()

            if cleanupstack[-1]:
                # This whole ifdef ... else ... endif is inside "CLEANUP"
                cleanupstack.append(True)
            else:
                if status:
                    # "ifdef" block is CLEANUP, so "else" block is KEEP
                    cleanupstack.append(False)
                else: # There are two possibilies: 
                    if dw in cleanWords:
                        cleanupstack.append(True)
                    else:
                        # we keep this "ELSE" line and the following block 
                        # because the def word is in in our list for cleanup
                        cleanupstack.append(False)
                        fou.write(line)  

            defstack.append(dw)  # do we need this????

        elif match_endif:

            status = cleanupstack.pop()
            dw = defstack.pop()

            if cleanupstack[-1]:
                # This whole ifdef ... else ... endif mraked by this END is inside "CLEANUP"
                pass
            else:
                if dw in cleanWords:
                    endofcleanedupblock = True
                else:
                    fou.write(line)

        elif match_ifndef:

            dw = match_ifndef.group('defword')
            defstack.append(dw)

            if cleanupstack[-1]:
                # We are inside "CLEANUP", so we always clean up this block
                cleanupstack.append(True)
            else:
                if dw in cleanWords:
                    cleanupstack.append(False)
                else:  
                    # Its word is not on out list
                    cleanupstack.append(False)
                    fou.write(line)

        else:
            
            if not cleanupstack[-1]:
                fou.write(line)

    fin.close()
    fou.close()


if __name__== "__main__":

    usage = """
    ppCleanup -c CLEANWORD_FILE -o outputfile inputfile
    """

    if len(argv) == 1:
        print usage
        sys.exit(2)

    try:
        opts, args = getopt.getopt(argv[1:], "c:o:", [])
    
    except getopt.GetoptError:
        print "invalid calling sequence"
        print usage
        sys.exit(2)

    #defauls
    cwfile = ""
    oufile = ""
    infile = ""

    for o,a in opts:
        if o == "-c":
            cwfile = a
        if o == "-o":
            oufile = a

    try:
        infile = args[0]
    except IndexError:
        print "ERROR: an input file to be cleaned up was not specified"
        print usage
        sys.exit(2)

    if cwfile == "":
        print "ERROR: must provide -c CLEANWORD_FILE"
        print usage
        sys.exit(2)

    if oufile == "":
        print "ERROR: must provide -o outputfile"
        print usage
        sys.exit(2)

    if infile == "":
        print "ERROR: must provide an input file"
        print usage
        sys.exit(2)

    ppCleanup(cwfile, infile, oufile)
