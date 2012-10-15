#!/usr/bin/env python

import re
import sys
import itertools

# these are the macros we search for
ifRE = re.compile(r"#if.*\n", re.IGNORECASE)
endifRE = re.compile(r"#endif.*\n", re.IGNORECASE)
elseRE = re.compile(r"#else.*\n", re.IGNORECASE)
ifndefRE = re.compile(r"#ifndef.*\n", re.IGNORECASE)

# some globals that we will use 
allMatches = {}
sortedMatchKeys = []
cleanWords = []
fileContents=""

######################################################################
# Exceptions
######################################################################
class CWError(Exception):
    """
    Some error raising when we donn't specify cleanWords.  We don't want to 
    just print this out, but rather throw an error because the output of the 
    ppCleanup script could be dumped to stdout.
    """
    def __init__(self,msg,parser):
        self.msg = msg
        self.parser = parser
    def __str__(self):
        msg = "\n"*2+"-"*30
        msg += "\n\n%s" %  self.msg
        msg += "\n"*2+"-"*30+"\n"*2
        if sys.version_info[1] < 7: # less than version 2.7, using optparse
            msg += self.parser.format_help()
        else: # using argparse
            msg += self.parser.print_help()
        return msg

class MacroHandleError(Exception):
    """
    Raise this when a macro isn't properly handled.
    """
    def __init__(self,msg):
        self.msg = msg
    def __str__(self):
        return repr(self.msg)


######################################################################
# Helper routines
######################################################################
def cleanupSlices(sList):
    """
the input slices are the parts of the fileContents that we want to get
rid of, so we swap them to bits that we want to keep on output.
"""

    # invert the slice info; we pad the slices with 1 in index space to remove
    # extra lines
    keep = []
    start = 0
    for thisSlice in sList:
        keep.append(slice(start,thisSlice[0]))
        start = thisSlice[1]
    # end bit
    keep.append(slice(start,-1))
    return keep


def parseOptions(argList):
    # for python version 2.7.X
    try:
        import argparse 

        parser = argparse.ArgumentParser(description=
                                         "Remove specified C Preprocessor " +
                                         "macros and associated code from a " +
                                         "source file.")
        parser.add_argument('-c','--cleanword_file',default=None,
                            help="file containing the cleanWords; you must "+
                            "specify at least one cleanWord on the command "+
                            "line or in the cleanWord file")
        parser.add_argument('fileToClean',
                            help="file to be cleaned")
        parser.add_argument('--cw', nargs='*', default=None,
                            help="cleanup words specified on the command " +
                            "line; you must specify at least one cleanWord "+
                            "on the command line or in the cleanWord file")
        parser.add_argument('-o','--output', default=None,
                            help="optionally specify the output file; " +
                            "default behaviour is stdout")
        args = parser.parse_args(argList)

        return args, args.fileToClean, parser

    except ImportError:
    # for older python versions
        import optparse

        usage = "Usage: %prog [options] fileToClean"
        
        parser = optparse.OptionParser(usage=usage)
        parser.add_option('-c','--cleanword_file',default=None,
                          dest="cleanword_file",
                          help="file containing the cleanWords; you must "+
                          "specify at least one cleanWord on the command "+
                          "line or in the cleanWord file")
        parser.add_option('--cw', nargs='*', default=None,
                          dest="cw",
                          help="cleanup words specified on the command " +
                          "line; you must specify at least one cleanWord "+
                          "on the command line or in the cleanWord file")
        parser.add_option('-o','--output', default=None,
                          dest="output",
                          help="optionally specify the output file; " +
                          "default behaviour is stdout")
        (options,args) = parser.parse_args(argList)
        
        return options, args[0], parser


def insideFalseBlock(mindex):
    # this is a block that we want to keep, but clean the edges
    
    mkey = sortedMatchKeys[mindex]; match = allMatches[mkey]

    # output
    outSlices = []
    outSkip = 0
    
    # locals
    lSlices = []
    lSkip = 0

    # slice this line
    outSlices.append((match.start(),match.end()))

    openBlocks = 1; skip = 0

    for key in sortedMatchKeys[mindex+1:]:
        match = allMatches[key]; pattern = match.re.pattern

        lSkip += 1

        if skip > 0:
            skip -= 1
            continue
        

        # clean the ends if needed
        if "#endif" in pattern:
            openBlocks -= 1
            
            if openBlocks is 0:
              outSlices.append((match.start(),match.end()))
              break
            continue
        # if this is a CW, then we swap to a trueblock
        elif any([cw in match.group() for cw in cleanWords]):
            lSlices,skip = insideTrueBlock(mindex+lSkip)
            for lSlice in lSlices:
                outSlices.append(lSlice)
            continue
        # just count open blocks
        elif "#if" in pattern:
            openBlocks += 1
            continue

        # this is equivalent to a trueblock
        elif "#else" in pattern:
            lSlices,skip = insideTrueBlock(mindex+lSkip)
            for lSlice in lSlices:
                outSlices.append(lSlice)
            continue

        # shouldn't get here
        raise MacroHandleError("ERROR in insideFalseBlock: %s "
                               + match.group() + 
                               " macro not handled properly.")

    outSkip = lSkip

    return outSlices,outSkip

def insideTrueBlock(mindex):
    # this is a "true" block; one that we want to eliminate

    mkey = sortedMatchKeys[mindex]; match = allMatches[mkey]

    # output
    outSlices = []
    outSkip = 0

    # locals
    lSlices = []
    lSkip = 0

    start = match.start()
    if fileContents[start-1] == "\n" and fileContents[start-2] == "\n": 
        start -= 1

    # start parsing subsequent matches
    openBlocks = 1
    for key in sortedMatchKeys[mindex+1:]:
        match = allMatches[key]; pattern = match.re.pattern
        lSkip += 1

        # just count the open blocks
        if "#if" in pattern:
            openBlocks += 1
            continue

        # we've hit a false block; this will grab the endif, so we should break
        elif "#else" in pattern:
            outSlices.append((start,match.start()))
            lSlices,skip = insideFalseBlock(mindex+lSkip)
            for slice in lSlices:
                outSlices.append(slice)
            lSkip += skip
            break

        # an endif; check to make sure it is the last
        elif "#endif" in pattern:
            openBlocks -= 1
            if openBlocks is 0:
                outSlices.append((start,match.end()))
                break
            continue

        # shouldn't get here
        raise MacroHandleError("ERROR in insideTrueBlock: %s " % match.group()
                               + " macro not handled properly.")

    outSkip = lSkip
    return outSlices,outSkip

        
def insideCWBlock(mindex):
    # a CW has been found in this match; see if it is a true or false block

    mkey = sortedMatchKeys[mindex]; match = allMatches[mkey]
    pattern = match.re.pattern 

    # output
    outSlices = []
    outSkip = 0

    if "#ifndef" in pattern:
        # false block
        outSlices,outSkip = insideFalseBlock(mindex)

    elif "#if" in pattern:
        # true block
        outSlices,outSkip = insideTrueBlock(mindex)
    else:
        # shouldn't get here
        raise MacroHandleError("ERROR in insideCWBlock: %s " % match.group()
                               + "macro not handled properly.")

    return outSlices, outSkip

######################################################################
# Main routine
######################################################################
def ppCleanup(cleanWordsFile,cmdLineCW,fileToClean,outputFile,parser):
    """
This script parses a text file for the C PreProcessor #if flags containing
the words in cleanWordsFile and removes them from the text file.  The
usecase is to remove additional complexity in software before
distributing it to the masses.
"""
    global allMatches, sortedMatchKeys, cleanWords, fileContents

    # get our cleanWords from the cleanWordsFile if it exists
    if cleanWordsFile:
        try:
            fh = open(cleanWordsFile)
            cleanWords=fh.read().split()
            fh.close()
        except IOError:
            print "Couldn't open and read the file containing clean words: %s"\
                % cleanWordsFile
            sys.exit()

    # append any command line-specified cleanWords, if they exist
    if cmdLineCW:
        for CW in cmdLineCW:
            cleanWords.append(CW)

    if len(cleanWords)==0:
        raise CWError("No cleanWords specified!",parser)

    # open and read our fileToClean
    try:
        fh = open(fileToClean)
        fileContents=fh.read()
        fh.close()
    except IOError:
        print "Couldn't open and read the file to clean: %s" % fileToClean
        sys.exit()

    # find all the #if's, #endif's, #else's, and #ifndef's
    ifMatches = ifRE.finditer(fileContents)
    endifMatches = endifRE.finditer(fileContents)
    elseMatches  = elseRE.finditer(fileContents)
    ifndefMatches = ifndefRE.finditer(fileContents)

    # build a dictionary of the matches
    # key: index in fileContents of start of the match
    # value: match object
    # we chain the #ifndef matches last so they overwrite any #if.* matches
    # that also would have been caught by the re.finditer
    for match in itertools.chain(ifMatches,
                                 endifMatches,
                                 elseMatches,
                                 ifndefMatches):
        allMatches[match.start()] = match

    # get the keys in sorted order
    sortedMatchKeys = sorted([int(key) for key in allMatches.keys()])

    # store the slices we want to get rid of as tuples
    # we'll later use this as a mask for printing, using actual slice objects
    omitSlices=[]

    skip = 0

    # process each match sequentially    
    for key in sortedMatchKeys:
        # if this match has already been processed by a helper routine, then
        # skip it
        if skip > 0:
            skip -=1
            continue

        # otherwise, see if it needs processed and do it
        match = allMatches[key]; pattern = match.re.pattern
        slices = []

        if any([cw in match.group() for cw in cleanWords]):
            slices,skip = insideCWBlock(sortedMatchKeys.index(key))

        # append any slices from insideCWBlock to our omitSlices
        for thisSlice in slices:
            omitSlices.append(thisSlice)

    # swap the slices to their inverse for what we want to keep
    keepSlices = cleanupSlices(omitSlices)

    outData = ""
    for thisSlice in keepSlices:
        outData += fileContents[thisSlice]

    # dump it
    if outputFile:
        outData += fileContents[-1]
        try:
            fh = open(outputFile,'w')
            fh.write(outData)
            fh.close()
        except IOError:
            print "Couldn't open and/or write to output file %s" % outputFile
            sys.exit()
    else:
        print outData



if  __name__ == "__main__":
    # parse the command line arguments
    args, fileToClean, parser = parseOptions(sys.argv[1:])

    if not args.cleanword_file and not args.cw:
        print "You must specify either a cleanWord file (via -c) or"
        print "cleanWords from the command line (via -cw)"
        parser.print_help()

    cleanWordFile = args.cleanword_file
    cmdLineCW = args.cw
    output = args.output

    # let's do this
    ppCleanup(cleanWordFile,cmdLineCW,fileToClean,output,parser)
