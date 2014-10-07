#!/usr/bin/env python

"""
A simple regression test framework for a BoxLib-based code

There are several major sections to this source: the runtime parameter
routines, the test suite routines, and the report generation routines.  They
are separated as such in this file.

This test framework understands source based out of the F_Src and C_Src BoxLib
frameworks.
"""

import os
import shutil
import sys
import getopt
import datetime
import time
import string
import tarfile
import subprocess
import smtplib
import email
import getpass
import socket
import time


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# T E S T   C L A S S E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class testObj:

    def __init__ (self, name):

        self.name = name

        self.buildDir = ""
        self.useExtraBuildDir = 0

        self.testSrcTree = ""

        self.inputFile = ""
        self.probinFile = ""
        self.auxFiles = []
        self.linkFiles = []

        self.dim = -1

        self.needsHelmEOS = 0

        self.restartTest = 0
        self.restartFileNum = -1

        self.compileTest = 0

        self.selfTest = 0
        self.stSuccessString = ""

        self.debug = 0

        self.useMPI = 0
        self.numprocs = -1

        self.useOMP = 0
        self.numthreads = -1

        self.doVis = 0
        self.visVar = ""

        self.analysisRoutine = ""
        self.analysisMainArgs = ""
        self.analysisOutputImage = ""

        self.outputFile = ""
        self.compareFile = ""

        self.diffDir = ""
        self.diffOpts = ""

        self.addToCompileString = ""

        self.reClean = 0    # set automatically, not by users

        self.wallTime = 0   # set automatically, not by users


    def __cmp__(self, other):
        return cmp(self.value(), other.value())

    def value(self):
        return self.name


class suiteObj:

    def __init__ (self):

        self.suiteName = "testDefault"

        self.sourceTree = ""
        self.boxLibDir = ""
        self.sourceDir = ""
        self.testTopDir = ""
        self.webTopDir = ""
        self.compareToolDir = ""
        self.helmeosDir = ""

        self.useExtSrc = 0     # set automatically -- not by users
        self.extSrcDir = ""
        self.extSrcCompString = ""

        self.boxLibGitBranch = "master"
        self.sourceGitBranch = "master"
        self.extSrcGitBranch = "master"

        # this will hold the # of extra build directories we need to worry 
        # about.  It is set automatically, not by users
        self.useExtraBuild = 0   
                            
        self.extraBuildDirs = []

        # this should be the environment variable name that should be
        # set so the builds in the extraBuildDir can see the main
        # source.  This environment variable will be set to the
        # sourceTree path and included on the make lines
        self.extraBuildDirCompString = ""

        # these are set automatically by the script -- they hold the
        # basename of the various source directories
        self.srcName = ""
        self.extSrcName = ""
        self.extraBuildNames = []

        self.MPIcommand = ""
        self.MPIhost = ""

        self.FCOMP = "gfortran"
        self.COMP = "g++"

        self.MAKE = "gmake"
        self.numMakeJobs = 1

        self.reportActiveTestsOnly = 0
        self.goUpLink = 0
        self.lenTestName = 0
        
        self.sendEmailWhenFail = 0
        self.emailFrom = ""
        self.emailTo = []
        self.emailSubject = ""
        self.emailBody = ""

        self.globalAddToExecString = ""



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# R U N T I M E   P A R A M E T E R   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import ConfigParser
import string, sys, re


#==============================================================================
# some utility functions to automagically determine what the data types are
#==============================================================================
def isInt(string):
    """ is the given string an interger? """
    try: int(string)
    except ValueError: return 0
    else: return 1


def isFloat(string):
    """ is the given string a float? """
    try: float(string)
    except ValueError: return 0
    else: return 1


def convertType(string):
    """ return an integer, float, or string from the input string """

    # check in turn whether this is an interger, float, or string
    if (isInt(string)):
        value = int(string)
    elif (isFloat(string)):
        value = float(string)
    else:
        value = string.strip()

    return value


#==============================================================================
# LoadParams
#==============================================================================
def LoadParams(file):
    """
    reads the parameter file and creates as list of test objects as well as
    the suite object
    """

    testList = []

    # check to see whether the file exists
    try: f = open(file, 'r')
    except IOError:
        fail('ERROR: parameter file does not exist: %s' % (file))
    else:
        f.close()


    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.read(file)


    # "main" is a special section containing the global suite parameters.
    mysuite = suiteObj()

    for opt in cp.options("main"):

        # get the value of the current option
        value = convertType(cp.get("main", opt))

        if (opt == "suiteName"):
            mysuite.suiteName = value

        elif (opt == "sourceTree"):
            if (not (value == "C_Src" or value == "F_Src" or value == "BoxLib")):
                fail("ERROR: invalid sourceTree")
            else:
                mysuite.sourceTree = value

        elif (opt == "sourceDir"):
            mysuite.sourceDir = checkTestDir(value)

        elif (opt == "boxLibDir"):
            mysuite.boxLibDir = checkTestDir(value)

        elif (opt == "testTopDir"):
            mysuite.testTopDir = checkTestDir(value)

        elif (opt == "webTopDir"):
            mysuite.webTopDir = os.path.normpath(value) + "/"

        elif (opt == "compareToolDir"):
            mysuite.compareToolDir = checkTestDir(value)

        elif (opt == "helmeosDir"):
            mysuite.helmeosDir = checkTestDir(value)

        elif (opt == "extSrcDir"):
            mysuite.extSrcDir = checkTestDir(value)
            mysuite.useExtSrc = 1

        elif (opt == "extSrcCompString"):
            mysuite.extSrcCompString = value

        elif (opt == "boxLibGitBranch"):
            mysuite.boxLibGitBranch = value

        elif (opt == "sourceGitBranch"):
            mysuite.sourceGitBranch = value

        elif (opt == "extSrcGitBranch"):
            mysuite.extSrcGitBranch = value

        # we will keep the extra build directories in a list -- we want them
        # in the correct order in the list so we can simply index later
        elif (opt == "extraBuildDir"):
            if len(mysuite.extraBuildDirs) == 0:
                mysuite.extraBuildDirs.append(checkTestDir(value))
            else:
                mysuite.extraBuildDirs.insert(0,checkTestDir(value))

        elif (opt == "extraBuildDir2"):
            mysuite.extraBuildDirs.append(checkTestDir(value))

        elif (opt == "extraBuildDirCompString"):
            mysuite.extraBuildDirCompString = value

        elif (opt == "MPIcommand"):
            mysuite.MPIcommand = value

        elif (opt == "MPIhost"):
            mysuite.MPIhost = value
            
        elif (opt == "FCOMP"):
            mysuite.FCOMP = value

        elif (opt == "COMP"):
            mysuite.COMP = value

        elif (opt == "MAKE"):
            mysuite.MAKE = value

        elif (opt == "numMakeJobs"):
            mysuite.numMakeJobs = value

        elif (opt == "reportActiveTestsOnly"):
            mysuite.reportActiveTestsOnly = value

        elif (opt == "goUpLink"):
            mysuite.goUpLink = value

        elif (opt == "sendEmailWhenFail"):
            mysuite.sendEmailWhenFail = value

        elif (opt == "emailFrom"):
            mysuite.emailFrom = value

        elif (opt == "emailTo"):
            mysuite.emailTo = value.split(",")

        elif (opt == "emailSubject"):
            mysuite.emailSubject = value

        elif (opt == "emailBody"):
            mysuite.emailBody = value

        elif (opt == "globalAddToExecString"):
            mysuite.globalAddToExecString = value

        else:
            warning("WARNING: suite parameter %s not valid" % (opt))


    mysuite.useExtraBuild = len(mysuite.extraBuildDirs)

    mysuite.srcName = os.path.basename(os.path.normpath(mysuite.sourceDir))

    # if there is an extra source directory, update the additional
    # compilation string
    if mysuite.useExtSrc:
        mysuite.extSrcName = os.path.basename(os.path.normpath(mysuite.extSrcDir))

        if mysuite.extSrcCompString != "":
            mysuite.extSrcCompString += "="+mysuite.extSrcDir

    # if there is an extra build directory (some problems defined there),
    # then update the additional compilation string
    if mysuite.useExtraBuild > 0:
        n = 0
        while n < len(mysuite.extraBuildDirs):
            mysuite.extraBuildNames.append(os.path.basename(os.path.normpath(mysuite.extraBuildDirs[n])))
            n += 1

        # since we are building in the extraBuildDir, we need to 
        # take make where the sourceDir is
        if mysuite.extraBuildDirCompString != "":
            mysuite.extraBuildDirCompString += "="+mysuite.sourceDir


    # BoxLib-only tests don't have a sourceDir
    if (mysuite.sourceTree == "BoxLib"):
        mysuite.sourceDir = mysuite.boxLibDir


    # checks
    if (mysuite.sendEmailWhenFail):
        if (mysuite.emailTo == [] or mysuite.emailBody == ""):
            fail("ERROR: when sendEmailWhenFail = 1, you must specify emailTo and emailBody\n")
            
        if (mysuite.emailFrom == ""):
            mysuite.emailFrom = '@'.join((getpass.getuser(), socket.getfqdn()))

        if (mysuite.emailSubject == ""):
            mysuite.emailSubject = mysuite.suiteName+" Regression Test Failed"

    if (mysuite.sourceTree == "" or mysuite.boxLibDir == "" or
        mysuite.sourceDir == "" or mysuite.testTopDir == "" or
        mysuite.compareToolDir == ""):
        fail("ERROR: required suite-wide directory not specified\n" + \
                 "(sourceTree, boxLibDir, sourceDir, testTopDir, compareToolDir)")


    # if no webTopDir was specified, use the default.  In either case, make
    # sure that the web directory is valid
    if (mysuite.webTopDir == ""):
        mysuite.webTopDir = "%s/%s-web/" % (mysuite.testTopDir, mysuite.suiteName)

    if (not (os.path.isdir(mysuite.webTopDir)) ):
        try: os.mkdir(mysuite.webTopDir)
        except: fail("ERROR: unable to create the web directory: %s\n" % 
                     (mysuite.webTopDir))


    # all other sections are tests
    print "\n"
    bold("finding tests and checking parameters...")

    for sec in cp.sections():

        if (sec == "main"): continue

        print "  %s" % (sec)

        
        # lenTestName is the maximum test name length -- used for HTML
        # formatting
        mysuite.lenTestName = max(mysuite.lenTestName, len(sec))


        # create the test object for this test
        mytest = testObj(sec)


        invalid = 0

        # set the test object data by looking at all the options in
        # the current section of the parameter file
        for opt in cp.options(sec):

            # get the value of the current option
            value = convertType(cp.get(sec, opt))

            if (opt == "buildDir"):
                mytest.buildDir = value

            elif (opt == "useExtraBuildDir"):
                mytest.useExtraBuildDir = value

            elif (opt == "testSrcTree"):
                mytest.testSrcTree = value

            elif (opt == "inputFile"):
                mytest.inputFile = value

            elif (opt == "probinFile"):
                mytest.probinFile = value
                
            elif (opt == "aux1File" or opt == "aux2File" or opt == "aux3File"):
                mytest.auxFiles.append(value)

            elif (opt == "link1File" or opt == "link2File" or opt == "link3File"):
                mytest.linkFiles.append(value)
                
            elif (opt == "dim"):
                mytest.dim = value
                
            elif (opt == "needsHelmEOS"):
                mytest.needsHelmEOS = value
                
            elif (opt == "restartTest"):
                mytest.restartTest = value

            elif (opt == "restartFileNum"):
                mytest.restartFileNum = value
                
            elif (opt == "compileTest"):
                mytest.compileTest = value
                
            elif (opt == "selfTest"):
                mytest.selfTest = value

            elif (opt == "stSuccessString"):
                mytest.stSuccessString = value
                
            elif (opt == "debug"):
                mytest.debug = value

            elif (opt == "useMPI"):
                mytest.useMPI = value
                
            elif (opt == "numprocs"):
                mytest.numprocs = value

            elif (opt == "useOMP"):
                mytest.useOMP = value
                
            elif (opt == "numthreads"):
                mytest.numthreads = value

            elif (opt == "doVis"):
                mytest.doVis = value

            elif (opt == "visVar"):
                mytest.visVar = value

            elif (opt == "analysisRoutine"):
                mytest.analysisRoutine = value

            elif (opt == "analysisMainArgs"):
                mytest.analysisMainArgs = value

            elif (opt == "analysisOutputImage"):
                mytest.analysisOutputImage = value

            elif (opt == "outputFile"):
                mytest.outputFile = value

            elif (opt == "compareFile"):
                mytest.compareFile = value

            elif (opt == "diffDir"):
                mytest.diffDir = value

            elif (opt == "diffOpts"):
                mytest.diffOpts = value

            elif (opt == "addToCompileString"):
                mytest.addToCompileString = value

            else:
                warning("   WARNING: unrecognized parameter %s for test %s" % (opt, sec))


        # make sure that the build directory actually exists
        if (mytest.useExtraBuildDir > 0):
            bDir = mysuite.extraBuildDirs[mytest.useExtraBuildDir-1] + mytest.buildDir
        else:
            bDir = mysuite.sourceDir + mytest.buildDir

        if (not os.path.isdir(bDir)):
            warning("   WARNING: invalid build directory: %s" % (bDir))
            invalid = 1



        # make sure all the require parameters are present
        if (mytest.compileTest):
            if (mytest.buildDir == ""):
                warning("   WARNING: mandatory runtime parameters for test %s not set" % (sec))
                invalid = 1

        else:
            if (mytest.buildDir == "" or mytest.inputFile == "" or
                (mysuite.sourceTree == "C_Src" and mytest.probinFile == "") or 
                mytest.dim == -1):
                warning("   WARNING: mandatory runtime parameters for test %s not set" % (sec))
                warning("            buildDir = %s" % (mytest.buildDir))
                warning("            inputFile = %s" % (mytest.inputFile))
                if (mysuite.sourceTree == "C_Src"):
                    warning("            probinFile = %s" % (mytest.probinFile))
                warning("            dim = %s" % (mytest.dim))

                invalid = 1


        # check the optional parameters
        if (mytest.restartTest and mytest.restartFileNum == -1):
            warning("   WARNING: test %s is a restart test but restartFileNum not set" % (sec))
            invalid = 1

        if (mytest.selfTest and mytest.stSuccessString == ""):
            warning("   WARNING: test %s is a self-test but stSuccessString not set" % (sec))
            invalid = 1

        if (mytest.useMPI and mytest.numprocs == -1):
            warning("   WARNING: test %s is an MPI parallel test but numprocs not set" % (sec))
            invalid = 1

        if (mytest.useOMP and mytest.numthreads == -1):
            warning("   WARNING: test %s is an OpenMP parallel test but numthreads not set" % (sec))
            invalid = 1

        if (mytest.doVis and mytest.visVar == ""):
            warning("   WARNING: test %s requested visualization visVar not set" % (sec))
            invalid = 1

        if (mysuite.sourceTree == "BoxLib" and mytest.testSrcTree == ""):
            warning("   WARNING: test %s is a BoxLib test but testSrcTree not set" % (sec))
            invalid = 1
        

        # add the current test object to the master list
        if (not invalid):
            testList.append(mytest)

        else:
            warning("   WARNING: test %s will be skipped" % (sec))

    
    # final checks
    
    # if any runs are parallel, make sure that the suite-wide
    # MPIcommand is defined
    anyMPI = 0
    for test in testList:
        if (test.useMPI):
            anyMPI = 1
            break

    if (anyMPI and mysuite.MPIcommand == ""):
        fail("ERROR: one or more tests are parallel, but MPIcommand is not defined")

    # if any runs use helmeos, make sure that the suite-wide
    # helmeosDir is defined
    anyhelmeos = 0
    for test in testList:
        if (test.needsHelmEOS == 1):
            anyhelmeos = 1
            break

    if (anyhelmeos and mysuite.helmeosDir == ""):
        fail("ERROR: one or more tests use helmeos, but helmeosDir is not defined")

    testList.sort()

    return mysuite, testList



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# O U T P U T   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# inspiration from 
# http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python 
# which in-turn cites the blender build scripts              
class termColors:
    WARNING = '\033[33m'
    SUCCESS = '\033[32m'
    FAIL = '\033[31m'
    BOLD = '\033[1m'
    ENDC = '\033[0m'


def fail(str):
    new_str = termColors.FAIL + str + termColors.ENDC
    print new_str
    sys.exit()

def testfail(str):
    new_str = termColors.FAIL + str + termColors.ENDC
    print new_str

def warning(str):
    new_str = termColors.WARNING + str + termColors.ENDC
    print new_str

def success(str):
    new_str = termColors.SUCCESS + str + termColors.ENDC
    print new_str

def bold(str):
    new_str = termColors.BOLD + str + termColors.ENDC
    print new_str



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# S Y S T E M   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def systemCall(string):    
    status = os.system('bash -c "' + string + '"')
    return status


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# T E S T   S U I T E   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#==============================================================================
# findBuildDirs
#==============================================================================
def findBuildDirs(testList):
    """ given the list of test objects, find the set of UNIQUE build 
        directories.  Note if we have the useExtraBuildDir flag set """
    
    buildDirs = []
    reClean = []

    for obj in testList:
        
        # be sneaky here.  We'll add a "+" to any of the tests that
        # are built in the first extraBuildDir and a "@" for an tests
        # that are build in the second extraBuild dir, instead of the
        # sourceDir.
        if obj.useExtraBuildDir == 1:
            prefix = "+"
        elif obj.useExtraBuildDir == 2:
            prefix = "@"
        else:
            prefix = ""


        # first find the list of unique build directories
        if (buildDirs.count(prefix + obj.buildDir) == 0):
            buildDirs.append(prefix + obj.buildDir)

        # sometimes a problem will specify an extra argument to the
        # compile line.  If this is the case, then we want to re-make
        # "clean" for ALL tests that use this build directory, just to
        # make sure that any unique build commands are seen.
        if (not obj.addToCompileString == ""):
            reClean.append(obj.buildDir)


    for bdir in reClean:
        for obj in testList:
            if (obj.buildDir == bdir):
                obj.reClean = 1

    return buildDirs


#==============================================================================
# getLastPlotfile
#==============================================================================
def getLastPlotfile(outputDir, test):
    """ given an output directory and the test name, find the last
        plotfile written """
        
    plotNum = -1
        
    # start by finding the last plotfile
    for file in os.listdir(outputDir):
       if (os.path.isdir(file) and file.startswith("%s_plt" % (test.name))):
           key = "_plt"
           index = string.rfind(file, key)
           plotNum = max(int(file[index+len(key):]), plotNum)

    if (plotNum == -1):
       warning("WARNING: test did not produce any output")
       compareFile = ""
    else:
       compareFile = "%s_plt%5.5d" % (test.name, plotNum)

    return compareFile


#==============================================================================
# getRecentFileName
#==============================================================================
def getRecentFileName(dir,base,extension):
    """ for F_Src builds, given the base and extension, find the
        most recent corresponding file"""
        
    # finding all files that are of type base.*.exe and store the
    # name of only the most recently created one
    ctime = -1
    executableFile = ""
    for file in os.listdir(dir):
       if (os.path.isfile(file) and
           file.startswith(base) and file.endswith(extension)):

          fileInfo = os.stat(file)
          fileCreationTime = fileInfo.st_ctime
          if (fileCreationTime > ctime):
             ctime = fileCreationTime
             executableFile = file

    return executableFile


#==============================================================================
# getTestCompDir
#==============================================================================
def checkTestDir(dirName):
   """ given a string representing a directory, check if it points to
       a valid directory.  If so, return the directory name """

   # make sure we end in a "/"
   if (not (string.rfind(dirName, "/") == len(dirName)-1)):
       dirName = dirName + "/"
           
   if (not os.path.isdir(dirName)):
       fail("ERROR: %s is not a valid directory" % (dirName))

   return dirName


#==============================================================================
# doGITUpdate
#==============================================================================
def doGITUpdate(topDir, root, outDir, gitbranch, githash):
   """ do a git update of the repository in topDir.  root is the name
       of the directory (used for labeling).  outDir is the full path
       to the directory where we will store the git output.  If githash
       is not empty, then we will check out that version instead of
       git-pulling."""

   os.chdir(topDir)


   # find out current branch so that we can go back later if we need.
   prog = ["git", "rev-parse", "--abbrev-ref", "HEAD"]
   p0 = subprocess.Popen(prog, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
   stdout0, stderr0 = p0.communicate()
   currentBranch = stdout0.rstrip('\n')
   p0.stdout.close()

   if currentBranch != gitbranch:
       print "\n"
       bold("git checkout %s in %s" % (gitbranch, topDir))
       prog = ["git", "checkout", gitbranch]
       p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
       stdout, stderr = p.communicate()
       p.stdout.close()
       p.stdin.close()

   if githash == "":

       print "\n"
       bold("'git pull' in %s" % (topDir))

       # we need to be tricky here to make sure that the stdin is
       # presented to the user to get the password.  Therefore, we use
       # the subprocess class instead of os.system
       prog = ["git", "pull"]
       p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
       stdout, stderr = p.communicate()
       p.stdout.close()
       p.stdin.close()

   else:

       prog = ["git", "checkout", githash]
       p = subprocess.Popen(prog, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
       stdout, stderr = p.communicate()
       p.stdout.close()


   try:
       cf = open("git.%s.out" % (root), 'w')
   
   except IOError:
       fail("  ERROR: unable to open file for writing")

   else:

       for line in stdout:
           cf.write(line)

       cf.close()


   if (stdout == ""):
       fail("  ERROR: git update was unsuccessful")
       
        
   shutil.copy("git.%s.out" % (root),  outDir)

   return currentBranch


#==============================================================================
def saveGITHEAD(topDir, root, outDir):

   os.chdir(topDir)

   print "\n"
   bold("saving git HEAD for %s/" % (root))

   systemCall("git rev-parse HEAD >& git.%s.HEAD" % (root) )

   shutil.copy("git.%s.HEAD" % (root),  outDir)


#==============================================================================
# doGITback
#==============================================================================
def doGITback(topDir, root, gitbranch):
   """ do a git checkout of gitbranch in topDir.  root is the name
       of the directory (used for labeling). """

   os.chdir(topDir)

   print "\n"
   bold("git checkout %s in %s" % (gitbranch, topDir))

   prog = ["git", "checkout", gitbranch]
   p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT)
   stdout, stderr = p.communicate()
   p.stdout.close()
   p.stdin.close()

   try:
       cf = open("git.%s.out" % (root), 'w')   
   except IOError:
       fail("  ERROR: unable to open file for writing")
   else:
       for line in stdout:
           cf.write(line)
       cf.close()

   if (stdout == ""):
       fail("  ERROR: git checkout was unsuccessful")


#==============================================================================
# makeGITChangeLog
#==============================================================================
def makeGITChangeLog(gitDir, root, outDir):
    """ generate a ChangeLog git repository named root.  outDir is the
        full path to the directory where we will store the git output"""

    os.chdir(gitDir)


    print "\n"
    bold("generating ChangeLog for %s/" % (root))
    
    systemCall("git log --name-only >& ChangeLog.%s" % (root) )
    
    shutil.copy("ChangeLog.%s" % (root), outDir)


#==============================================================================
# allAreCompileTests
#==============================================================================
def allAreCompileTests(testList):
    """ return 1 if all the tests in the list are compile tests """

    allCompile = 1
    for test in testList:
        if (not test.compileTest):
            allCompile = 0
            break
            
    return allCompile


#==============================================================================
# getLastRun
#==============================================================================
def getLastRun(suite):
    """ return the name of the directory corresponding to the previous
        run of the test suite """

    outdir = suite.testTopDir + suite.suiteName + "-tests/"
    
    dirs = []
    for dir in os.listdir(outdir):
        # this will work through 2099
        if os.path.isdir(outdir + dir) and dir.startswith("20"):
            dirs.append(dir)

    dirs.sort()

    return dirs[-1]

#==============================================================================
# getTestFailures
#==============================================================================
def getTestFailures(suite, testDir):
    """ look at the test run in testDir and return the list of tests that 
        failed """

    cwd = os.getcwd()

    outdir = suite.testTopDir + suite.suiteName + "-tests/"

    os.chdir(outdir + testDir)
            
    failed = []

    for test in os.listdir("."):
        if not os.path.isdir(test): continue

        # the status files are in the web dir
        statusFile = suite.webTopDir + testDir + "/%s.status" % (test)
        sf = open(statusFile, "r")
        lines = sf.readlines()
        for line in lines:
            if string.find(line, "FAILED") >= 0:
                failed.append(test)

    os.chdir(cwd)

    return failed


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# test
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def testSuite(argv):

    usage = """
    ./testnew.py [--make_benchmarks comment,
                  --no_update  none or all or a list of codes excluded from update,
                  --single_test test
                  --tests "test1 test2 test3 ..."
                  --do_temp_run
                  --boxLibGitHash boxlibhash
                  --sourceGitHash sourcehash
                  --extSrcGitHash extsourcehash
                  --note note
                  -d dimensionality
                  --redo_failed
                  --complete_report_from_crash testdir]
        testfile.ini


    arguments:

      testfile.ini 
          This is the input file that defines the tests that the
          suite knows about.  It has the format

            [main]
            boxLibDir      = < directory to the BoxLib/ directory >
            sourceDir      = < directory to the main Source directory >
            testTopDir     = < full path to test output directory >
            webTopDir      = < full path to test web directory >
            compareToolDir = < full path to the AmrPostprocessing/F_Src directory >
            helmeosDir     = < full path to helm_table.dat, 
                               no need to set this if not using helmeos>
            extSrcDir      = < directory to an extra source directory other than
                               BoxLib and the main source, if there is one >
            extSrcCompString = < a string.  If both extSrcCompString and extSrcDir
                                 are set, they will be used by make as follows,
                                 make extSrcCompString=extSrcDir > 

            sourceTree = < C_Src, F_Src, or BoxLib -- what type is it? >

            suiteName = < descriptive name (i.e. Castro) >


            FCOMP = < name of Fortran compiler >
            COMP  = < name of C/C++ compiler >

            MAKE        = < name of make > 
            numMakeJobs = < number of make jobs >

            MPIcommand = < MPI run command, with placeholders for host,
                           # of proc, and program command.  Should look
                           something like :
                           mpiexec -host @host@ -n @nprocs@ @command@ >
                           
            MPIhost = < host for MPI job -- depends on MPI implementation >

            reportActiveTestsOnly = <If 1, inactive tests will not be include in the web page.>

            goUpLink = <If 1, add "Go UP" link at top of the web page.>

            sendEmailWhenFail = < If 1, send email when any tests fail>
            emailTo           = < list of email addresses separated by commas, 
                                  such as,
                                  foo@example.com, bar@example.com >
            emailBody         = < email body >

            [Sod-x]
            buildDir = < relative path (from sourceDir) for this problem >
            inputFile = < input file name >
            probinFile = < probin file name >
            dim = < dimensionality: 1, 2, or 3 >
            needs_helmeos = < need Helmholtz eos? 0 for no, 1 for yes >

            aux1File = < name of additional file needed by the test >
            link1File = < name of additional file needed by the test >

            restartTest = < is this a restart test? 0 for no, 1 for yes >
            restartFileNum = < # of file to restart from (if restart test) >

            useMPI = <is this a parallel (MPI) job? 0 for no, 1 for yes) >
            numprocs = < # of processors to run on (if parallel job) >

            useOMP = <is this an OpenMP job? 0 for no, 1 for yes) >
            numthreads = < # of threads to us with OpenMP (if OpenMP job) >

            debug = < 0 for normal run, 1 if we want debugging options on >

            compileTest = < 0 for normal run, 1 if we just test compilation >

            selfTest = < 0 for normal run, 1 if test self-diagnoses if it 
                         succeeded >
            stSuccessString = < string to search for in self-test output to 
                         determine success >

            doVis = < 0 for no visualization, 1 if we do visualization >
            visVar = < string of the variable to visualize >

            analysisRoutine = < name of the script to run on the output.  The
                                script is run as:

                                analysisRoutine [options] plotfile >

            analysisMainArgs = < commandline arguments to pass to the 
                                 analysisRoutine -- these should refer to
                                 options from the [main] block >

            analysisOutputImage = < result of the analysis to show on the 
                                    test results page >

            compareFile = < explicit output file to do the comparison with >

            diffDir = < directory or file to do a plain text diff on 
                       (recursive, if directory) >

            diffOpts = < options to use with the diff command for the diffDir
                        comparison >

          Here, [main] lists the parameters for the test suite as a
          whole and [Sod-x] is a single test.  There can be many more
          tests, each with their own unique name, following the format
          of [Sod-x] above.

          In [main],

            boxLibDir is the full path to the BoxLib/ directory.
            
            sourceDir should be the full path to the directory where
            the application source directories live.  For example, 
            this could point to Castro/

            testTopDir is the full path to the directory that will
            contain the test run directories, the test web
            directories, and the test benchmark directories.  It can
            be the same as sourceDir.  The necessary sub- directories
            will be created at runtime if they don't already exist.

            compareToolDir is the full path to the directory
            containing the comparison tool.  This resides in
            AmrPostprocesing/F_Src

            helmeosDir is the full path to the directory containing
            the helm_table.dat file needed by the general EOS.

            sourceTree is either C_Src for applications in the C++
            BoxLib framework or F_Src for applications in the F90
            BoxLib framework.  For internal BoxLib tests, sourceTree
            is set to BoxLib.  In this case, each test will need to
            specify what source tree to use for itself (through the
            test parameter testSrcTree).

            suiteName is a descriptive name for the test run.  The
            output directories will be named suiteName-tests/ and
            suiteName-web/.  The benchmark directory will be
            suiteName-benchmarks/.

            FCOMP is the name of the Fortran compiler -- this should be
            a name that the Makefiles of the code recognize.

            COMP is the name of the C++ compiler -- this should be a 
            name that the Makefils of the code recognize.

            MAKE is the name of the make utility -- this should be name that
            the operating system recongnizes.

            numMakeJobs is the number of make jobs to run simultaneously.

            To run jobs in Parallel, the following need to be set:

               MPIcommand is the full MPI run command, with placeholders
               for the hosts to run on, number of processors, and executable.
               For MPICH-2, this would look like:

                 MPIcommand = mpiexec -host @host@ -n @nprocs@ @command@

               MPIhost is the name of the host to run on -- you may or may
               not need this, depending on the machine you are running on.

               The number of processors is problem dependent and will be
               set in the problem specific blocks described below.

          For a test problem (e.g. [Sod-x]),  

            buildDir is the relative path (wrt sourceDir) where the
            make command will be issued for the test problem.
            
            inputFile and probinFile are the names of the input file
            and associated probin file, respectively.  Note, a
            probinFile is only required for a C_Src (not F_Src)
            sourceTree run.

            dim is the dimensionality for the problem.

            needs_helmeos is set to 1 if the Helmholtz EOS is used.
            This will ensure that the helm_table.dat file is copied
            into the run directory.
            
            aux1File (also aux2File and aux3File) is the name of 
            any additional file needed by the test.  This will be 
            COPIED into the run directory

            link1File (also link2File and link3File) is the name of
            any additional file needed by the test.  This can be
            used instead of the aux files.  The difference is that the
            link files are SYMLINKed into the run directory, instead
            of copied.

            restartTest = 1 means that this is a restart test.  Instead of
            comparing to a stored benchmark, we will run the test and then
            restart from restartFileNum and run until the end.  The last
            file from the original run and the last from the restart will
            be compared.

            useMPI is set to 1 if the job is parallel (with MPI).  In
            this case, you also must specify the number of processors,
            via numprocs.

            useOMP is set to 1 if the job uses OpenMP.  In this case,
            you also must specify the number of threads, via
            numthreads.

            debug is set to 1 if we want to turn on the debugging
            options defined in the problem's makefile.

            compileTest is set to 0 for a normal test.  Setting to 1 will
            just test compilation, and not any subsequent output.

            selfTest is set to 0 for a normal test, in which case the
            test output will be compared to the stored benchmark.  Setting 
            to 1 will determine success by searching for the string
            stSuccessString in the execution output.

            compareFile is the name of the plotfile that is output that
            should be used for the comparison.  Normally this is not
            specified and the suite uses the last plotfile output by the
            test.

            diffDir is the optional name of a directory (or single
            file) output by the test that you wish to do a plain test
            comparison on with a stored benchmark version of the
            directory.  This is just a straight diff (recursively,
            within the directory).  This is used, for example, for
            particle output from some of the codes.  diffOpts is
            a string providing options to the diff command.  For
            example: '-I "^#"' to ignore command lines in 
            Maestro diag files.

          Each test problem should get its own [testname] block
          defining the problem.  The name between the [..] will be how
          the test is referred to on the webpages.


    options:
    
       --make_benchmarks \"comment\"
          run the test suite and make the current output the new
          benchmarks for comparison.  When run in this mode, no
          comparison is done.  This is useful for the first time
          the test suite is run.

          \"comment\" describes the reason for the benchmark
          update and will be appended to the web output for
          future reference.

       --no_update ?
          None or All or a list of codes seperated by "," to be excluded
          from update. Default is None.

       --single_test mytest
          run only the test named mytest

       --tests \"test1 test2 test3\"
          run only the tests listsed

       --do_temp_run
          Temporary run without updating the web.

       --boxLibGitHash boxlibhash
          Git hash of a version of BoxLib.  If provided, this version
          will be used to run tests.

       --sourceGitHash sourcehash
          Git hash of a version of the source code.  For BoxLib tests,
          this will be ignored.

       --extSrcGitHash extsourcehash
          Git hash of a version of the source code.  For BoxLib tests,
          this will be ignored.

       --note \"note\"
          print the note on the resulting test webpages
          
       --complete_report_from_crash \"testdir\"
          if a test suite run crashed (or the terminal you were running 
          from disconnected) and the file report was not generated, this
          option will just generate the report for the test suite run and
          the overall report for all runs.  testdir is the name of the
          testdir (e.g. 20XX-XX-XX) that was running when the suite crashed.

       --redo_failed
          only run the tests that failed the last time the suite was run

    Getting started:

      To set up a test suite, it is probably easiest to write the
      testfile.ini as described above and then run the test routine
      with the --make_benchmarks option to create the benchmark
      directory.  Subsequent runs can be done as usual, and will
      compare to the newly created benchmarks.  If differences arise
      in the comparisons due to (desired) code changes, the benchmarks
      can be updated using --make_benchmarks to reflect the new
      ``correct'' solution.
      
    """

            
    #--------------------------------------------------------------------------
    # parse the commandline arguments
    #--------------------------------------------------------------------------
    if len(argv) == 1:
        print usage
        sys.exit(2)

    try:
        opts, next = getopt.getopt(argv[1:], "d:",
                                   ["make_benchmarks=",
                                    "no_update=",
                                    "single_test=",
                                    "tests=",
                                    "do_temp_run",
                                    "boxLibGitHash=",
                                    "sourceGitHash=",
                                    "extSrcGitHash=",
                                    "note=",
                                    "complete_report_from_crash=",
                                    "redo_failed"])

    except getopt.GetoptError:
        print "invalid calling sequence"
        print usage
        sys.exit(2)


    # defaults
    dimensionality = -1
    make_benchmarks = 0
    no_update = "None"
    single_test = ""
    tests = ""
    comment = ""
    do_temp_run = False
    boxLibGitHash = ""
    sourceGitHash = ""
    extSrcGitHash = ""
    note = ""
    completeReportFromCrash = ""
    redo_failed = 0

    for o, a in opts:

        if o == "-d":
            dimensionality = int(a)

        if o == "--make_benchmarks":
            make_benchmarks = 1
            comment = a

        if o == "--redo_failed":
            redo_failed = 1

        if o == "--no_update":
            no_update = a

        if o == "--single_test":
            single_test = a

        if o == "--tests":
            tests = a

        if o == "--do_temp_run":
            do_temp_run = True

        if o == "--boxLibGitHash":
            boxLibGitHash = a
            
        if o == "--sourceGitHash":
            sourceGitHash = a

        if o == "--extSrcGitHash":
            extSrcGitHash = a

        if o == "--note":
            note = a

        if o == "--complete_report_from_crash":
            completeReportFromCrash = a

    try:
        testFile = next[0]

    except IndexError:
        print "ERROR: a test file was not specified"
        print usage
        sys.exit(2)

    #--------------------------------------------------------------------------
    # read in the test information
    #--------------------------------------------------------------------------
    bold("loading " + testFile)

    suite, testList = LoadParams(testFile)

    defined_tests = testList[:]

    # if we only want to run the tests that failed previously, remove the
    # others
    if redo_failed == 1:
        last_run = getLastRun(suite)
        failed = getTestFailures(suite, last_run)

        testListold = testList[:]
        for t in testListold:
            if not t.name in failed:
                testList.remove(t)


    # if we only want to run tests of a certain dimensionality, remove
    # the others
    if dimensionality >= 1 and dimensionality <= 3:
        testListold = testList[:]
        for t in testListold:
            if not t.dim == dimensionality:
                testList.remove(t)

    activeTestList = [t.name for t in defined_tests]

    # store the full path to the testFile
    testFilePath = os.getcwd() + '/' + testFile

    if (len(testList) == 0):
        fail("No valid tests defined")


    if not completeReportFromCrash == "":

        # make sure the web directory from the crash run exists
        fullWebDir = "%s/%s/"  % (suite.webTopDir, completeReportFromCrash)
        if not os.path.isdir(fullWebDir):
            fail("Crash directory does not exist")

        # find all the tests that completed in that web directory
        tests = []
        testFile = ""
        wasBenchmarkRun = 0
        for file in os.listdir(fullWebDir):
            if os.path.isfile(file) and file.endswith(".status"):
                index = string.rfind(file, ".status")
                tests.append(file[:index])

                f = open(fullWebDir + file, "r")
                for line in f:
                    if string.find(line, "benchmarks updated") > 0:
                        wasBenchmarkRun = 1

            if os.path.isfile(file) and file.endswith(".ini"):
                testFile = file


        # create the report for this test run
        numFailed = reportThisTestRun(suite, wasBenchmarkRun, "", 
                                      "recreated report after crash of suite", 
                                      "",  0, 0, 0,
                                      tests, completeReportFromCrash, testFile, fullWebDir)


        # create the suite report
        bold("creating suite report...")
        tableHeight = min(max(suite.lenTestName, 4), 16)
        reportAllRuns(suite, activeTestList, suite.webTopDir, tableHeight=tableHeight)
        sys.exit("done")


    #--------------------------------------------------------------------------
    # figure out which git repos we will update
    #--------------------------------------------------------------------------
    no_update_low = no_update.lower()

    if no_update_low == "none":
        updateBoxLib = True
        updateSource = True
        updateExtSrc = True
        updateExtraBuild = [True]*suite.useExtraBuild

    elif no_update_low == "all":
        updateBoxLib = False
        updateSource = False
        updateExtSrc = False
        updateExtraBuild = [False]*suite.useExtraBuild

    else:
        nouplist = no_update_low.split(",")

        if "boxlib" in nouplist:
            updateBoxLib = False
        else:
            updateBoxLib = True

        if suite.srcName.lower() in nouplist:
            updateSource = False
        else:
            updateSource = True

        if suite.extSrcName.lower() in nouplist:
            updateExtSrc = False
        else:
            updateExtSrc = True

        # each extra build directory has its own update flag
        updateExtraBuild = []
        for e in suite.extraBuildNames:
            if e.lower() in nouplist:
                updateExtraBuild.append(False)
            else:
                updateExtraBuild.append(True)

    if not suite.useExtSrc:
        updateExtSrc = False

    if not suite.useExtraBuild:
        updateExtraBuild = [False]
    
    if suite.sourceTree == "BoxLib":
        updateSource = False # to avoid updating BoxLib twice.
                             # The update of BoxLib is controlled by updateBoxLib
        sourceGitHash = ""


    if boxLibGitHash:
        updateBoxLib = False 

    if sourceGitHash:
        updateSource = False 

    if extSrcGitHash:
        updateExtSrc = False 



    #--------------------------------------------------------------------------
    # if we are doing a single test, remove all other tests
    # if we specified a list of tests, check each one
    # if we did both --single_test and --tests, complain
    #--------------------------------------------------------------------------
    if (not single_test == "" and not tests == ""):
        fail("ERROR: specify tests either by --single_test or --tests, not both")

    if (not single_test == ""):
        found = 0
        for obj in testList:
            if (obj.name == single_test):
                found = 1
                newTestList = [obj]
                break

        if (not found):
            fail("ERROR: %s is not a valid test" % (single_test))
        else:
            testList = newTestList
        
    elif (not tests == ""):
        testsFind = list(set(string.split(tests)))

        newTestList = []
        for test in testsFind:
            found = 0
            for obj in testList:
                if (obj.name == test):
                    found = 1
                    newTestList.append(obj)
                    break
            
            if (not found):
                fail("ERROR: %s is not a valid test" % (test))
        
        testList = newTestList
    

    #--------------------------------------------------------------------------
    # get the name of the benchmarks directory
    #--------------------------------------------------------------------------
    allCompile = allAreCompileTests(testList)

    if (not allCompile):
        benchDir = suite.testTopDir + suite.suiteName + "-benchmarks/"
        if (not os.path.isdir(benchDir)):
            if (make_benchmarks):
                os.mkdir(benchDir)
            else:
                fail("ERROR: benchmark directory, %s, does not exist" % (benchDir))
    

    #--------------------------------------------------------------------------
    # create the output directories
    #--------------------------------------------------------------------------
    os.chdir(suite.testTopDir)

    todayDate = datetime.date.today()
    today = todayDate.__str__()


    # figure out what the current output directory should be
    maxRuns = 100      # maximum number of tests in a given day

    testDir = today + "/"

    # test output stored in a directory suiteName-tests/2007-XX-XX/
    # make sure that the suiteName-tests directory exists
    if (not os.path.isdir(suite.testTopDir + suite.suiteName + "-tests/")):
        os.mkdir(suite.testTopDir + suite.suiteName + "-tests/")
        
    fullTestDir = suite.testTopDir + suite.suiteName + "-tests/" + testDir

    if do_temp_run:
        testDir = "TEMP_RUN/"
        fullTestDir = suite.testTopDir + suite.suiteName + "-tests/" + testDir
        systemCall("rm -rf %s" % (fullTestDir))
    else:
        i = 0
        while (i < maxRuns-1 and os.path.isdir(fullTestDir)):
            i = i + 1
            testDir = today + "-" + ("%3.3d" % i) + "/"
            fullTestDir = suite.testTopDir + suite.suiteName + "-tests/" + testDir

    print "\n"
    bold("testing directory is: " + testDir)
    os.mkdir(fullTestDir)

    # make the web directory -- this is where all the output and HTML will be
    # put, so it is easy to move the entire test website to a different disk
    fullWebDir = "%s/%s/"  % (suite.webTopDir, testDir)

    if do_temp_run:
        systemCall("rm -rf %s" % (fullWebDir))
        
    os.mkdir(fullWebDir)

    # copy the test file into the web output directory
    shutil.copy(testFilePath, fullWebDir)


    #--------------------------------------------------------------------------
    # do the source updates
    #--------------------------------------------------------------------------
    now = time.localtime(time.time())
    updateTime = time.strftime("%Y-%m-%d %H:%M:%S %Z", now)

    os.chdir(suite.testTopDir)

    if updateSource or sourceGitHash:

        # main suite
        sourceGitBranch_Orig = doGITUpdate(suite.sourceDir, 
                                           suite.srcName, fullWebDir,
                                           suite.sourceGitBranch,
                                           sourceGitHash)
    
    if updateExtSrc or extSrcGitHash:

        # extra source
        if (suite.useExtSrc):
            extSrcGitBranch_Orig = doGITUpdate(suite.extSrcDir, 
                                               suite.extSrcName, fullWebDir,
                                               suite.extSrcGitBranch,
                                               extSrcGitHash)

    if any(updateExtraBuild):

        # extra build directory
        n = 0
        while n < suite.useExtraBuild:
            if updateExtraBuild[n]:
                extSrcGitBranch_Orig = doGITUpdate(suite.extraBuildDirs[n], 
                                                   suite.extraBuildNames[n], fullWebDir,
                                                   "master",
                                                   "")
            n += 1

    if updateBoxLib or boxLibGitHash:

        # BoxLib
        boxLibGitBranch_Orig = doGITUpdate(suite.boxLibDir, 
                                           "BoxLib", fullWebDir,
                                           suite.boxLibGitBranch,
                                           boxLibGitHash)

    #--------------------------------------------------------------------------
    # Save git HEADs
    #--------------------------------------------------------------------------
    saveGITHEAD(suite.boxLibDir, "BoxLib", fullWebDir)
    
    if suite.sourceTree != "BoxLib":
        saveGITHEAD(suite.sourceDir, suite.srcName, fullWebDir)

    if suite.useExtSrc:
        saveGITHEAD(suite.extSrcDir, suite.extSrcName, fullWebDir)

    n = 0
    while n < suite.useExtraBuild:
        saveGITHEAD(suite.extraBuildDirs[n], suite.extraBuildNames[n], fullWebDir)
        n += 1


    #--------------------------------------------------------------------------
    # generate the ChangeLogs
    #--------------------------------------------------------------------------
    if updateSource:

        # main suite
        makeGITChangeLog(suite.sourceDir, suite.srcName, fullWebDir)
        
    if updateExtSrc:

        # extra source
        if (suite.useExtSrc):
            makeGITChangeLog(suite.extSrcDir, suite.extSrcName, fullWebDir)


    # extra build directories            
    n = 0
    while n < suite.useExtraBuild:
        if updateExtraBuild[n]:
            makeGITChangeLog(suite.extraBuildDirs[n], suite.extraBuildNames[n], fullWebDir)
        n += 1

    if updateBoxLib:

        # BoxLib
        makeGITChangeLog(suite.boxLibDir, "BoxLib", fullWebDir)


    #--------------------------------------------------------------------------
    # build the comparison and visualization tools
    #--------------------------------------------------------------------------
    print "\n"
    bold("building the comparison tools...")

    os.chdir(suite.compareToolDir)

    compString = "%s BOXLIB_HOME=%s COMP=%s realclean >& /dev/null" % \
                   (suite.MAKE, suite.boxLibDir, suite.FCOMP)
    print "  " + compString
    systemCall(compString)

    compString = "%s -j%s BOXLIB_HOME=%s programs=fcompare NDEBUG=t MPI= COMP=%s  >& fcompare.make.out" % \
                   (suite.MAKE, suite.numMakeJobs, suite.boxLibDir, suite.FCOMP)
    print "  " + compString
    systemCall(compString)
    compareExecutable = getRecentFileName(suite.compareToolDir,"fcompare",".exe")

    shutil.copy(compareExecutable, fullTestDir + "/fcompare.exe")


    anyDoVis = {'2D':0, '3D':0}
    for test in testList:
        if test.doVis:
            if test.dim == 2:
                anyDoVis['2D'] += 1
            elif test.dim == 3:
                anyDoVis['3D'] += 1

    if anyDoVis['2D'] or anyDoVis['3D']:
        print "\n"
        bold("building the visualization tools...")

        if anyDoVis['2D']:
            compString = "%s -j%s BOXLIB_HOME=%s programs=fsnapshot2d NDEBUG=t MPI= COMP=%s  2>&1 > fsnapshot2d.make.out" % \
                         (suite.MAKE, suite.numMakeJobs, suite.boxLibDir, suite.FCOMP)
            print "  " + compString
            systemCall(compString)
            vis2dExecutable = getRecentFileName(suite.compareToolDir,"fsnapshot2d",".exe")

        if anyDoVis['3D']:
            compString = "%s -j%s BOXLIB_HOME=%s programs=fsnapshot3d NDEBUG=t MPI= COMP=%s  2>&1 > fsnapshot3d.make.out" % \
                         (suite.MAKE, suite.numMakeJobs, suite.boxLibDir, suite.FCOMP)
            print "  " + compString
            systemCall(compString)
            vis3dExecutable = getRecentFileName(suite.compareToolDir,"fsnapshot3d",".exe")
    

    #--------------------------------------------------------------------------
    # output test list
    #--------------------------------------------------------------------------
    print "\n"
    bold("running tests: ")
    for obj in testList:
        print "  %s " % obj.name
    

    #--------------------------------------------------------------------------
    # do a make clean, only once per build directory
    #--------------------------------------------------------------------------
    allBuildDirs = findBuildDirs(testList)

    print "\n"
    bold("make clean in...")

    for dir in allBuildDirs:

        # if the buildDir has a "+" at the start, that means it is 
        # in the extraBuildDir path.
        if (dir.find("+") == 0):
            inExtra = 1
            dir = dir[1:]
        elif (dir.find("@") == 0):
            inExtra = 2
            dir = dir[1:]
        else:
            inExtra = 0

        if (inExtra > 0):
            print "  %s in %s" % (dir,suite.extraBuildNames[inExtra-1])
            os.chdir(suite.extraBuildDirs[inExtra-1] + dir)
        else:
            print "  %s" % (dir)
            os.chdir(suite.sourceDir + dir)
            
        systemCall("%s BOXLIB_HOME=%s %s %s realclean >& /dev/null" % 
                   (suite.MAKE, suite.boxLibDir, 
                    suite.extSrcCompString, suite.extraBuildDirCompString))

            
    os.chdir(suite.testTopDir)
    

    #--------------------------------------------------------------------------
    # main loop over tests
    #--------------------------------------------------------------------------
    for test in testList:


        print "\n"
        bold("working on test: %s" % (test.name))

        if (test.restartTest and make_benchmarks):
            warning("  WARNING: test %s is a restart test -- " % (test.name))
            warning("           no benchmarks are stored.")
            warning("           skipping\n")
            continue
     
        if (test.compileTest and make_benchmarks):
            warning("  WARNING: test %s is a compile test -- " % (test.name))
            warning("           no benchmarks are stored.")
            warning("           skipping\n")
            continue            

        if (test.selfTest and make_benchmarks):
            warning("  WARNING: test %s is a self-test -- " % (test.name))
            warning("           no benchmarks are stored.")
            warning("           skipping\n")
            continue            



        #----------------------------------------------------------------------
        # make the run directory
        #----------------------------------------------------------------------
        outputDir = fullTestDir + test.name + '/'
        os.mkdir(outputDir)
    

        #----------------------------------------------------------------------
        # compile the code
        #----------------------------------------------------------------------
        if test.useExtraBuildDir > 0:
            bDir = suite.extraBuildDirs[test.useExtraBuildDir-1] + test.buildDir
        else:
            bDir = suite.sourceDir + test.buildDir

        os.chdir(bDir)

        if (test.reClean == 1):
            # for one reason or another, multiple tests use different
            # build options, make clean again to be safe
            print "  re-making clean..."

            systemCall("%s BOXLIB_HOME=%s %s %s realclean >& /dev/null" % 
                       (suite.MAKE, suite.boxLibDir, 
                        suite.extSrcCompString, suite.extraBuildDirCompString))

        
        print "  building..."

        if (suite.sourceTree == "C_Src" or test.testSrcTree == "C_Src"):

            buildOptions = ""
            exeSuffix = ""

            if (test.debug):
                buildOptions += "DEBUG=TRUE "
                exeSuffix += ".DEBUG"
            else:
                buildOptions += "DEBUG=FALSE "

            if (test.useMPI):
                buildOptions += "USE_MPI=TRUE "
                exeSuffix += ".MPI"
            else:
                buildOptions += "USE_MPI=FALSE "

            if (test.useOMP):
                buildOptions += "USE_OMP=TRUE "
                exeSuffix += ".OMP"
            else:
                buildOptions += "USE_OMP=FALSE "

            if test.useExtraBuildDir > 0:
                buildOptions += suite.extraBuildDirCompString + " "

            executable = "%s%dd" % (suite.suiteName, test.dim) + exeSuffix + ".ex"

            compString = "%s -j%s BOXLIB_HOME=%s %s %s DIM=%d %s COMP=%s FCOMP=%s executable=%s  >& %s/%s.make.out" % \
                (suite.MAKE, suite.numMakeJobs, suite.boxLibDir, 
                 suite.extSrcCompString, test.addToCompileString, 
                 test.dim, buildOptions, suite.COMP, suite.FCOMP, 
                 executable, outputDir, test.name)

            print "    " + compString
            systemCall(compString)            
	       
            
        elif (suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src"):

            buildOptions = ""

            if (test.debug):
                buildOptions += "NDEBUG= "
            else:
                buildOptions += "NDEBUG=t "

            if (test.useMPI):
                buildOptions += "MPI=t "
            else:
                buildOptions += "MPI= "

            if (test.useOMP):
                buildOptions += "OMP=t "
            else:
                buildOptions += "OMP= "
            
            if test.useExtraBuildDir > 0:
                buildOptions += suite.extraBuildDirCompString + " "

            compString = "%s -j%s BOXLIB_HOME=%s %s %s %s COMP=%s >& %s/%s.make.out" % \
                (suite.MAKE, suite.numMakeJobs, suite.boxLibDir, 
                 suite.extSrcCompString, test.addToCompileString, 
                 buildOptions, suite.FCOMP, outputDir, test.name)

            print "    " + compString
            systemCall(compString)


            # we need a better way to get the executable name here
            executable = getRecentFileName(bDir,"main",".exe")

        

        if (test.compileTest):
            
            # compilation tests are done now -- just make the report and ...
            shutil.copy("%s/%s.make.out"    % (outputDir, test.name), fullWebDir)

            print "  creating problem test report ..."
            reportSingleTest(suite, test, compString, "", testDir, fullWebDir)

            # ... skip to the next test in the loop
            continue
            
            
        #----------------------------------------------------------------------
        # copy the necessary files over to the run directory
        #----------------------------------------------------------------------
        print "  copying files to run directory..."

        try: shutil.copy(executable, outputDir)
        except IOError:

            # compilation failed.  First copy the make.out into the
            # web directory and then report
            shutil.copy("%s/%s.make.out" % (outputDir, test.name), fullWebDir)

            errorMsg = "    ERROR: compilation failed"
            reportTestFailure(errorMsg, test, testDir, fullWebDir, compString=compString)
            continue

        try: shutil.copy(test.inputFile, outputDir)
        except IOError:
            errorMsg = "    ERROR: unable to copy input file: %s" % test.inputFile
            reportTestFailure(errorMsg, test, testDir, fullWebDir)
            continue

	# sometimes the input file was in a subdirectory under the
	# build directory.  Keep only the input file for latter
	index = string.rfind(test.inputFile, "/")
	if (index > 0):
	   test.inputFile = test.inputFile[index+1:]


        # if we are a "C_Src" build, we need the probin file
        if (suite.sourceTree == "C_Src" or \
                (test.testSrcTree == "C_Src" and test.probinFile != "")):
            try: shutil.copy(test.probinFile, outputDir)
            except IOError:
                errorMsg = "    ERROR: unable to copy probin file: %s" % test.probinFile
                reportTestFailure(errorMsg, test, testDir, fullWebDir)
                continue

            # sometimes the probin file was in a subdirectory under the
            # build directory.  Keep only the probin file for latter
            index = string.rfind(test.probinFile, "/")
            if (index > 0):
               test.probinFile = test.probinFile[index+1:]	   


        # if we are using the Helmholtz EOS, we need the input table
        if (test.needsHelmEOS):
            helmeosFile = suite.helmeosDir + "helm_table.dat"
            try: os.symlink(helmeosFile, outputDir + "helm_table.dat")
            except IOError:
                errorMsg = "    ERROR: unable to links helmeos file: %s" % helmeosFile
                reportTestFailure(errorMsg, test, testDir, fullWebDir)
                continue


        # python doesn't allow labelled continue statements, so we
        # use skip_to_next_test to decide if we need to skip to 
        # the next test
        skip_to_next_test = 0
        for file in test.auxFiles:
            try: shutil.copy(file, outputDir)
            except IOError:
                errorMsg = "    ERROR: unable to copy aux file: %s" % file
                reportTestFailure(errorMsg, test, testDir, fullWebDir)
                skip_to_next_test = 1
                break
            
        if (skip_to_next_test):
            continue


        # python doesn't allow labelled continue statements, so we
        # use skip_to_next_test to decide if we need to skip to 
        # the next test
        skip_to_next_test = 0
        for file in test.linkFiles:
            if (not os.path.isfile(file)):
                errorMsg = "    ERROR: link file %s does not exist" % file
                reportTestFailure(errorMsg, test, testDir, fullWebDir)
                skip_to_next_test = 1
                break

            else:
                if os.path.isabs(file):
                    link_source = file
                    link_name = outputDir + os.path.basename(file)
                else:
                    link_source = os.path.abspath(file)
                    link_name = outputDir + file
                try: os.symlink(link_source, link_name)
                except IOError:
                    errorMsg = "    ERROR: unable to symlink link file: %s" % file
                    reportTestFailure(errorMsg, test, testDir, fullWebDir)
                    skip_to_next_test = 1
                    break
            
        if (skip_to_next_test):
            continue


        #----------------------------------------------------------------------
        # run the test
        #----------------------------------------------------------------------
        print "  running the test..."

        os.chdir(outputDir)

        test.wallTime = time.time()

        if (suite.sourceTree == "C_Src" or test.testSrcTree == "C_Src"):

	    if (test.useMPI and test.useOMP):

	       # create the MPI executable
	       testRunCommand = "OMP_NUM_THREADS=%s %s" % (test.numthreads, suite.MPIcommand)
	       testRunCommand = testRunCommand.replace("@host@", suite.MPIhost)
	       testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (test.numprocs))

               # keep around the checkpoint files only for the restart runs
               if (test.restartTest):
                   command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=1 amr.check_int=%d >&  %s.run.out < /dev/null" % \
                       (executable, test.inputFile, test.name, test.name, test.restartFileNum, test.name)
               else:
                   command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 >&  %s.run.out < /dev/null" % \
                       (executable, test.inputFile, test.name, test.name, test.name)
                   
               testRunCommand = testRunCommand.replace("@command@", command)
                   
               print "    " + testRunCommand
               systemCall(testRunCommand)

	    elif (test.useMPI):

	       # create the MPI executable
	       testRunCommand = suite.MPIcommand
	       testRunCommand = testRunCommand.replace("@host@", suite.MPIhost)
	       testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (test.numprocs))

                # keep around the checkpoint files only for the restart runs
	       if (test.restartTest):               
                   command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=1 amr.check_int=%d >&  %s.run.out < /dev/null" % \
                       (executable, test.inputFile, test.name, test.name, test.restartFileNum, test.name)
               else:
                   command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 >&  %s.run.out < /dev/null" % \
                       (executable, test.inputFile, test.name, test.name, test.name)
	       
	       testRunCommand = testRunCommand.replace("@command@", command)

	       print "    " + testRunCommand
               systemCall(testRunCommand)

	    elif (test.useOMP):

                # keep around the checkpoint files only for the restart runs
                if (test.restartTest):
                    testRunCommand = "OMP_NUM_THREADS=%s ./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=1 amr.check_int=%d >&  %s.run.out < /dev/null" % \
                        (test.numthreads, executable, test.inputFile, test.name, test.name, test.restartFileNum, test.name)
                else:
                    testRunCommand = "OMP_NUM_THREADS=%s ./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 >&  %s.run.out < /dev/null" % \
                        (test.numthreads, executable, test.inputFile, test.name, test.name, test.name)
	       
                print "    " + testRunCommand
                systemCall(testRunCommand)
	       
            else:

                # keep around the checkpoint files only for the restart runs
                if (test.restartTest):
                    testRunCommand = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=1 amr.check_int=%d >&  %s.run.out" % \
                        (executable, test.inputFile, test.name, test.name, test.restartFileNum, test.name)
                else:
                    testRunCommand = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 >&  %s.run.out" % \
                        (executable, test.inputFile, test.name, test.name, test.name)

                print "    " + testRunCommand
                systemCall(testRunCommand)

        elif (suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src"):

            if (test.useMPI and test.useOMP):

                # create the MPI executable
                testRunCommand = "OMP_NUM_THREADS=%s %s" % (test.numthreads, suite.MPIcommand)
                testRunCommand = testRunCommand.replace("@host@", suite.MPIhost)
                testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (test.numprocs))

                # keep around the checkpoint files only for the restart runs
                if (test.restartTest):
                    command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk %s >& %s.run.out" % \
                        (executable, test.inputFile, test.name, test.name, suite.globalAddToExecString, test.name)
                else:
                    command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 %s >& %s.run.out" % \
                        (executable, test.inputFile, test.name, test.name, suite.globalAddToExecString, test.name)

                testRunCommand = testRunCommand.replace("@command@", command)

                print "    " + testRunCommand
                systemCall(testRunCommand)

            elif (test.useMPI):

	       # create the MPI executable
	       testRunCommand = suite.MPIcommand
	       testRunCommand = testRunCommand.replace("@host@", suite.MPIhost)
	       testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (test.numprocs))

                # keep around the checkpoint files only for the restart runs
	       if (test.restartTest):
                   command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk %s >& %s.run.out" % \
                       (executable, test.inputFile, test.name, test.name, suite.globalAddToExecString, test.name)
               else:
                   command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 %s >& %s.run.out" % \
                       (executable, test.inputFile, test.name, test.name, suite.globalAddToExecString, test.name)

	       testRunCommand = testRunCommand.replace("@command@", command)

	       print "    " + testRunCommand
               systemCall(testRunCommand)

            elif (test.useOMP):

                # keep around the checkpoint files only for the restart runs
                if (test.restartTest):
                    testRunCommand = "OMP_NUM_THREADS=%s ./%s %s --plot_base_name %s_plt --check_base_name %s_chk %s >& %s.run.out" % \
                        (test.numthreads, executable, test.inputFile, test.name, test.name, suite.globalAddToExecString, test.name)
                else:
                    testRunCommand = "OMP_NUM_THREADS=%s ./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 %s >& %s.run.out" % \
                        (test.numthreads, executable, test.inputFile, test.name, test.name, suite.globalAddToExecString, test.name)

                print "    " + testRunCommand
                systemCall(testRunCommand)


            else:

                # keep around the checkpoint files only for the restart runs
                if (test.restartTest):
                    testRunCommand = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk %s >& %s.run.out" % \
                        (executable, test.inputFile, test.name, test.name, suite.globalAddToExecString, test.name)
                else:
                    testRunCommand = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 %s >& %s.run.out" % \
                        (executable, test.inputFile, test.name, test.name, suite.globalAddToExecString, test.name)

                print "    " + testRunCommand
                systemCall(testRunCommand)


        # if it is a restart test, then rename the final output file and
        # restart the test
        if (test.restartTest):
            lastFile = getLastPlotfile(outputDir, test)

            if (lastFile == ""):
                errorMsg = "ERROR: test did not produce output.  Restart test not possible"
                reportTestFailure(errorMsg, test, testDir, fullWebDir)
                continue

            origLastFile = "orig_%s" % (lastFile)
            shutil.move(lastFile, origLastFile)

            if test.diffDir:
                origDiffDir = "orig_%s" % (test.diffDir)
                shutil.move(test.diffDir, origDiffDir)

            # get the file number to restart from
            restartFile = "%s_chk%5.5d" % (test.name, test.restartFileNum)

            print "  restarting from %s ... " % (restartFile)
           
            if (suite.sourceTree == "C_Src" or test.testSrcTree == "C_Src"):
  
                if (test.useMPI and test.useOMP):

                    # create the MPI executable
                    testRunCommand = "OMP_NUM_THREADS=%s %s" % (test.numthreads, suite.MPIcommand)
                    testRunCommand = testRunCommand.replace("@host@", suite.MPIhost)
                    testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (test.numprocs))
                    
                    command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 amr.restart=%s >> %s.run.out 2>&1" % \
                        (executable, test.inputFile, test.name, test.name, restartFile, test.name)
                    
                    testRunCommand = testRunCommand.replace("@command@", command)
                    
                    print "    " + testRunCommand
                    systemCall(testRunCommand)

                elif (test.useMPI):

                    # create the MPI executable
                    testRunCommand = suite.MPIcommand
                    testRunCommand = testRunCommand.replace("@host@", suite.MPIhost)
                    testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (test.numprocs))
                    
                    command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 amr.restart=%s >> %s.run.out 2>&1" % \
                        (executable, test.inputFile, test.name, test.name, restartFile, test.name)
                    
                    testRunCommand = testRunCommand.replace("@command@", command)

                    print "    " + testRunCommand
                    systemCall(testRunCommand)

                elif (test.useOMP):

                    testRunCommand = "OMP_NUM_THREADS=%s ./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 amr.restart=%s >>  %s.run.out 2>&1" % \
                        (test.numthreads, executable, test.inputFile, test.name, test.name, restartFile, test.name)
	       
                    print "    " + testRunCommand
                    systemCall(testRunCommand)
                    
                else:
                    
                    testRunCommand = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 amr.restart=%s >> %s.run.out 2>&1" % \
                        (executable, test.inputFile, test.name, test.name, restartFile, test.name)
                    
                    print "    " + testRunCommand
                    systemCall(testRunCommand)

            elif (suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src"):

                if (test.useMPI and test.useOMP):

                    # create the MPI executable
                    testRunCommand = "OMP_NUM_THREADS=%s %s" % (test.numthreads, suite.MPIcommand)
                    testRunCommand = testRunCommand.replace("@host@", suite.MPIhost)
                    testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (test.numprocs))

                    command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 --restart %d %s >> %s.run.out 2>&1" % \
                        (executable, test.inputFile, test.name, test.name, test.restartFileNum, suite.globalAddToExecString, test.name)

                    testRunCommand = testRunCommand.replace("@command@", command)
                    
                    print "    " + testRunCommand
                    systemCall(testRunCommand)

                elif (test.useMPI):

                    # create the MPI executable
                    testRunCommand = suite.MPIcommand
                    testRunCommand = testRunCommand.replace("@host@", suite.MPIhost)
                    testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (test.numprocs))
                    
                    command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 --restart %d %s >> %s.run.out 2>&1" % \
                        (executable, test.inputFile, test.name, test.name, test.restartFileNum, suite.globalAddToExecString, test.name)

                    testRunCommand = testRunCommand.replace("@command@", command)
                    
                    print "    " + testRunCommand
                    systemCall(testRunCommand)

                elif (test.useOMP):

                    testRunCommand = "OMP_NUM_THREADS=%s ./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 --restart %d %s >> %s.run.out 2>&1" % \
                        (test.numthreads, executable, test.inputFile, test.name, test.name, test.restartFileNum, suite.globalAddToExecString, test.name)

                    print "    " + testRunCommand
                    systemCall(testRunCommand)


                else:

                    testRunCommand = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 --restart %d %s >> %s.run.out 2>&1" % \
                        (executable, test.inputFile, test.name, test.name, test.restartFileNum, suite.globalAddToExecString, test.name)

                    print "    " + testRunCommand
                    systemCall(testRunCommand)


        test.wallTime = time.time() - test.wallTime
           
            
        #----------------------------------------------------------------------
        # do the comparison
        #----------------------------------------------------------------------
        if (not test.selfTest):

            if (test.outputFile == ""):
                if (test.compareFile == ""):
                    compareFile = getLastPlotfile(outputDir, test)
                else:
                    compareFile = test.compareFile
                outputFile = compareFile
            else:
                outputFile = test.outputFile
                compareFile = test.name+'_'+outputFile


            if (not make_benchmarks):

                print "  doing the comparison..."
                print "    comparison file: ", outputFile

                if (not test.restartTest):
                    benchFile = benchDir + compareFile
                else:
                    benchFile = origLastFile

                # see if it exists
                # note, with BoxLib, the plotfiles are actually directories
            
                if (not os.path.isdir(benchFile)):
                    warning("    WARNING: no corresponding benchmark found")
                    benchFile = ""

                    cf = open("%s.compare.out" % (test.name), 'w')
                    cf.write("WARNING: no corresponding benchmark found\n")
                    cf.write("         unable to do a comparison\n")
                    cf.close()
                    
                else:
                    if (not compareFile == ""):

                        print "    benchmark file: ", benchFile
                   
                        command = "../fcompare.exe -n 0 --infile1 %s --infile2 %s >> %s.compare.out 2>&1" % (benchFile, outputFile, test.name)

                        cf = open("%s.compare.out" % (test.name), 'w')
                        cf.write(command)
                        cf.write("\n")
                        cf.close()
                    
                        systemCall(command)

                    else:
                        warning("WARNING: unable to do a comparison")

                        cf = open("%s.compare.out" % (test.name), 'w')
                        cf.write("WARNING: run did not produce any output\n")
                        cf.write("         unable to do a comparison\n")
                        cf.close()

                if (not test.diffDir == ""):
                    if not test.restartTest:
                        diffDirBench = benchDir + '/' + test.name + '_' + test.diffDir
                    else:
                        diffDirBench = origDiffDir
                    
                    print "  doing the diff..."
                    print "    diff dir: ", test.diffDir

                    command = "diff %s -r %s %s >> %s.compare.out 2>&1" \
                        % (test.diffOpts, diffDirBench, test.diffDir, test.name)

                    cf = open("%s.compare.out" % (test.name), 'a')
                    cf.write("\n\n")
                    cf.write(command)
                    cf.write("\n")
                    cf.close()

                    diffstatus = systemCall(command)

                    if (diffstatus == 0):
                        cf = open("%s.compare.out" % (test.name), 'a')
                        cf.write("diff was SUCCESSFUL\n")
                        cf.close()


            else:   # make_benchmarks

                print "  storing output of %s as the new benchmark..." % (test.name)
                print "     new benchmark file: ", compareFile

                if (not compareFile == ""):
                    if outputFile != compareFile:
                        systemCall("rm -rf %s/%s" % (benchDir, compareFile))
                        systemCall("cp -rf %s %s/%s" % (outputFile, benchDir, compareFile))
                    else:
                        systemCall("cp -rf %s %s" % (compareFile, benchDir))

                    cf = open("%s.status" % (test.name), 'w')
                    cf.write("benchmarks updated.  New file:  %s\n" % (compareFile) )
                    cf.close()
                else:
                    cf = open("%s.status" % (test.name), 'w')
                    cf.write("benchmarks failed")
                    cf.close()

                if (not test.diffDir == ""):
                    diffDirBench = benchDir + '/' + test.name + '_' + test.diffDir
                    systemCall("rm -rf %s" % (diffDirBench))
                    print "     new diffDir: ", test.name + '_' + test.diffDir
                    systemCall("cp -r %s %s" % (test.diffDir, diffDirBench))


        else:   # selfTest

            if (not make_benchmarks):

                print "  looking for selfTest success string: %s ..." % test.stSuccessString

                try:
                    of = open("%s.run.out" % (test.name), 'r')

                except IOError:
                    warning("WARNING: no output file found")
                    compareSuccessful = 0
                    outLines = ['']

                else:
                    outLines = of.readlines()

                    # successful comparison is indicated by PLOTFILES AGREE
                    compareSuccessful = 0
    
                    for line in outLines:
                        if (string.find(line, test.stSuccessString) >= 0):
                            compareSuccessful = 1
                            break
    
                    of.close()

                cf = open("%s.compare.out" % (test.name), 'w')

                if (compareSuccessful):
                    cf.write("SELF TEST SUCCESSFUL\n")
                else:
                    cf.write("SELF TEST FAILED\n")

                cf.close()


                
        #----------------------------------------------------------------------
        # do any requested visualization (2- and 3-d only)
        #---------------------------------------------------------------------- 
        if (test.doVis and not make_benchmarks):

            if (not outputFile == ""):

                print "  doing the visualization..."

                if (test.dim == 2):
                    systemCall('%s/%s --palette %s/Palette -cname "%s" -p "%s" >& /dev/null' %
                              (suite.compareToolDir, vis2dExecutable, suite.compareToolDir, 
                               test.visVar, outputFile) )
                elif (test.dim == 3):
                    systemCall('%s/%s --palette %s/Palette -n 1 -cname "%s" -p "%s" >& /dev/null' %
                              (suite.compareToolDir, vis3dExecutable, suite.compareToolDir, 
                               test.visVar, outputFile) )
                else:
                    print "    Visualization not supported for dim = %d" % (test.dim)
                    

                # convert the .ppm files into .png files
                ppmFile = getRecentFileName(outputDir, "", ".ppm")

                systemCall("convert %s `basename %s .ppm`.png" % 
                          (ppmFile, ppmFile) )
        
            else:
                warning("    WARNING: no output file.  Skipping visualization")
        


        #----------------------------------------------------------------------
        # do any analysis
        #---------------------------------------------------------------------- 
        if (not test.analysisRoutine == "" and not make_benchmarks):

            if (not outputFile == ""):

                print "  doing the analysis..."
                if test.useExtraBuildDir > 0:
                    shutil.copy("{}/{}".format(suite.extraBuildDirs[test.useExtraBuildDir-1], 
                                               test.analysisRoutine),
                                os.getcwd())
                else:
                    shutil.copy("{}/{}".format(suite.sourceDir, test.analysisRoutine),
                                os.getcwd())
                            

                option = eval("suite.{}".format(test.analysisMainArgs))
                
                systemCall("{} {} {}".format(os.path.basename(test.analysisRoutine), 
                                             option, outputFile))

            else:
                warning("    WARNING: no output file.  Skipping visualization")
        
                       
        #----------------------------------------------------------------------
        # move the output files into the web directory
        #----------------------------------------------------------------------
        if (not make_benchmarks):
            shutil.copy("%s.run.out"     % (test.name), fullWebDir)
            shutil.copy("%s.make.out"    % (test.name), fullWebDir)
            shutil.copy("%s.compare.out" % (test.name), fullWebDir)

            shutil.copy(test.inputFile, "%s/%s.%s" % (fullWebDir, test.name, test.inputFile) )

            if (suite.sourceTree == "C_Src"):
                shutil.copy(test.probinFile, "%s/%s.%s" % (fullWebDir, test.name, test.probinFile) )

            for file in test.auxFiles:

                # sometimes the auxFile was in a subdirectory under the
                # build directory.  
                index = string.rfind(file, "/")
                if (index > 0):
                    file = file[index+1:]	   

                shutil.copy(file, "%s/%s.%s" % (fullWebDir, test.name, file) )

            if (test.doVis):
               pngFile = getRecentFileName(outputDir, "", ".png")
               try: shutil.copy(pngFile, fullWebDir)
               except IOError:
                   # visualization was not successful.  Reset doVis
                   test.doVis = 0

            if (not test.analysisRoutine == ""):
                try: shutil.copy(test.analysisOutputImage, fullWebDir)
                except IOError:
                    # analysis was not successful.  Reset the output image
                    test.analysisOutputImage = ""

               
        else:
            shutil.copy("%s.status" % (test.name), fullWebDir)
            

        #----------------------------------------------------------------------
        # archive the output
        #----------------------------------------------------------------------
        print "  archiving the output..."
        for file in os.listdir(outputDir):
            if (os.path.isdir(file) and 
                (file.startswith("%s_plt" % (test.name)) or 
                 file.startswith("%s_chk" % (test.name)) ) ):

                try:
                    tar = tarfile.open("%s.tgz" % (file), "w:gz")
                    tar.add("%s" % (file))
                    tar.close()

                except:
                    warning("    WARNING: unable to tar output file %s" % (file))

                else:
                    shutil.rmtree(file)
                    

        #----------------------------------------------------------------------
        # write the report for this test
        #----------------------------------------------------------------------
        if (not make_benchmarks):
            print "  creating problem test report ..."
            reportSingleTest(suite, test, compString, testRunCommand, testDir, fullWebDir)

        
    #--------------------------------------------------------------------------
    # write the report for this instance of the test suite
    #--------------------------------------------------------------------------
    print "\n"
    bold("creating new test report...")
    numFailed = reportThisTestRun(suite, make_benchmarks, comment, note, 
                                  updateTime,  updateBoxLib,
                                  updateSource, updateExtSrc,
                                  testList, testDir, testFile, fullWebDir)


    # make sure that all of the files in the web directory are world readable
    for file in os.listdir(fullWebDir):    
       currentFile = fullWebDir + file

       if (os.path.isfile(currentFile)):
          os.chmod(currentFile, 0644)
          
    if updateBoxLib or boxLibGitHash:
        doGITback(suite.boxLibDir, "BoxLib", boxLibGitBranch_Orig)

    if updateSource or sourceGitHash:
        doGITback(suite.sourceDir, suite.srcName, sourceGitBranch_Orig)

    if updateExtSrc or extSrcGitHash:
        doGITback(suite.extSrcDir, suite.extSrcName, extSrcGitBranch_Orig)
            
    #--------------------------------------------------------------------------
    # For temporary run, return now without creating suote report.
    #--------------------------------------------------------------------------
    if do_temp_run:
        return numFailed


    #--------------------------------------------------------------------------
    # generate the master report for all test instances 
    #--------------------------------------------------------------------------
    print "\n"
    bold("creating suite report...")
    tableHeight = min(max(suite.lenTestName, 4), 16)
    reportAllRuns(suite, activeTestList, suite.webTopDir, tableHeight=tableHeight)

    def emailDevelopers():
        msg = email.message_from_string(suite.emailBody)
        msg['From'] = suite.emailFrom
        msg['To'] = ",".join(suite.emailTo)
        msg['Subject'] = suite.emailSubject

        server = smtplib.SMTP('localhost')
        server.sendmail(suite.emailFrom, suite.emailTo, msg.as_string())
        server.quit()

    if (numFailed > 0 and suite.sendEmailWhenFail):
        print "\n"
        bold("sending email...")
        emailDevelopers()


    return numFailed


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# R E P O R T   W R I T I N G   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

cssContents = \
r"""
body {font-family: "Arial", san-serif;}

h1 {font-family: "Tahoma","Arial", sans-serif;
    color: #333333;}

h3 {display: inline;}

h3.passed {text-decoration: none; display: inline;
           color: black; background-color: lime; padding: 2px}

a.passed:link {color: black; text-decoration: none;}
a.passed:visited {color: black; text-decoration: none;}
a.passed:hover {color: #ee00ee; text-decoration: underline;}

h3.failed {text-decoration: none; display: inline; 
           color: black; background-color: red; padding: 2px}

a.failed:link {color: yellow; text-decoration: none;}
a.failed:visited {color: yellow; text-decoration: none;}
a.failed:hover {color: #00ffff; text-decoration: underline;}


h3.benchmade {text-decoration: none; display: inline; 
              color: black; background-color: orange; padding: 2px}

a.benchmade:link {color: black; text-decoration: none;}
a.benchmade:visited {color: black; text-decoration: none;}
a.benchmade:hover {color: #00ffff; text-decoration: underline;}


span.nobreak {white-space: nowrap;}

a.main:link {color: yellow; text-decoration: none;}
a.main:visited {color: yellow; text-decoration: none;}
a.main:hover {color: #00ffff; text-decoration: underline;}

th {background-color: black;
    padding: 4px;
    color: yellow;
    border-width: 0px;}

td {border-width: 0px;
    padding: 5px;
    background-color: white;
    vertical-align: middle;}

td.passed {background-color: lime; opacity: 0.8;}
td.failed {background-color: red; color: yellow; opacity: 0.8;}
td.benchmade {background-color: orange; opacity: 0.8;}
td.date {background-color: #666666; color: white; opacity: 0.8; font-weight: bold;}

.maintable tr:hover {background-color: blue;}


table {border-collapse: separate;
       border-spacing: 2px;
       margin-left: auto;
       margin-right: auto;
       border-width: 1px;
       border-color: gray;
       border-style: solid;
       box-shadow: 10px 10px 5px #888888;}

/* http://blog.petermares.com/2010/10/27/vertical-text-in-html-table-headers-for-webkitmozilla-browsers-without-using-images/ */

div.verticaltext {text-align: center;
                  vertical-align: middle;
                  width: 20px;
                  margin: 0px;
                  padding: 0px;
                  padding-left: 3px;
                  padding-right: 3px;
                  padding-top: 10px;
                  white-space: nowrap;
                  -webkit-transform: rotate(-90deg); 
                  -moz-transform: rotate(-90deg);}

th {background-color: grey;
    color: yellow;
    text-align: center;
    vertical-align: bottom;
    height: @TABLEHEIGHT@;
    padding-bottom: 3px;
    padding-left: 5px;
    padding-right: 5px;}

"""

HTMLHeader = \
r"""
<HTML>
<HEAD>
<TITLE>@TESTDIR@ / @TESTNAME@</TITLE>
<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=ISO-8859-1">
<LINK REL="stylesheet" TYPE="text/css" HREF="tests.css">
</HEAD>
<BODY>
"""

MainHeader = \
r"""
<HTML>
<HEAD>
<TITLE>@TITLE@</TITLE>
<META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=ISO-8859-1">
<LINK REL="stylesheet" TYPE="text/css" HREF="tests.css">
</HEAD>
<BODY>
<!--GOUPLINK-->
<CENTER><H1>@TITLE@</H1></CENTER>
"""


#==============================================================================
# create_css
#==============================================================================
def create_css(tableHeight=16):
    """ write the css file for the webpages """
    
    cssC = cssContents.replace("@TABLEHEIGHT@", "%sem" % (tableHeight))

    cssFile = "tests.css"
    cf = open(cssFile, 'w')
    cf.write(cssC)
    cf.close()


#==============================================================================
# reportSingleTest
#==============================================================================
def reportSingleTest(suite, test, compileCommand, runCommand, testDir, fullWebDir):
    """ generate a single problem's test result page """
    
    # get the current directory
    currentDir = os.getcwd()
    
    # switch to the web directory and open the report file
    os.chdir(fullWebDir)

    #--------------------------------------------------------------------------
    # parse the compilation report and determine if we compiled
    #--------------------------------------------------------------------------
    compileFile = "%s.make.out" % (test.name)

    try:
        cf = open(compileFile, 'r')

    except IOError:
        warning("WARNING: no compilation file found")
        compileSuccessful = 0

    else:
        makeLines = cf.readlines()

        # successful compilation be indicated by SUCCESS or
        # Nothing to be done for `all'.  Look for both
        compileSuccessful = 0
    
        for line in makeLines:
            if (string.find(line, "SUCCESS") >= 0 or
                string.find(line, "is up to date.") or
                string.find(line, "Nothing to be done") >= 0):
                compileSuccessful = 1
                break
        
        cf.close()
            
    
    #--------------------------------------------------------------------------
    # parse the compare report and determine if we passed
    #--------------------------------------------------------------------------
    if (not test.compileTest):
        compareFile = "%s.compare.out" % (test.name)
    
        try:
            cf = open(compareFile, 'r')

        except IOError:
            warning("WARNING: no comparison file found")
            compareSuccessful = 0
            diffLines = ['']

        else:
            diffLines = cf.readlines()

            # successful comparison is indicated by PLOTFILES AGREE
            compareSuccessful = 0
    
            for line in diffLines:
                if (string.find(line, "PLOTFILES AGREE") >= 0 or
                    string.find(line, "SELF TEST SUCCESSFUL") >= 0):
                    compareSuccessful = 1
                    break

            if (compareSuccessful):
                if (not test.diffDir == ""):
                    compareSuccessful = 0
                    for line in diffLines:
                        if (string.find(line, "diff was SUCCESSFUL") >= 0):
                            compareSuccessful = 1
                            break
    
            cf.close()


    #--------------------------------------------------------------------------
    # write out the status file for this problem, with either
    # PASSED or FAILED
    #--------------------------------------------------------------------------
    statusFile = "%s.status" % (test.name)
    sf = open(statusFile, 'w')

    if (compileSuccessful and 
        (test.compileTest or
         (not test.compileTest and compareSuccessful))):
        sf.write("PASSED\n")
        success("    %s PASSED" % (test.name))
    else:
        sf.write("FAILED\n")
        testfail("    %s FAILED" % (test.name))

    sf.close()


    #--------------------------------------------------------------------------
    # generate the HTML page for this test
    #--------------------------------------------------------------------------

    # write the css file
    create_css()

    
    htmlFile = "%s.html" % (test.name)
    hf = open(htmlFile, 'w')

    newHead = HTMLHeader + r"""<CENTER><H1><A HREF="index.html">@TESTDIR@</A> / @TESTNAME@</H1></CENTER>"""

    newHead = newHead.replace("@TESTDIR@", os.path.normpath(testDir))
    newHead = newHead.replace("@TESTNAME@", test.name)

    hf.write(newHead)

    hf.write("<P><b>build directory:</b> %s" % (test.buildDir) )
    if test.useExtraBuildDir > 0:
        hf.write(" in %s\n" % (suite.extraBuildDirs[test.useExtraBuildDir-1]))
    else:
        hf.write("\n")


    hf.write("<P>&nbsp;\n")

    if (not test.compileTest):

        if (test.debug):
            hf.write("<p><b>Debug test</b>\n")
            hf.write("<p>&nbsp;\n")

        if (test.useMPI):

            hf.write("<P><b>Parallel (MPI) Run</b><br>numprocs = %d\n" % (test.numprocs) )
            hf.write("<P>&nbsp;\n")


        if (test.useOMP):

            hf.write("<P><b>OpenMP Run</b><br>numthreads = %d\n" % (test.numthreads) )
            hf.write("<P>&nbsp;\n")


        hf.write("<p><b>Execution Time</b> (seconds) = %f\n" % (test.wallTime))


        # is this a restart test?
        if (test.restartTest):

            hf.write("<P><b>Restart Test</b><br>Job was run as normal and then restarted from checkpoint # %d, and the two final outputs were compared\n" % (test.restartFileNum) )

        hf.write("<P>&nbsp;\n")       

        # write out the information about the test
        hf.write("<P><b>input file:</b> <A HREF=\"%s.%s\">%s</A>\n" %
                 (test.name, test.inputFile, test.inputFile) )


        if (suite.sourceTree == "C_Src"):
            hf.write("<P><b>probin file:</b> <A HREF=\"%s.%s\">%s</A>\n" %
                     (test.name, test.probinFile, test.probinFile) )    


        i = 1
        for file in test.auxFiles:
            hf.write("<P><b>aux%dFile:</b> <A HREF=\"%s.%s\">%s</A>\n" %
                     (i, test.name, file, file) )
            i = i + 1
        

        hf.write("<P><b>dimensionality:</b> %s\n" % (test.dim) )

        hf.write("<P>&nbsp;\n")

    
    # write out the compilation report
    if (compileSuccessful):
        hf.write("<P><H3 CLASS=\"passed\">Compilation Successful</H3></P>\n")
    else:
        hf.write("<P><H3 CLASS=\"failed\">Compilation Failed</H3></P>\n")

    hf.write("<P>Compliation command:<br>\n&nbsp; %s\n" % (compileCommand) )
    hf.write("<P><A HREF=\"%s.make.out\">make output</A>\n" % (test.name) )

    hf.write("<P>&nbsp;\n")

    
    if (not test.compileTest):

        # write out the comparison report
        if (compareSuccessful):
            hf.write("<P><H3 CLASS=\"passed\">Comparison Successful</H3></P>\n")
        else:
            hf.write("<P><H3 CLASS=\"failed\">Comparison Failed</H3></P>\n")

        hf.write("<P>Execution command:<br>\n&nbsp; %s\n" % (runCommand) )
        hf.write("<P><A HREF=\"%s.run.out\">execution output</A>\n" % (test.name) )


        hf.write("<P>&nbsp;\n")
        hf.write("<PRE>\n")
    
        for line in diffLines:
            hf.write(line.replace('<','&lt;').replace('>','&gt;'))

        hf.write("</PRE>\n")

        # show any visualizations
        if (test.doVis):
            pngFile = getRecentFileName(fullWebDir, test.name, ".png")
            hf.write("<P>&nbsp;\n")
            hf.write("<P><IMG SRC='%s' BORDER=0>" % (pngFile) )

        # show any analysis
        if (not test.analysisOutputImage == ""):
            hf.write("<P>&nbsp;\n")
            hf.write("<P><IMG SRC='%s' BORDER=0>" % (test.analysisOutputImage) )
    

    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")    

    hf.close()
    
    
    # switch back to the original directory
    os.chdir(currentDir)
	

#==============================================================================
# reportTestAbort
#==============================================================================
def reportTestFailure(message, test, testDir, fullWebDir, compString=None):
    """ generate a simple report for an error encountered while performing 
        the test """
    
    testfail("    aborting test")
    testfail(message)

    # get the current directory
    currentDir = os.getcwd()
    
    # switch to the web directory and open the report file
    os.chdir(fullWebDir)

    #--------------------------------------------------------------------------
    # write out the status file for this problem -- FAILED
    #--------------------------------------------------------------------------
    statusFile = "%s.status" % (test.name)
    sf = open(statusFile, 'w')
    sf.write("FAILED\n")

    sf.close()

    testfail("    %s FAILED" % (test.name))


    #--------------------------------------------------------------------------
    # generate the HTML page for this test
    #--------------------------------------------------------------------------

    # check to see if the CSS file is present, if not, write it
    if (not os.path.isfile("tests.css")):
        create_css()

    
    htmlFile = "%s.html" % (test.name)
    hf = open(htmlFile, 'w')

    newHead = HTMLHeader + r"""<CENTER><H1><A HREF="index.html">@TESTDIR@</A> / @TESTNAME@</H1></CENTER>"""

    newHead = newHead.replace("@TESTDIR@", testDir)
    newHead = newHead.replace("@TESTNAME@", test.name)

    hf.write(newHead)

    # write out the information about the test
    hf.write("<P><b>build directory:</b> %s\n" % (test.buildDir) )

    hf.write("<P><H3 CLASS=\"failed\">Test Failed</H3></P>\n")
    hf.write("<P>%s</P>\n" % (message) )

    if (not compString == None):
        hf.write("<P>compliation command:\n %s\n" % (compString) )
        hf.write("<P><A HREF=\"%s.make.out\">make output</A>\n" % (test.name) )


    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")    

    hf.close()
    
    
    # switch back to the original directory
    os.chdir(currentDir)
	

#==============================================================================
# reportThisTestRun
#==============================================================================
def reportThisTestRun(suite, make_benchmarks, comment, note, updateTime,
                      updateBoxLib, updateSource, updateExtSrc,
                      testList, testDir, testFile, fullWebDir):
    """ generate the master page for a single run of the test suite """
    
    # get the current directory
    currentDir = os.getcwd()
    
    # switch to the web directory and open the report file
    os.chdir(fullWebDir)


    # keep track of the number of tests that passed and the number that failed
    numFailed = 0
    numPassed = 0

    
    #--------------------------------------------------------------------------
    # generate the HTML page for this run of the test suite
    #--------------------------------------------------------------------------

    # always create the css (in case it changes)
    create_css()


    # create the master filename
    htmlFile = "index.html"
    
    hf = open(htmlFile, 'w')

    newHead = HTMLHeader + r"""<CENTER><H1><A HREF="../">@TESTDIR@</A> / @TESTNAME@</H1></CENTER>"""

    newHead = newHead.replace("@TESTDIR@", suite.suiteName)
    newHead = newHead.replace("@TESTNAME@", testDir)

    hf.write(newHead)

    if (not note == ""):
       hf.write("<p><b>Test run note:</b><br><font color=\"gray\">%s</font>\n" % (note) )
       
    if (make_benchmarks):
       hf.write("<p><b>Benchmarks updated</b><br>comment: <font color=\"gray\">%s</font>\n" % (comment) )
       hf.write("<p>&nbsp;\n")

    
    hf.write("<p><b>test input parameter file:</b> <A HREF=\"%s\">%s</A>\n" %
             (testFile, testFile) )

    if updateBoxLib or updateSource or updateExtSrc:
        hf.write("<p>&nbsp;\n")
        hf.write("<p><b>Git update was done at: </b>%s\n" % (updateTime) )

        if updateSource:
            hf.write("<p>&nbsp;&nbsp;<b>%s ChangeLog:</b> <A HREF=\"%s\">%s</A>\n" %
                     (suite.srcName, "ChangeLog."+suite.srcName, "ChangeLog."+suite.srcName) )        

        if updateBoxLib:
            hf.write("<p>&nbsp;&nbsp;<b>BoxLib ChangeLog:</b> <A HREF=\"%s\">%s</A>\n" %
                     ("ChangeLog.BoxLib", "ChangeLog.BoxLib") )        

        if updateExtSrc:
            hf.write("<p>&nbsp;&nbsp;<b>%s ChangeLog:</b> <A HREF=\"%s\">%s</A>\n" %
                     (suite.extSrcName, "ChangeLog."+suite.extSrcName, "ChangeLog."+suite.extSrcName) )        
    else:
        hf.write("<p>No git update done\n")

    hf.write("<p>&nbsp;\n")    

    hf.write("<P><TABLE>\n")
    
    # loop over the tests and add a line for each
    for test in testList:

        if (not make_benchmarks):
            
            # check if it passed or failed
            statusFile = "%s.status" % (test.name)

            sf = open(statusFile, 'r')
            lines = sf.readlines()

            testPassed = 0
        
            for line in lines:
                if (string.find(line, "PASSED") >= 0):
                    testPassed = 1
                    numPassed += 1
                    break

            if (not testPassed):
                numFailed += 1

        
            sf.close()

            # write out this test's status
            hf.write("<TR><TD><A HREF=\"%s.html\">%s</A></TD><TD>&nbsp;</TD>" %
                     (test.name, test.name) )
        
            if (testPassed):
                hf.write("<TD><H3 class=\"passed\">PASSED</H3></TD></TR>\n")
            else:
                hf.write("<TD><H3 class=\"failed\">FAILED</H3></TD></TR>\n")

        
            hf.write("<TR><TD>&nbsp;</TD></TR>\n")


        else:

            if (test.restartTest):
                continue

            if (test.compileTest):
                continue

            if (test.selfTest):
                continue

                
            # the benchmark was updated -- find the name of the new benchmark file
            benchStatusFile = "%s.status" % (test.name)

            bf = open(benchStatusFile, 'r')
            lines = bf.readlines()

            benchFile = "none"

            for line in lines:
                index = string.find(line, "file:")
                if (index >= 0):
                    benchFile = line[index+5:]
                    break
                

            if (not benchFile == "none"):
                    
                 hf.write("<TR><TD>%s</TD><TD>&nbsp;</TD><TD><H3 class=\"benchmade\">BENCHMARK UPDATED</H3></TD><TD>&nbsp;</TD><TD>(new benchmark file is %s)</TD></TR>\n" %
                          (test.name, benchFile) )
            else:
                 hf.write("<TR><TD>%s</TD><TD>&nbsp;</TD><TD><H3 class=\"failed\">BENCHMARK NOT UPDATED</H3></TD><TD>&nbsp;</TD><TD>(compilation or execution failed)</TD></TR>\n" %
                          (test.name) )
                 
            hf.write("<TR><TD>&nbsp;</TD></TR>\n")


    hf.write("</TABLE>\n")    

    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")    

    hf.close()


    #--------------------------------------------------------------------------
    # write out a status file for all the tests
    #--------------------------------------------------------------------------
    
    index = string.find(testDir, "/")
    if index > 0:
        statusFile = testDir[0:index] + ".status"
    else:
        statusFile = testDir + ".status"

    sf = open(statusFile, 'w')

    if (not make_benchmarks):
        if (numFailed == 0):
            sf.write("ALL PASSED\n")
        elif (numFailed > 0 and numPassed > 0):
            sf.write("SOME FAILED\n")
        else:
            sf.write("ALL FAILED\n")

    else:
        sf.write("BENCHMARKS UPDATED\n")
        
    sf.close()
    
    
    # switch back to the original directory
    os.chdir(currentDir)

    return numFailed

	
#==============================================================================
# reportAllRuns
#==============================================================================
def reportAllRuns(suite, activeTestList, webTopDir, tableHeight=16):

    os.chdir(webTopDir)

    create_css(tableHeight=tableHeight)

    validDirs = []
    allTests = []
    

    #--------------------------------------------------------------------------
    # start by finding the list of valid test directories
    #--------------------------------------------------------------------------
    for file in os.listdir(webTopDir):

        # look for a directory of the form 20* (this will work up until 2099
        if (file.startswith("20") and os.path.isdir(file)):

            # look for the status file
            statusFile = file + '/' + file + '.status'

            if (os.path.isfile(statusFile)):
                validDirs.append(file)


    validDirs.sort()
    validDirs.reverse()
    

    #--------------------------------------------------------------------------
    # now find all of the unique problems in the test directories
    #--------------------------------------------------------------------------
    for dir in validDirs:

        for file in os.listdir(webTopDir + dir):

            if (file.endswith(".status") and not file.startswith("20")):

                index = string.rfind(file, ".status")
                testName = file[0:index]

                if (allTests.count(testName) == 0):
                    if (not suite.reportActiveTestsOnly) or (testName in activeTestList):
                        allTests.append(testName)


    allTests.sort()
    

    #--------------------------------------------------------------------------
    # generate the HTML
    #--------------------------------------------------------------------------
    numTests = len(allTests)

    htmlFile = "index.html"

    title = "%s regression tests" % (suite.suiteName)
    
    hf = open(htmlFile, 'w')

    header = MainHeader.replace("@TITLE@", title)
    if (suite.goUpLink):
        header2 = header.replace("<!--GOUPLINK-->", '<a href="../">GO UP</a>')
        hf.write(header2)
    else:
        hf.write(header)

    hf.write("<P><TABLE class='maintable'>\n")

    # write out the header
    hf.write("<TR><TH ALIGN=CENTER>date</TH>\n")
    for test in allTests:
        hf.write("<TH><div class='verticaltext'>%s</div></TH>\n" % (test))
    
    hf.write("</TR>\n")


    # loop over all the test runs 
    for dir in validDirs:

        # first look to see if there are any valid tests at all --
        # otherwise we don't do anything for this date
        valid = 0
        for test in allTests:
            statusFile = "%s/%s/%s.status" % (webTopDir, dir, test)
            if (os.path.isfile(statusFile)):
                valid = 1
                break

        if not valid:
            continue

        # write out the directory (date)
        hf.write("<TR><TD class='date'><SPAN CLASS='nobreak'><A class='main' HREF=\"%s/index.html\">%s&nbsp;</A></SPAN></TD>\n" %
                 (dir, dir) )
        
        for test in allTests:
            
            # look to see if the current test was part of this suite run
            statusFile = "%s/%s/%s.status" % (webTopDir, dir, test)

            status = 0

            if (os.path.isfile(statusFile)):            

                sf = open(statusFile, 'r')
                lines = sf.readlines()

                # status = -1 (failed); 1 (passed); 10 (benchmark update)
                status = -1  
                                    
                for line in lines:
                    if (string.find(line, "PASSED") >= 0):
                        status = 1
                        break
                    elif (string.find(line, "FAILED") >= 0):
                        status = -1
                        break
                    elif (string.find(line, "benchmarks updated") >= 0):
                        status = 10
                        break
                    
                sf.close()

                
            # write out this test's status
            if (status == 1):
                hf.write("<TD ALIGN=CENTER class=\"passed\"><H3><a href=\"%s/%s.html\" class=\"passed\">:)</a></H3></TD>\n" % (dir, test))
            elif (status == -1):
                hf.write("<TD ALIGN=CENTER class=\"failed\"><H3><a href=\"%s/%s.html\" class=\"failed\">&nbsp;!&nbsp;</a></H3></TD>\n" % (dir, test))
            elif (status == 10):
                hf.write("<TD ALIGN=CENTER class=\"benchmade\"><H3>U</H3></TD>\n")
            else:
                hf.write("<TD>&nbsp;</TD>\n")


        hf.write("</TR>\n\n")
        
    hf.write("</TABLE>\n")    

    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")    

    hf.close()



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# m a i n
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
if __name__== "__main__":

    testSuite(sys.argv)
