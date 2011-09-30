#!/usr/bin/env python

"""
A simple regression test framework for a BoxLib-based code

There are several major sections to this source: the runtime parameter
routines, the test suite routines, and the report generation routines.  They
are separated as such in this file.

This test framework understands source based out of the Parallel/ and
fParallel/ frameworks.

2011-09-23
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


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# T E S T   C L A S S E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class testObj:

    def __init__ (self, name):

        self.name = name
        self.buildDir = ""

        self.inputFile = ""
        self.probinFile = ""
        self.auxFiles = []
        self.linkFiles = []

        self.dim = -1

        self.needsHelmEOS = -1

        self.restartTest = 0
        self.restartFileNum = -1

        self.compileTest = 0

        self.selfTest = 0
        self.stSuccessString = ""

        self.useMPI = 0
        self.numprocs = -1

        self.doVis = 0
        self.visVar = ""



class suiteObj:

    def __init__ (self):

        self.suiteName = "testDefault"

        self.sourceTree = ""
        self.boxLibDir = ""
        self.sourceDir = ""
        self.testTopDir = ""
        self.compareToolDir = ""
        self.helmeosDir = ""

        self.MPIcommand = ""
        self.MPIhost = ""

        self.FCOMP = "gfortran"
        self.COMP = "g++"

        self.MAKE = "gmake"
        self.numMakeJobs = 1

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
            if (not (value == "Parallel" or value == "fParallel")):
                fail("ERROR: invalid sourceTree")
            else:
                mysuite.sourceTree = value

        elif (opt == "sourceDir"):
            mysuite.sourceDir = checkTestDir(value)

        elif (opt == "boxLibDir"):
            mysuite.boxLibDir = checkTestDir(value)

        elif (opt == "testTopDir"):
            mysuite.testTopDir = checkTestDir(value)

        elif (opt == "compareToolDir"):
            mysuite.compareToolDir = checkTestDir(value)

        elif (opt == "helmeosDir"):
            mysuite.helmeosDir = checkTestDir(value)

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

        else:
            warning("WARNING: suite parameter %s not valid" % (opt))
            

    # checks
    if (mysuite.sourceTree == "" or mysuite.boxLibDir == "" or
        mysuite.sourceDir == "" or mysuite.testTopDir == "" or
        mysuite.compareToolDir == "" or mysuite.helmeosDir == ""):
        fail("ERROR: required suite-wide directory not specified\n" + \
                 "(sourceTree, boxLibDir, sourceDir, testTopDir, compareToolDir, helmeosDir)")

    # all other sections are tests
    print "\n"
    bold("finding tests and checking parameters...")

    for sec in cp.sections():

        if (sec == "main"): continue

        print "  %s" % (sec)

        # create the test object for this test
        mytest = testObj(sec)

        # set the test object data by looking at all the options in
        # the current section of the parameter file
        for opt in cp.options(sec):

            # get the value of the current option
            value = convertType(cp.get(sec, opt))

            if (opt == "buildDir"):
                mytest.buildDir = value

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
                
            elif (opt == "useMPI"):
                mytest.useMPI = value
                
            elif (opt == "numprocs"):
                mytest.numprocs = value

            elif (opt == "doVis"):
                mytest.doVis = value

            elif (opt == "visVar"):
                mytest.visVar = value

            else:
                warning("   WARNING: unrecognized parameter %s for test %s" % (opt, sec))



        invalid = 0

        # make sure all the require parameters are present
        if (mytest.buildDir == "" or mytest.inputFile == "" or
            (mysuite.sourceTree == "Parallel" and mytest.probinFile == "") or 
            mytest.dim == -1 or mytest.needsHelmEOS == -1):
            warning("   WARNING: manditory runtime parameters for test %s not set" % (sec))
            invalid = 1


        # check the optional parameters
        if (mytest.restartTest and mytest.restartFileNum == -1):
            warning("   WARNING: test %s is a restart test but restartFileNum not set" % (sec))
            invalid = 1

        if (mytest.selfTest and mytest.stSuccessString == ""):
            warning("   WARNING: test %s is a self-test but stSuccessString not set" % (sec))
            invalid = 1

        if (mytest.useMPI and mytest.numprocs == -1):
            warning("   WARNING: test %s is a parallel test but numprocs not set" % (sec))
            invalid = 1

        if (mytest.doVis and mytest.visVar == ""):
            warning("   WARNING: test %s requested visualization visVar not set" % (sec))
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
        directories"""
    
    buildDirs = []
    for obj in testList:

        currentBuildDir = obj.buildDir
        if (buildDirs.count(currentBuildDir) == 0):
            buildDirs.append(currentBuildDir)

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
    """ for fParallel builds, given the base and extension, find the
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
# doCVSUpdate
#==============================================================================
def doCVSUpdate(topDir, root, outDir):
   """ do a CVS update of the repository named root.  topDir is the full path
       to the directory containing root.  outDir is the full path to the 
       directory where we will store the CVS output """

   os.chdir(topDir)
 
   print "\n"
   bold("cvs update %s" % (root))

   # we need to be tricky here to make sure that the stdin is presented to      
   # the user to get the password.  Therefore, we use the subprocess            
   # class instead of os.system                                                 
   prog = ["cvs", "update", "%s" % (root)]
   p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT)
   stdout, stderr = p.communicate()


   # check if CVS was successful and if so, write stdout (and error
   # which was combined with stdout) into a log file
   cvsFailed = 0
   for line in stdout:
       if (string.find(line, "update aborted") >= 0):
           cvsFailed = 1
           break

   try:
       cf = open("cvs.%s.out" % (root), 'w')
   
   except IOError:
       fail("  ERROR: unable to open file for writing")

   else:
       for line in stdout:
           cf.write(line)

       cf.close()


   if (cvsFailed or stdout == ""):
       fail("  ERROR: cvs update was aborted. See cvs.%s.out for details" % (root))
       
        
   shutil.copy("cvs.%s.out" % (root),  outDir)



#==============================================================================
# doGITUpdate
#==============================================================================
def doGITUpdate(topDir, root, outDir):
   """ do a git update of the repository in topDir.  root is the name
       of the directory (used for labeling).  outDir is the full path
       to the directory where we will store the git output"""

   os.chdir(topDir)
 
   print "\n"
   bold("'git pull' in %s" % (topDir))

   # we need to be tricky here to make sure that the stdin is presented to      
   # the user to get the password.  Therefore, we use the subprocess            
   # class instead of os.system                                                 
   prog = ["git", "pull"]
   p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT)
   stdout, stderr = p.communicate()


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


#==============================================================================
# makeCVSChangeLog
#==============================================================================
def makeCVSChangeLog(suite, root, outDir):
   """ generate a ChangeLog for the CVS repository named root.  outDir
       is the full path to the directory where we will store the CVS
       output"""

   os.chdir(suite.sourceDir)

   have_cvs2cl = 0 

   print "\n"
   bold("generating ChangeLog for %s/" % (root))
    
   if (not os.path.isfile("%s/cvs2cl.pl" % (root) )):
       if (not os.path.isfile("%s/Tools/F_scripts/cvs2cl.pl" % 
                              (suite.boxLibDir) )):
           warning("  WARNING: unable to locate cvs2cl.pl script.")
           warning("           no ChangeLog will be generated")
       else:
           shutil.copy("%s/Tools/F_scripts/cvs2cl.pl" % (suite.boxLibDir), "%s/" % (root) )
           have_cvs2cl = 1
   else:
       have_cvs2cl = 1

   os.chdir("%s/" % (root) )

   if (have_cvs2cl):
       systemCall("./cvs2cl.pl -f ChangeLog.%s >& /dev/null" % (root) )
   else:
       cf = open("ChangeLog.%s" % (root), 'w')
       cf.write("unable to generate ChangeLog")
       cf.close()

   os.chdir(suite.sourceDir)

   shutil.copy("%s/ChangeLog.%s" % (root, root), outDir)


#==============================================================================
# makeGITChangeLog
#==============================================================================
def makeGITChangeLog(suite, root, outDir):
   """ generate a ChangeLog git repository named root.  outDir is the
       full path to the directory where we will store the git output"""

   os.chdir(suite.boxLibDir)


   print "\n"
   bold("generating ChangeLog for %s/" % (root))
    
   systemCall("git log --name-only >& ChangeLog.%s" % (root) )

   shutil.copy("ChangeLog.%s" % (root), outDir)



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# test
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def testSuite(argv):

    usage = """
    ./test.py [--make_benchmarks comment,
               --no_update,
               --single_test test,
               --note note]
        testfile.ini


    arguments:

      testfile.ini 
          This is the input file that defines the tests that the
          suite knows about.  It has the format

            [main]
            boxLibDir      = < directory to the main BoxLib/ directory >
            sourceDir      = < directory above Parallel/ and fParallel/ >
            testTopDir     = < full path to test output directory >
            compareToolDir = < full path to the fParallel/data_processing/ directory >
            helmeosDir     = < full path to helm_table.dat >

            sourceTree = < Parallel or fParallel -- what type is it? >

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

            useMPI = <is this a parallel job? 0 for no, 1 for yes) >
            numprocs = < # of processors to run on (if parallel job) >

            compileTest = < 0 for normal run, 1 if we just test compilation >

            selfTest = < 0 for normal run, 1 if test self-diagnoses if it 
                         succeeded >
            stSuccessString = < string to search for in self-test output to 
                         determine success >

            doVis = < 0 for no visualization, 1 if we do visualization >
            visVar = < string of the variable to visualize >


          Here, [main] lists the parameters for the test suite as a
          whole and [Sod-x] is a single test.  There can be many more
          tests, each with their own unique name, following the format
          of [Sod-x] above.

          In [main],

            boxLibDir is the full path to the BoxLib/ directory.
            
            sourceDir should be the full path to the directory where
            the Parallel/ and fParallel/ source directories live.

            testTopDir is the full path to the directory that will
            contain the test run directories, the test web
            directories, and the test benchmark directories.  It can
            be the same as sourceDir.  The necessary sub- directories
            will be created at runtime if they don't already exist.

            compareToolDir is the full path to the directory
            containing the comparison tool.  This resides in
            Parallel/util/Convergence/

            helmeosDir is the full path to the directory containing
            the helm_table.dat file needed by the general EOS.

            sourceTree is either Parallel (for applications that live
            in the Parallel/ tree) or fParallel (for applications that
            live in the fParallel/ tree).

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
            probinFile is only required for a Parallel (not fParallel)
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

            useMPI is set to 1 if the job is parallel.  In this case, you
            also must specify the number of processors, via numprocs.

            compileTest is set to 0 for a normal test.  Setting to 1 will
            just test compilation, and not any subsequent output.

            selfTest is set to 0 for a normal test, in which case the
            test output will be compared to the stored benchmark.  Setting 
            to 1 will determine success by searching for the string
            stSuccessString in the execution output.

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

       --no_update
          skip the cvs and git updates and run the suite on the code as it
          exists now.

       --single_test mytest
          run only the test named mytest

       --note \"note\"
          print the note on the resulting test webpages
          

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
    if len(sys.argv) == 1:
        print usage
        sys.exit(2)

    try:
        opts, next = getopt.getopt(argv[1:], "",
                                   ["make_benchmarks=",
                                    "no_update",
                                    "single_test=",
                                    "note="])

    except getopt.GetoptError:
        print "invalid calling sequence"
        print usage
        sys.exit(2)


    # defaults
    make_benchmarks = 0
    no_update = 0
    single_test = ""
    comment = ""
    note = ""
    
    for o, a in opts:

        if o == "--make_benchmarks":
            make_benchmarks = 1
            comment = a

        if o == "--no_update":
            no_update = 1

        if o == "--single_test":
            single_test = a
            
        if o == "--note":
            note = a
            
    try:
        testFile = next[0]

    except IndexError:
        print "ERROR: a test file was not specified"
        print usage
        sys.exit(2)

        
    #--------------------------------------------------------------------------
    # read in the test information
    #--------------------------------------------------------------------------
    bold("loading" + testFile)

    suite, testList = LoadParams(testFile)

    if (len(testList) == 0):
        fail("No valid tests defined")

    # store the full path to the testFile
    testFilePath = os.getcwd() + '/' + testFile


    #--------------------------------------------------------------------------
    # get the name of the benchmarks directory
    #--------------------------------------------------------------------------
    benchDir = suite.testTopDir + suite.suiteName + "-benchmarks/"
    if (not os.path.isdir(benchDir)):
        if (make_benchmarks):
            os.mkdir(benchDir)
        else:
            fail("ERROR: benchmark directory, %s, does not exist" % (benchDir))

    

    #--------------------------------------------------------------------------
    # if we are doing a single test, remove all other tests
    #--------------------------------------------------------------------------
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
    webTopDir = "%s/%s-web/" % (suite.testTopDir, suite.suiteName)
    fullWebDir = "%s/%s-web/%s/"  % (suite.testTopDir, suite.suiteName, testDir)
        

    if (not (os.path.isdir(webTopDir)) ):
        os.mkdir(webTopDir)

    os.mkdir(fullWebDir)

    # copy the test file into the web output directory
    shutil.copy(testFilePath, fullWebDir)


    #--------------------------------------------------------------------------
    # do the CVS updates
    #--------------------------------------------------------------------------
    now = time.localtime(time.time())
    cvsTime = time.strftime("%Y-%m-%d %H:%M:%S %Z", now)

    os.chdir(suite.testTopDir)
    
    if (not no_update):

       # Parallel
       if (suite.sourceTree == "Parallel"):
          doCVSUpdate(suite.sourceDir, "Parallel", fullWebDir)

       # fParallel
       doCVSUpdate(suite.sourceDir, "fParallel", fullWebDir)
    
       # BoxLib
       doGITUpdate(suite.boxLibDir, "BoxLib", fullWebDir)


    #--------------------------------------------------------------------------
    # generate the ChangeLogs
    #--------------------------------------------------------------------------
    if (not no_update):

       # Parallel
       if (suite.sourceTree == "Parallel"):
          makeCVSChangeLog(suite, "Parallel", fullWebDir)

       # fParallel
       makeCVSChangeLog(suite, "fParallel", fullWebDir)
    
       # BoxLib
       makeGITChangeLog(suite, "BoxLib", fullWebDir)


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

    bold("building the visualization tools...")

    compString = "%s -j%s BOXLIB_HOME=%s programs=fsnapshot2d NDEBUG=t MPI= COMP=%s  2>&1 > fsnapshot2d.make.out" % \
                   (suite.MAKE, suite.numMakeJobs, suite.boxLibDir, suite.FCOMP)
    print "  " + compString
    systemCall(compString)
    vis2dExecutable = getRecentFileName(suite.compareToolDir,"fsnapshot2d",".exe")
    
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

        print "  %s" % (dir)
        os.chdir(suite.sourceDir + dir)

        if (suite.sourceTree == "Parallel"):
            systemCall("%s BOXLIB_HOME=%s realclean >& /dev/null" % (suite.MAKE, suite.boxLibDir))
        else:
            systemCall("%s BOXLIB_HOME=%s realclean >& /dev/null" % (suite.MAKE, suite.boxLibDir))
            
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
        os.chdir(suite.sourceDir + test.buildDir)
        
        print "  building..."

        if (suite.sourceTree == "Parallel"):

	    if (test.useMPI):
	       executable = "%s%dd.MPI.ex" % (suite.suiteName, test.dim)
               compString = "%s -j%s BOXLIB_HOME=%s DIM=%d USE_MPI=TRUE COMP=%s FCOMP=%s executable=%s  >& %s/%s.make.out" % \
                   (suite.MAKE, suite.numMakeJobs, suite.boxLibDir, test.dim, suite.COMP, suite.FCOMP, executable, outputDir, test.name)
               print "    " + compString
               systemCall(compString)

            else:
	       executable = "%s%dd.ex" % (suite.suiteName, test.dim)
               compString = "%s -j%s BOXLIB_HOME=%s DIM=%d USE_MPI=false COMP=%s FCOMP=%s executable=%s  >& %s/%s.make.out" % \
                   (suite.MAKE, suite.numMakeJobs, suite.boxLibDir, test.dim, suite.COMP, suite.FCOMP, executable, outputDir, test.name)
               print "    " + compString
               systemCall(compString)
	       
            
        elif (suite.sourceTree == "fParallel"):

            if (test.useMPI):
                compString = "%s -j%s BOXLIB_HOME=%s MPI=t NDEBUG=t COMP=%s >& %s/%s.make.out" % \
                    (suite.MAKE, suite.numMakeJobs, suite.boxLibDir, suite.FCOMP, outputDir, test.name)
                print "    " + compString
                systemCall(compString)

            else:
                compString = "%s -j%s BOXLIB_HOME=%s MPI= NDEBUG=t COMP=%s >& %s/%s.make.out" % \
                    (suite.MAKE, suite.numMakeJobs, suite.boxLibDir, suite.FCOMP, outputDir, test.name)
                print "    " + compString
                systemCall(compString)


            # we need a better way to get the executable name here
            executable = getRecentFileName(suite.sourceDir + test.buildDir,"main",".exe")

        

        if (test.compileTest):
            
            # compilation tests are done now -- just make the report and ...
            shutil.copy("%s/%s.make.out"    % (outputDir, test.name), fullWebDir)

            print "  creating problem test report ..."
            reportSingleTest(suite, test, testDir, fullWebDir)

            # ... skip to the next test in the loop
            continue
            
            
        #----------------------------------------------------------------------
        # copy the necessary files over to the run directory
        #----------------------------------------------------------------------
        print "  copying files to run directory..."

        try: shutil.copy(executable, outputDir)
        except IOError:
           errorMsg = "    ERROR: compilation failed"
           reportTestFailure(errorMsg, test, testDir, fullWebDir)
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


        # if we are a "Parallel" build, we need the probin file
        if (suite.sourceTree == "Parallel"):
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

        if (suite.sourceTree == "Parallel"):

	    if (test.useMPI):

	       # create the MPI executable
	       testRunCommand = suite.MPIcommand
	       testRunCommand = testRunCommand.replace("@host@", suite.MPIhost)
	       testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (test.numprocs))

	       command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk >&  %s.run.out < /dev/null" % \
			 (executable, test.inputFile, test.name, test.name, test.name)
	       
	       testRunCommand = testRunCommand.replace("@command@", command)

	       print "    " + testRunCommand
               systemCall(testRunCommand)
	       
            else:
               testRunCommand = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk >&  %s.run.out" % \
                   (executable, test.inputFile, test.name, test.name, test.name)

               print "    " + testRunCommand
               systemCall(testRunCommand)

        elif (suite.sourceTree == "fParallel"):

            if (test.useMPI):

	       # create the MPI executable
	       testRunCommand = suite.MPIcommand
	       testRunCommand = testRunCommand.replace("@host@", suite.MPIhost)
	       testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (test.numprocs))

                # keep around the checkpoint files only for the restart runs
	       if (test.restartTest):
                   command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk >& %s.run.out" % \
                       (executable, test.inputFile, test.name, test.name, test.name)
               else:
                   command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 >& %s.run.out" % \
                       (executable, test.inputFile, test.name, test.name, test.name)

	       testRunCommand = testRunCommand.replace("@command@", command)

	       print "    " + testRunCommand
               systemCall(testRunCommand)

            else:

                # keep around the checkpoint files only for the restart runs
                if (test.restartTest):
                    testRunCommand = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk >& %s.run.out" % \
                        (executable, test.inputFile, test.name, test.name, test.name)
                else:
                    testRunCommand = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 >& %s.run.out" % \
                        (executable, test.inputFile, test.name, test.name, test.name)

                print "    " + testRunCommand
                systemCall(testRunCommand)


        # if it is a restart test, then rename the final output file and
        # restart the test

        if (test.restartTest):
            lastFile = getLastPlotfile(outputDir, test)
            origLastFile = "orig_%s" % (lastFile)
            shutil.move(lastFile, origLastFile)

            # get the file number to restart from
            restartFile = "%s_chk%5.5d" % (test.name, test.restartFileNum)

            print "  restarting from %s ... " % (restartFile)
           
            if (suite.sourceTree == "Parallel"):
                testRunCommand = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.restart=%s >>  %s.run.out 2>&1" % \
                    (executable, test.inputFile, test.name, test.name, restartFile, test.name)
                print "    " + testRunCommand
                systemCall(testRunCommand)               

            elif (suite.sourceTree == "fParallel"):
                testRunCommand = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --restart %d >> %s.run.out 2>&1" % \
                    (executable, test.inputFile, test.name, test.name, test.restartFileNum, test.name)
                print "    " + testRunCommand
                systemCall(testRunCommand)                
           
            
        #----------------------------------------------------------------------
        # do the comparison
        #----------------------------------------------------------------------
        if (not test.selfTest):

            compareFile = getLastPlotfile(outputDir, test)

            if (not make_benchmarks):

                print "  doing the comparison..."
                print "    comparison file: ", compareFile

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
                   
                        command = "../fcompare.exe -n 0 --infile1 %s --infile2 %s >> %s.compare.out 2>&1" % (benchFile, compareFile, test.name)

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

            else:   # make_benchmarks

                print "  storing output of %s as the new benchmark..." % (test.name)
                print "     new benchmark file: ", compareFile

                if (not compareFile == ""):
                    systemCall("cp -rf %s %s" % (compareFile, benchDir))

                    cf = open("%s.status" % (test.name), 'w')
                    cf.write("benchmarks updated.  New file:  %s\n" % (compareFile) )
                    cf.close()
                else:
                    cf = open("%s.status" % (test.name), 'w')
                    cf.write("benchmarks failed")
                    cf.close()

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

            if (not compareFile == ""):

                print "  doing the visualization..."

                if (test.dim == 2):
                    systemCall('%s/%s --palette %s/Palette -cname "%s" -p "%s" >& /dev/null' %
                              (suite.compareToolDir, vis2dExecutable, suite.compareToolDir, 
                               test.visVar, compareFile) )
                elif (test.dim == 3):
                    systemCall('%s/%s --palette %s/Palette -n 1 -cname "%s" -p "%s" >& /dev/null' %
                              (suite.compareToolDir, vis3dExecutable, suite.compareToolDir, 
                               test.visVar, compareFile) )
                else:
                    print "    Visualization not supported for dim = %d" % (test.dim)
                    

                # convert the .ppm files into .png files
                ppmFile = getRecentFileName(outputDir, "", ".ppm")

                systemCall("convert %s `basename %s .ppm`.png" % 
                          (ppmFile, ppmFile) )
        
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

            if (suite.sourceTree == "Parallel"):
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
            reportSingleTest(suite, test, testDir, fullWebDir)


        
    #--------------------------------------------------------------------------
    # write the report for this instance of the test suite
    #--------------------------------------------------------------------------
    print "\n"
    bold("creating new test report...")
    reportThisTestRun(suite, make_benchmarks, comment, note, cvsTime, 
                      testList, testDir, testFile, fullWebDir)


    # make sure that all of the files in the web directory are world readable
    for file in os.listdir(fullWebDir):    
       currentFile = fullWebDir + file

       if (os.path.isfile(currentFile)):
          os.chmod(currentFile, 0644)
          

    #--------------------------------------------------------------------------
    # generate the master report for all test instances 
    #--------------------------------------------------------------------------
    print "\n"
    bold("creating suite report...")
    reportAllRuns(suite, webTopDir)



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# R E P O R T   W R I T I N G   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

cssContents = \
r"""
h3.passed {text-decoration: none; display: inline; 
           color: black; background-color: lime; padding: 2px}
h3.failed {text-decoration: none; display: inline; 
           color: black; background-color: red; padding: 2px}
h3.benchmade {text-decoration: none; display: inline; 
              color: black; background-color: orange; padding: 2px}
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
<CENTER><H1>@TITLE@</H1></CENTER>
"""


#==============================================================================
# create_css
#==============================================================================
def create_css():
    """ write the css file for the webpages """
    
    cssFile = "tests.css"
    cf = open(cssFile, 'w')
    cf.write(cssContents)
    cf.close()


#==============================================================================
# reportSingleTest
#==============================================================================
def reportSingleTest(suite, test, testDir, fullWebDir):
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

    # check to see if the CSS file is present, if not, write it
    if (not os.path.isfile("tests.css")):
        create_css()

    
    htmlFile = "%s.html" % (test.name)
    hf = open(htmlFile, 'w')

    newHead = HTMLHeader + r"""<CENTER><H1><A HREF="index.html">@TESTDIR@</A> / @TESTNAME@</H1></CENTER>"""

    newHead = newHead.replace("@TESTDIR@", testDir)
    newHead = newHead.replace("@TESTNAME@", test.name)

    hf.write(newHead)


    if (not test.compileTest):

        if (test.useMPI):

            hf.write("<P><b>Parallel Run</b><br>numprocs = %d\n" % (test.numprocs) )
            hf.write("<P>&nbsp;\n")
       

        # is this a restart test?
        if (test.restartTest):

            hf.write("<P><b>Restart Test</b><br>Job was run as normal and then restarted from checkpoint # %d, and the two final outputs were compared\n" % (test.restartFileNum) )

        hf.write("<P>&nbsp;\n")       

        # write out the information about the test
        hf.write("<P><b>build directory:</b> %s\n" % (test.buildDir) )

        hf.write("<P><b>input file:</b> <A HREF=\"%s.%s\">%s</A>\n" %
                 (test.name, test.inputFile, test.inputFile) )

        if (suite.sourceTree == "Parallel"):
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

    hf.write("<A HREF=\"%s.make.out\">make output</A>\n" % (test.name) )

    hf.write("<P>&nbsp;\n")

    
    if (not test.compileTest):

        # write out the comparison report
        if (compareSuccessful):
            hf.write("<P><H3 CLASS=\"passed\">Comparison Successful</H3></P>\n")
        else:
            hf.write("<P><H3 CLASS=\"failed\">Comparison Failed</H3></P>\n")

        hf.write("<A HREF=\"%s.run.out\">execution output</A>\n" % (test.name) )


        hf.write("<P>&nbsp;\n")
        hf.write("<PRE>\n")
    
        for line in diffLines:
            hf.write(line)

        hf.write("</PRE>\n")

        # show any visualizations
        if (test.doVis):
            pngFile = getRecentFileName(fullWebDir, test.name, ".png")
            hf.write("<P>&nbsp;\n")
            hf.write("<P><IMG SRC='%s' BORDER=0>" % (pngFile) )
    

    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")    

    hf.close()
    
    
    # switch back to the original directory
    os.chdir(currentDir)
	

#==============================================================================
# reportTestAbort
#==============================================================================
def reportTestFailure(message, test, testDir, fullWebDir):
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

    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")    

    hf.close()
    
    
    # switch back to the original directory
    os.chdir(currentDir)
	

#==============================================================================
# reportThisTestRun
#==============================================================================
def reportThisTestRun(suite, make_benchmarks, comment, note, cvsTime, 
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

    # check to see if the CSS file is present, if not, write it
    if (not os.path.isfile("tests.css")):
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

    hf.write("<p>&nbsp;\n")
    hf.write("<p><b>CVS update was done at: </b>%s\n" % (cvsTime) )

    if (suite.sourceTree == "Parallel"):
        hf.write("<p>&nbsp;&nbsp;<b>cvs update on Parallel/:</b> <A HREF=\"%s\">%s</A>\n" %
                 ("cvs.Parallel.out", "cvs.Parallel.out") )

    hf.write("<p>&nbsp;&nbsp;<b>cvs update on fParallel/:</b> <A HREF=\"%s\">%s</A>\n" %
             ("cvs.fParallel.out", "cvs.fParallel.out") )        
    hf.write("<p>&nbsp;\n")


    if (suite.sourceTree == "Parallel"):
            hf.write("<p>&nbsp;&nbsp;<b>Parallel/ ChangeLog:</b> <A HREF=\"%s\">%s</A>\n" %
                     ("ChangeLog.Parallel", "ChangeLog.Parallel") )
            
    hf.write("<p>&nbsp;&nbsp;<b>fParallel/ ChangeLog:</b> <A HREF=\"%s\">%s</A>\n" %
             ("ChangeLog.fParallel", "ChangeLog.fParallel") )        
    hf.write("<p>&nbsp;\n")    

    hf.write("<P><TABLE BORDER=0 CELLPADDING=3>\n")
    
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
    statusFile = testDir[0:index] + ".status"

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

	
#==============================================================================
# reportAllRuns
#==============================================================================
def reportAllRuns(suite, webTopDir):

    os.chdir(webTopDir)

    if (not os.path.isfile("tests.css")):
        create_css()

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
    hf.write(header)

    hf.write("<P><TABLE BORDER=0 CELLPADDING=5>\n")

    # write out the header
    hf.write("<TR><TD ALIGN=CENTER>date</TD><TD>&nbsp;</TD>")
    for test in allTests:
        hf.write("<TD ALIGN=CENTER>%s</TD>" % (test))
    
    hf.write("</TR>\n")


    # loop over all the test runs 
    for dir in validDirs:

        hf.write("<TR><TD><A HREF=\"%s/index.html\">%s</A></TD><TD>&nbsp;</TD>" %
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
                hf.write("<TD ALIGN=CENTER><H3 class=\"passed\">PASSED</H3></TD>")
            elif (status == -1):
                hf.write("<TD ALIGN=CENTER><H3 class=\"failed\">FAILED</H3></TD>")
            elif (status == 10):
                hf.write("<TD ALIGN=CENTER><H3 class=\"benchmade\">new benchmark</H3></TD>")
            else:
                hf.write("<TD>&nbsp;</TD>")


        hf.write("</TR>\n")
        
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
