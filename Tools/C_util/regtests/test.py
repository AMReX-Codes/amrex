#!/usr/bin/env python

"""
A simple regression test framework for a BoxLib-based code

There are several major sections to this source: the runtime parameter
routines, the test suite routines, and the report generation routines.  They
are separated as such in this file.

This test framework understands source based out of the Parallel/ and
fParallel/ frameworks.

2007-11-14
"""

import os
import shutil
import sys
import getopt
import datetime
import time
import string



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# R U N T I M E   P A R A M E T E R   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import ConfigParser
import string, sys, re

# we will keep track of the parameters globally
globalParams = {}
globalSections = []


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


#==============================================================================
# setParamDefaults
#==============================================================================
def setParamDefaults(defaults):
    global globalParams
    
    globalParams = defaults.copy()


#==============================================================================
# LoadParams
#==============================================================================
def LoadParams(file):
    """
    returns a dictionary with keys of the form
    <section>.<option> and the corresponding values
    """
    global globalParams
    global globalSections

    # check to see whether the file exists
    try: f = open(file, 'r')
    except IOError:
        print 'ERROR: parameter file does not exist: ', file
        sys.exit()
    else:
        f.close()


    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.read(file)

    globalSections = cp.sections()

    for sec in cp.sections():

        for opt in cp.options(sec):

            value = cp.get(sec, opt)
            
            # check in turn whether this is an interger, float, or string
            if (isInt(value)):
                globalParams[sec + "." + opt] = int(value)
            elif (isFloat(value)):
                globalParams[sec + "." + opt] = float(value)
            else:
                globalParams[sec + "." + opt] = value.strip()


#==============================================================================
# getParam
#==============================================================================
def getParam(key):
    """ return the value of the runtime parameter corresponding to key """
    
    if globalParams == {}:
        print "WARNING: runtime parameters not yet initialized"
        LoadParams("_defaults")
        
    if key in globalParams.keys():
        return globalParams[key]
    else:
        raise ValueError()
        

#==============================================================================
# keyIsValid
#==============================================================================
def keyIsValid(key):
    """ check whether the key is contained in the global parameter list """

    isValid = 1
    
    try:
        temp = getParam(key)

    except ValueError:
        isValid = 0
        print "   WARNING: %s not set" % (key)

    return isValid


#==============================================================================
# getSections
#==============================================================================
def getSections():
    """ return the sections """
    sections = globalSections
    sections.sort()
    return sections


#==============================================================================
# PrintAllParams
#==============================================================================
def PrintAllParams():
    keys = globalParams.keys()
    keys.sort()

    for key in keys:
        print key, "=", globalParams[key]

    print " "
    
            

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# T E S T   S U I T E   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#==============================================================================
# abortTests
#==============================================================================
def abortTests(message):
    print message
    sys.exit(2)

    
#==============================================================================
# findBuildDirs
#==============================================================================
def findBuildDirs(tests):
    """ given the list of tests, find the set of UNIQUE build directories"""
    
    buildDirs = []
    for test in tests:

        currentBuildDir = getParam(test + ".buildDir")
        if (buildDirs.count(currentBuildDir) == 0):
            buildDirs.append(currentBuildDir)

    return buildDirs


#==============================================================================
# getValidTests
#==============================================================================
def getValidTests(sourceTree):
    """ get the list of tests and make sure that the tests have all the
        required parameters.  The returned list contains only those that
        have valid data """

    tests = getSections()
    newTests = tests[:]

    # [main] is reserved for test suite parameters
    newTests.remove("main")

    removeList = []
    
    for test in newTests:

        print " "
        print "checking parameters for test %s" % (test)
        
        # check for the manditory parameters
        if (not (keyIsValid("%s.buildDir" % (test))      and
                 keyIsValid("%s.inputFile" % (test))     and
                 (sourceTree == "fParallel" or 
                  keyIsValid("%s.probinFile" % (test)) ) and
                 keyIsValid("%s.needs_helmeos" % (test)) and
                 keyIsValid("%s.dim" % (test)) ) ):
            print "ERROR: manditory runtime parameters for test %s not set\n" % (test)
            removeList.append(test)
            continue


        # check for optional parameters

        # restartTest
        if (not keyIsValid("%s.restartTest" % (test)) ):
            print "   Assuming normal (not restart) run.\n" 
            globalParams["%s.restartTest" % (test)] = 0
        else:

           if (getParam("%s.restartTest" % (test)) ):

              # make sure that the file number to restart from is defined
              if (not keyIsValid("%s.restartFileNum" % (test)) ):
                 print "ERROR: test %s is a restart test, but is missing the restartFileNum parameter.\n" % (test)
                 removeList.append(test)
                 continue
      

        # compileTest
        if (not keyIsValid("%s.compileTest" % (test)) ):
            print "   Assuming normal (not compilation test) run.\n" 
            globalParams["%s.compileTest" % (test)] = 0


        # useMPI
        if (not keyIsValid("%s.useMPI" % (test)) ):
            print "   Assuming normal (not MPI) run.\n"
            globalParams["%s.useMPI" % (test)] = 0
        else:

           if (getParam("%s.useMPI" % (test)) ):

              # make sure that the number of processors is defined
              if (not keyIsValid("%s.numprocs" % (test)) ):
                 print "ERROR: test %s is a parallel test, but did not specify the numprocs parameter.\n" % (test)
                 removeList.append(test)
                 continue
      

        # doVis
        if (not keyIsValid("%s.doVis" % (test)) ):
            print "   Assuming no visualization.\n"
            globalParams["%s.doVis" % (test)] = 0
        else:

           if (getParam("%s.doVis" % (test)) ):

              # find out what variable to plot
              if (not keyIsValid("%s.visVar" % (test)) ):
                 print "ERROR: test %s requested visualization but did not specify the visVar parameter.\n" % (test)
                 removeList.append(test)
                 continue
      

    # remove the invalid tests
    for test in removeList:
       newTests.remove(test)
       
    return newTests


#==============================================================================
# getAuxFiles
#==============================================================================
def getAuxFiles(test):
    """ given a test, get all the list of all the auxillary files it needs """
    
    auxFiles = []

    MAX_AUX_FILES = 9
    i = 1
    while (i < MAX_AUX_FILES):

        # the auxillary files will have names like aux1File, aux2File, ...
        try:
            aux = getParam("%s.aux%dFile" % (test, i) )

        except ValueError:
            i = i + 1
            continue

        auxFiles.append(aux)
        i = i + 1
        
    return auxFiles


#==============================================================================
# getLastPlotfile
#==============================================================================
def getLastPlotfile(outputDir, test):
    """ given an output directory and the test name, find the last
        plotfile written """
        
    plotNum = -1
        
    # start by finding the last plotfile
    for file in os.listdir(outputDir):
       if (os.path.isdir(file) and file.startswith("%s_plt" % (test))):
           key = "_plt"
           index = string.rfind(file, key)
           plotNum = max(int(file[index+len(key):]), plotNum)

    if (plotNum == -1):
       print "WARNING: test did not produce any output"
       compareFile = ""
    else:
       compareFile = "%s_plt%5.5d" % (test, plotNum)

    return compareFile


#==============================================================================
# getRecentFileName
#==============================================================================
def getRecentFileName(dir,base,extension):
    """ for fParallel builds, given the base and extension, find the most recent
        corresponding file """
        
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
def getTestCompDir(dirKey):
   """ given a key, dirKey, from the inputs file, check (1) if the key is 
       valid, and (2) if it points to a vaid directory.  If so, return the
       directory name """

   if (not keyIsValid(dirKey)):
       abortTests("ERROR: %s is not defined" % (dirKey))

   else:
       dirName = getParam(dirKey)

       # make sure we end in a "/"
       if (not (string.rfind(dirName, "/") == len(dirName)-1)):
           dirName = dirName + "/"
           
       if (not os.path.isdir(dirName)):
           abortTests("ERROR: %s is not a valid directory" % (dirName))

       return dirName


#==============================================================================
# doCVSUpdate
#==============================================================================
def doCVSUpdate(topDir, root, outDir):
   """ do a CVS update of the repository named root.  topDir is the full path
       to the directory containing root.  outDir is the full path to the 
       directory where we will store the CVS output """

   os.chdir(topDir)
 
   print "cvs update %s" % (root)
   os.system("cvs update %s >& cvs.%s.out" % (root, root))
   print "\n"
        
   # make sure that the cvs update was not aborted -- this can happen, for
   # instance, if the CVSROOT was not defined
   try:
       cf = open("cvs.%s.out" % (root), 'r')

   except IOError:
       abortTests("ERROR: no CVS output")

   else:
       cvsLines = cf.readlines()

       cvsFailed = 0
    
       for line in cvsLines:
           if (string.find(line, "update aborted") >= 0):
               cvsFailed = 1
               break
        
       cf.close()

   if (cvsFailed):
       abortTests("ERROR: cvs update was aborted. See cvs.%s.out for details" % (root) )
        
   shutil.copy("cvs.%s.out" % (root),  outDir)


#==============================================================================
# makeChangeLog
#==============================================================================
def makeChangeLog(topDir, root, outDir):
   """ generate a ChangeLog for the CVS repository named root.  topDir is 
       the full path to the directory containing root.  outDir is the full 
       path to the directory where we will store the CVS output """

   os.chdir(topDir)

   have_cvs2cl = 0 

   print "Generating ChangeLog for %s/" % (root)
    
   if (not os.path.isfile("%s/cvs2cl.pl" % (root) )):
       if (not os.path.isfile("fParallel/scripts/cvs2cl.pl")):
           print "WARNING: unable to locate cvs2cl.pl script."
           print "         no ChangeLog will be generated"
       else:
           shutil.copy("fParallel/scripts/cvs2cl.pl", "%s/" % (root) )
           have_cvs2cl = 1
   else:
       have_cvs2cl = 1

   os.chdir("%s/" % (root) )

   if (have_cvs2cl):
       os.system("./cvs2cl.pl -f ChangeLog.%s" % (root) )
   else:
       cf = open("ChangeLog.%s" % (root), 'w')
       cf.write("unable to generate ChangeLog")
       cf.close()

   os.chdir(topDir)

   shutil.copy("%s/ChangeLog.%s" % (root, root), outDir)



#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# test
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def test(argv):

    usage = """
    ./test.py [--make_benchmarks comment,
               --no_cvs_update,
               --single_test test,
               --note note]
        testfile.ini


    arguments:

      testfile.ini 
          This is the input file that defines the tests that the
          suite knows about.  It has the format

            [main]
            sourceDir      = < directory above Parallel/ and fParallel/ >
            testTopDir     = < full path to test output directory >
            compareToolDir = < full path to the fParallel/data_processing/ directory >
            helmeosDir     = < full path to helm_table.dat >

            sourceTree = < Parallel or fParallel -- what type is it? >

            suiteName = < descriptive name (i.e. Castro) >

            MPIcommand = < MPI run command, with placeholders for host,
                           # of proc, and program command.  Should look
                           something like :
                           mpiexec -host @host@ -n @nprocs@ @command@ >
                           
            MPIhost = < host for MPI job -- depends on MPI implementation >

            

            [Sod-x]
            buildDir = < relative path (from sourceDir) for this problem >
            inputFile = < input file name >
            probinFile = < probin file name >
            needs_helmeos = < need Helmholtz eos? 0 for no, 1 for yes >
            restartTest = < is this a restart test? 0 for no, 1 for yes >
            restartFileNum = < # of file to restart from (if restart test) >
            dim = < dimensionality: 1, 2, or 3 >
            useMPI = <is this a parallel job? 0 for no, 1 for yes) >
            numprocs = < # of processors to run on (if parallel job) >
            compileTest = < 0 for normal run, 1 if we just test compilation >            

          Here, [main] lists the parameters for the test suite as a
          whole and [Sod-x] is a single test.  There can be many more
          tests, each with their own unique name, following the format
          of [Sod-x] above.

          In [main],

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

            needs_helmeos is set to 1 if the Helmholtz EOS is used.
            This will ensure that the helm_table.dat file is copied
            into the run directory.
            
            dim is the dimensionality for the problem.

            restartTest = 1 means that this is a restart test.  Instead of
            comparing to a stored benchmark, we will run the test and then
            restart from restartFileNum and run until the end.  The last
            file from the original run and the last from the restart will
            be compared.

            useMPI is set to 1 if the job is parallel.  In this case, you
            also must specify the number of processors, via numprocs.

            compileTest is set to 0 for a normal test.  Setting to 1 will
            just test compilation, and not any subsequent output.

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

       --no_cvs_update
          skip the cvs update and run the suite on the code as it
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
                                    "no_cvs_update",
                                    "single_test=",
                                    "note="])

    except getopt.GetoptError:
        print "invalid calling sequence"
        print usage
        sys.exit(2)


    # defaults
    make_benchmarks = 0
    no_cvs_update = 0
    single_test = ""
    comment = ""
    note = ""
    
    for o, a in opts:

        if o == "--make_benchmarks":
            make_benchmarks = 1
            comment = a

        if o == "--no_cvs_update":
            no_cvs_update = 1

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
    
    LoadParams(testFile)

    # store the full path to the testFile
    testFilePath = os.getcwd() + '/' + testFile


    #--------------------------------------------------------------------------
    # get the test directory parameters
    #--------------------------------------------------------------------------

    # get the source directory
    sourceDir = getTestCompDir("main.sourceDir")
        
    # get the top-level testing directory
    testTopDir = getTestCompDir("main.testTopDir")

    # get the directory that contains the comparison tools
    compareToolDir = getTestCompDir("main.compareToolDir")

    # get the directory that contains the Helmholtz EOS
    helmeosDir = getTestCompDir("main.helmeosDir")


    # get the name of the test suite
    if (not keyIsValid("main.suiteName")):
        suiteName = "testdefault"

    else:
        suiteName = getParam("main.suiteName")


    # get the type of source we are dealing with: Parallel or fParallel
    # these have different Makefile styles, so we'll need to do things
    # slightly differently for each.
    if (not keyIsValid("main.sourceTree")):
        print "WARNING: sourceTree not set, assuming Parallel"
        sourceTree = "Parallel"

    else:
        sourceTree = getParam("main.sourceTree")

        if (not (sourceTree == "Parallel" or sourceTree == "fParallel")):
           abortTests("ERROR: invalid sourceTree")
           

    # get the name of the benchmarks directory
    benchDir = testTopDir + suiteName + "-benchmarks/"
    if (not os.path.isdir(benchDir)):
        if (make_benchmarks):
            os.mkdir(benchDir)
        else:
            abortTests("ERROR: benchmark directory, %s, does not exist" % (benchDir))

    
    #--------------------------------------------------------------------------
    # get the general MPI parameters (if available)
    #--------------------------------------------------------------------------
    if (not keyIsValid("main.MPIcommand")):
       MPIcommand = ""

    else:
       MPIcommand = getParam("main.MPIcommand")


    if (not keyIsValid("main.MPIhost")):
       MPIhost = ""

    else:
       MPIhost = getParam("main.MPIhost")


    #--------------------------------------------------------------------------
    # get the Fortran compiler for fParallel tests and utilities (if available)
    #--------------------------------------------------------------------------
    if (not keyIsValid("main.FCOMP")):

       print "WARNING: FCOMP not set; assuming Intel"
       FCOMP = "Intel"

    else:
       FCOMP = getParam("main.FCOMP")

       
    #--------------------------------------------------------------------------
    # get the name of the individual tests 
    #--------------------------------------------------------------------------
    tests = getValidTests(sourceTree)

    if (not single_test == ""):
        if (single_test in tests):
            print "running only test: %s " % (single_test)
            tests = [single_test]
        else:
            abortTests("ERROR: %s is not a valid test" % (single_test))
        
    print tests
    
    
    #--------------------------------------------------------------------------
    # create the output directories
    #--------------------------------------------------------------------------
    os.chdir(testTopDir)

    todayDate = datetime.date.today()
    today = todayDate.__str__()


    # figure out what the current output directory should be
    maxRuns = 100      # maximum number of tests in a given day

    testDir = today + "/"

    # test output stored in a directory suiteName-tests/2007-XX-XX/
    # make sure that the suiteName-tests directory exists
    if (not os.path.isdir(testTopDir + suiteName + "-tests/")):
        os.mkdir(testTopDir + suiteName + "-tests/")
        
    fullTestDir = testTopDir + suiteName + "-tests/" + testDir

    i = 0
    while (i < maxRuns-1 and os.path.isdir(fullTestDir)):
        i = i + 1
        testDir = today + "-" + ("%3.3d" % i) + "/"
        fullTestDir = testTopDir + suiteName + "-tests/" + testDir

    print "testing directory is: ", testDir
    os.mkdir(fullTestDir)


    # make the web directory -- this is where all the output and HTML will be
    # put, so it is easy to move the entire test website to a different disk
    webTopDir = "%s/%s-web/" % (testTopDir, suiteName)
    fullWebDir = "%s/%s-web/%s/"  % (testTopDir, suiteName, testDir)
        

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

    os.chdir(testTopDir)
    
    if (not no_cvs_update):

       # Parallel
       if (sourceTree == "Parallel"):
          doCVSUpdate(testTopDir, "Parallel", fullWebDir)

       # fParallel
       doCVSUpdate(testTopDir, "fParallel", fullWebDir)
    

    #--------------------------------------------------------------------------
    # generate the ChangeLogs
    #--------------------------------------------------------------------------
    if (not no_cvs_update):

       # Parallel
       if (sourceTree == "Parallel"):
          makeChangeLog(testTopDir, "Parallel", fullWebDir)

       # fParallel
       makeChangeLog(testTopDir, "fParallel", fullWebDir)
    

    #--------------------------------------------------------------------------
    # build the comparison and visualization tools
    #--------------------------------------------------------------------------
    print "building the comparison tools..."

    os.chdir(compareToolDir)

    os.system("gmake COMP=%s realclean" % (FCOMP))

    os.system("gmake -j 4 programs=fcompare NDEBUG=t MPI= COMP=%s  >& fcompare.make.out" % (FCOMP))
    compareExecutable = getRecentFileName(compareToolDir,"fcompare",".exe")

    shutil.copy(compareExecutable, fullTestDir + "/fcompare.exe")

    print "building the visualization tools..."
    os.system("gmake -j 4 programs=fsnapshot2d NDEBUG=t MPI= COMP=%s  >& fsnapshot2d.make.out" % (FCOMP))
    vis2dExecutable = getRecentFileName(compareToolDir,"fsnapshot2d",".exe")
    
    os.system("gmake -j 4 programs=fsnapshot3d NDEBUG=t MPI= COMP=%s  >& fsnapshot3d.make.out" % (FCOMP))
    vis3dExecutable = getRecentFileName(compareToolDir,"fsnapshot3d",".exe")
    
    print "\n"


    #--------------------------------------------------------------------------
    # do a make clean, only once per build directory
    #--------------------------------------------------------------------------
    allBuildDirs = findBuildDirs(tests)

    print "make clean in ..."

    for dir in allBuildDirs:

        print "  %s..." % (dir)
        os.chdir(sourceDir + dir)

        if (sourceTree == "Parallel"):
	   os.system("gmake realclean")
        else:
            os.system("gmake MPI= COMP=%s realclean NDEBUG=t" % (FCOMP))
            
    print "\n"

    os.chdir(testTopDir)
    

    #--------------------------------------------------------------------------
    # main loop over tests
    #--------------------------------------------------------------------------
    for test in tests:

        print "working on " + test + " test"

        if (getParam(test + ".restartTest") and make_benchmarks):
            print "  WARNING: test %s is a restart test -- " % (test)
            print "           no benchmarks are stored." 
            print "           skipping\n"
            continue
     
        if (getParam(test + ".compileTest") and make_benchmarks):
            print "  WARNING: test %s is a compile test -- " % (test)
            print "           no benchmarks are stored."
            print "           skipping\n"
            continue            


        #----------------------------------------------------------------------
        # make the run directory
        #----------------------------------------------------------------------
        outputDir = fullTestDir + test + '/'
        os.mkdir(outputDir)
    

        #----------------------------------------------------------------------
        # compile the code
        #----------------------------------------------------------------------
        buildDir = getParam(test + ".buildDir")
        os.chdir(sourceDir + buildDir)
        
        print "  building..."

        if (sourceTree == "Parallel"):

            dim = getParam(test + ".dim")

	    useMPI = getParam(test + ".useMPI")
	    
	    if (useMPI):
	       executable = "%s%dd.MPI.ex" % (suiteName, dim)
	       os.system("gmake -j 4 DIM=%d USE_MPI=TRUE executable=%s  >& %s/%s.make.out" % 
			 (dim, executable, outputDir, test))
            else:
	       executable = "%s%dd.ex" % (suiteName, dim)
	       os.system("gmake -j 4 DIM=%d USE_MPI=false executable=%s  >& %s/%s.make.out" % 
			 (dim, executable, outputDir, test))
	       
            
        elif (sourceTree == "fParallel"):
            os.system("gmake -j 4 MPI= NDEBUG=t COMP=%s >& %s/%s.make.out" % (FCOMP, outputDir, test))

            # we need a better way to get the executable name here
            executable = getRecentFileName(sourceDir + buildDir,"main",".exe")

        
        compileTest = getParam(test + ".compileTest")

        if (compileTest):
            
            # compilation tests are done now -- just make the report and 
            shutil.copy("%s/%s.make.out"    % (outputDir, test), fullWebDir)

            print "  creating problem test report ..."
            reportSingleTest(sourceTree, test, testDir, fullWebDir)

            # skip to the next test in the loop
            continue
            
            
        #----------------------------------------------------------------------
        # copy the necessary files over to the run directory
        #----------------------------------------------------------------------
        print "  copying files to run directory..."

        try: shutil.copy(executable, outputDir)
        except IOError:
           print '  ERROR: compilation failed'

           statusFile = fullWebDir + "%s.status" % (test)
           sf = open(statusFile, 'w')

           sf.write("FAILED\n")
           sf.close()
           
           print "\n"           
           continue

        inputsFile = getParam(test + ".inputFile")
        shutil.copy(inputsFile, outputDir)

	# sometimes the input file was in a subdirectory under the
	# build directory.  Keep only the input file for latter
	index = string.rfind(inputsFile, "/")
	if (index > 0):
	   inputsFile = inputsFile[index+1:]


        # if we are a "Parallel" build, we need the probin file
        if (sourceTree == "Parallel"):
            probinFile = getParam(test + ".probinFile")
            shutil.copy(probinFile, outputDir)

            # sometimes the probin file was in a subdirectory under the
            # build directory.  Keep only the probin file for latter
            index = string.rfind(probinFile, "/")
            if (index > 0):
               probinFile = probinFile[index+1:]	   


        # if we are using the Helmholtz EOS, we need the input table
        needs_helmeos = getParam(test + ".needs_helmeos")
        if (needs_helmeos):
            helmeosFile = helmeosDir + "helm_table.dat"
            shutil.copy(helmeosFile, outputDir)

        # copy any other auxillary files that are needed at runtime
        auxFiles = getAuxFiles(test)

        for file in auxFiles:
            shutil.copy(file, outputDir)
            

        #----------------------------------------------------------------------
        # run the test
        #---------------------------------------------------------------------- 
        print "  running the test..."

        restart = getParam(test + ".restartTest")

        os.chdir(outputDir)

	useMPI = getParam(test + ".useMPI")
	
        if (sourceTree == "Parallel"):

	    if (useMPI):
	       numprocs = useMPI = getParam(test + ".numprocs")

	       # create the MPI executable
	       testRunCommand = MPIcommand
	       testRunCommand = testRunCommand.replace("@host@", MPIhost)
	       testRunCommand = testRunCommand.replace("@nprocs@", "%s" % (numprocs))

	       command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk >&  %s.run.out < /dev/null" % \
			 (executable, inputsFile, test, test, test)
	       
	       testRunCommand = testRunCommand.replace("@command@", command)

	       print testRunCommand
	       
               os.system(testRunCommand)
	       
            else:
               os.system("./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk >&  %s.run.out" %
			 (executable, inputsFile, test, test, test))	       


        elif (sourceTree == "fParallel"):

            # keep around the checkpoint files only for the restart runs
            if (restart):
                os.system("./%s %s --plot_base_name %s_plt --check_base_name %s_chk >& %s.run.out" %
                          (executable, inputsFile, test, test, test))
            else:
                os.system("./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 >& %s.run.out" %
                          (executable, inputsFile, test, test, test))


        # if it is a restart test, then rename the final output file and
        # restart the test

        if (restart):
           lastFile = getLastPlotfile(outputDir, test)
           origLastFile = "orig_%s" % (lastFile)
           shutil.move(lastFile, origLastFile)

           # get the file number to restart from
           restartFileNum = getParam(test + ".restartFileNum")
	   restartFile = "%s_chk%5.5d" % (test, restartFileNum)

           print "    restarting from %s ... " % (restartFile)
           
           if (sourceTree == "Parallel"):
              os.system("./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.restart=%s >>  %s.run.out 2>&1" %
                      (executable, inputsFile, test, test, restartFile, test))

           elif (sourceTree == "fParallel"):
              os.system("./%s %s --plot_base_name %s_plt --check_base_name %s_chk --restart %d >> %s.run.out 2>&1" %
                      (executable, inputsFile, test, test, restartFileNum, test))
           
            
        #----------------------------------------------------------------------
        # do the comparison
        #----------------------------------------------------------------------
        compareFile = getLastPlotfile(outputDir, test)
        dim = getParam(test + ".dim")

        if (not make_benchmarks):

            print "  doing the comparison..."
            print "    comparison file: ", compareFile

            if (not restart):
               benchFile = benchDir + compareFile
            else:
               benchFile = origLastFile

            # see if it exists
            # note, with BoxLib, the plotfiles are actually directories
            
            if (not os.path.isdir(benchFile)):
                print "WARNING: no corresponding benchmark found"
                benchFile = ""

                cf = open("%s.compare.out" % (test), 'w')
                cf.write("WARNING: no corresponding benchmark found\n")
                cf.write("         unable to do a comparison\n")
                cf.close()
                    
            else:
                if (not compareFile == ""):

                    print "    benchmark file: ", benchFile
                   
                    command = "../fcompare.exe -n 0 --infile1 %s --infile2 %s >> %s.compare.out 2>&1" % (benchFile, compareFile, test)

                    cf = open("%s.compare.out" % (test), 'w')
                    cf.write(command)
                    cf.write("\n")
                    cf.close()
                    
                    os.system(command)

                else:
                    print "WARNING: unable to do a comparison"

                    cf = open("%s.compare.out" % (test), 'w')
                    cf.write("WARNING: run did not produce any output\n")
                    cf.write("         unable to do a comparison\n")
                    cf.close()

        else:

            print "  storing output of %s as the new benchmark..." % (test)
            print "     new benchmark file: ", compareFile

            if (not compareFile == ""):
                os.system("cp -rf %s %s" % (compareFile, benchDir))

                cf = open("%s.status" % (test), 'w')
                cf.write("benchmarks updated.  New file:  %s\n" % (compareFile) )
                cf.close()
            else:
                cf = open("%s.status" % (test), 'w')
                cf.write("benchmarks failed")
                cf.close()

                
        #----------------------------------------------------------------------
        # do any requested visualization (2- and 3-d only)
        #---------------------------------------------------------------------- 
        doVis = getParam(test + ".doVis")

        if (doVis and not make_benchmarks):

           visVar = getParam(test + ".visVar")

           if (dim == 2):
              os.system('%s/%s --palette %s/Palette -cname "%s" -p "%s"' %
                        (compareToolDir, vis2dExecutable, compareToolDir, 
			 visVar, compareFile) )
           elif (dim == 3):
              os.system('%s/%s --palette %s/Palette -n 1 -cname "%s" -p "%s"' %
                        (compareToolDir, vis3dExecutable, compareToolDir, 
			 visVar, compareFile) )
           else:
              print "Visualization not supported for dim = %d" % (dim)


           # convert the .ppm files into .png files
           ppmFile = getRecentFileName(outputDir, "", ".ppm")

           os.system("convert %s `basename %s .ppm`.png" % (ppmFile, ppmFile) )
        
                       
        #----------------------------------------------------------------------
        # move the output files into the web directory
        #----------------------------------------------------------------------
        if (not make_benchmarks):
            shutil.copy("%s.run.out"     % (test), fullWebDir)
            shutil.copy("%s.make.out"    % (test), fullWebDir)
            shutil.copy("%s.compare.out" % (test), fullWebDir)

            shutil.copy(inputsFile, "%s/%s.%s" % (fullWebDir, test, inputsFile) )

            if (sourceTree == "Parallel"):
                shutil.copy(probinFile, "%s/%s.%s" % (fullWebDir, test, probinFile) )

            for file in auxFiles:
                shutil.copy(file, "%s/%s.%s" % (fullWebDir, test, file) )

            if (doVis):
               pngFile = getRecentFileName(outputDir, "", ".png")
               shutil.copy(pngFile, fullWebDir)
               
        else:
            shutil.copy("%s.status" % (test), fullWebDir)
            

        #----------------------------------------------------------------------
        # archive the output
        #----------------------------------------------------------------------
        for file in os.listdir(outputDir):
            if (os.path.isdir(file) and 
                (file.startswith("%s_plt" % (test)) or 
                 file.startswith("%s_chk" % (test)) ) ):
                err = os.system("tar czf %s.tgz %s >& %s.tar.out" % (file, file, file))
                if (err == 0):
                    shutil.rmtree(file)
                    os.remove("%s.tar.out" % (file))
                    

        #----------------------------------------------------------------------
        # write the report for this test
        #----------------------------------------------------------------------
        if (not make_benchmarks):
            print "  creating problem test report ..."
            reportSingleTest(sourceTree, test, testDir, fullWebDir)


        print "\n"
    
        
    #--------------------------------------------------------------------------
    # write the report for this instance of the test suite
    #--------------------------------------------------------------------------
    print "creating new test report..."
    reportThisTestRun(suiteName, make_benchmarks, comment, note, cvsTime, tests, testDir, testFile, fullWebDir)


    # make sure that all of the files in the web directory are world readable
    for file in os.listdir(fullWebDir):    
       currentFile = fullWebDir + file

       if (os.path.isfile(currentFile)):
          os.chmod(currentFile, 0644)
          

    #--------------------------------------------------------------------------
    # generate the master report for all test instances 
    #--------------------------------------------------------------------------
    reportAllRuns(suiteName, webTopDir)



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
def reportSingleTest(sourceTree, testName, testDir, fullWebDir):
    """ generate a single problem's test result page """
    
    # get the current directory
    currentDir = os.getcwd()
    
    # switch to the web directory and open the report file
    os.chdir(fullWebDir)

    #--------------------------------------------------------------------------
    # parse the compilation report and determine if we compiled
    #--------------------------------------------------------------------------
    compileFile = "%s.make.out" % (testName)

    try:
        cf = open(compileFile, 'r')

    except IOError:
        print "WARNING: no compilation file found"
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
    compileTest = getParam(testName + ".compileTest")

    if (not compileTest):
        compareFile = "%s.compare.out" % (testName)
    
        try:
            cf = open(compareFile, 'r')

        except IOError:
            print "WARNING: no comparison file found"
            compareSuccessful = 0
            diffLines = ['']

        else:
            diffLines = cf.readlines()

        # successful comparison is indicated by PLOTFILES AGREE
            compareSuccessful = 0
    
            for line in diffLines:
                if (string.find(line, "PLOTFILES AGREE") >= 0):
                    compareSuccessful = 1
                    break
    
            cf.close()


    #--------------------------------------------------------------------------
    # write out the status file for this problem, with either
    # PASSED or FAILED
    #--------------------------------------------------------------------------
    statusFile = "%s.status" % (testName)
    sf = open(statusFile, 'w')

    if (compileSuccessful and 
        (compileTest or
         (not compileTest and compareSuccessful))):
        sf.write("PASSED\n")
        print "    %s PASSED" % (testName)
    else:
        sf.write("FAILED\n")
        print "    %s FAILED" % (testName)

    sf.close()


    #--------------------------------------------------------------------------
    # generate the HTML page for this test
    #--------------------------------------------------------------------------

    # check to see if the CSS file is present, if not, write it
    if (not os.path.isfile("tests.css")):
        create_css()

    
    htmlFile = "%s.html" % (testName)
    hf = open(htmlFile, 'w')

    newHead = HTMLHeader + r"""<CENTER><H1><A HREF="index.html">@TESTDIR@</A> / @TESTNAME@</H1></CENTER>"""

    newHead = newHead.replace("@TESTDIR@", testDir)
    newHead = newHead.replace("@TESTNAME@", testName)

    hf.write(newHead)


    if (not compileTest):

        useMPI = getParam("%s.useMPI" % (testName))
        if (useMPI):

            numprocs = getParam("%s.numprocs" % (testName))

            hf.write("<P><b>Parallel Run</b><br>numprocs = %d\n" % (numprocs) )
            hf.write("<P>&nbsp;\n")
       

        # is this a restart test?
        restart = getParam("%s.restartTest" % (testName) )
        if (restart):

            restartFileNum = getParam("%s.restartTest" % (testName) )
       
            hf.write("<P><b>Restart Test</b><br>Job was run as normal and then restarted from checkpoint # %d, and the two final outputs were compared\n" % (restartFileNum) )

        hf.write("<P>&nbsp;\n")       

        # write out the information about the test
        buildDir = getParam("%s.buildDir" % (testName) )
        hf.write("<P><b>build directory:</b> %s\n" % (buildDir) )

        inputFile = getParam("%s.inputFile" % (testName) )
        hf.write("<P><b>input file:</b> <A HREF=\"%s.%s\">%s</A>\n" %
                 (testName, inputFile, inputFile) )

        if (sourceTree == "Parallel"):
            probinFile = getParam("%s.probinFile" % (testName) )
            hf.write("<P><b>probin file:</b> <A HREF=\"%s.%s\">%s</A>\n" %
                     (testName, probinFile, probinFile) )    


        auxFiles = getAuxFiles(testName)

        i = 1
        for file in auxFiles:
            hf.write("<P><b>aux%dFile:</b> <A HREF=\"%s.%s\">%s</A>\n" %
                     (i, testName, file, file) )
            i = i + 1
        

        dim = getParam("%s.dim" % (testName) )
        hf.write("<P><b>dimensionality:</b> %s\n" % (dim) )

        hf.write("<P>&nbsp;\n")

    
    # write out the compilation report
    if (compileSuccessful):
        hf.write("<P><H3 CLASS=\"passed\">Compilation Successful</H3></P>\n")
    else:
        hf.write("<P><H3 CLASS=\"failed\">Compilation Failed</H3></P>\n")

    hf.write("<A HREF=\"%s.make.out\">make output</A>\n" % (testName) )

    hf.write("<P>&nbsp;\n")

    
    if (not compileTest):

        # write out the comparison report
        if (compareSuccessful):
            hf.write("<P><H3 CLASS=\"passed\">Comparison Successful</H3></P>\n")
        else:
            hf.write("<P><H3 CLASS=\"failed\">Comparison Failed</H3></P>\n")

        hf.write("<A HREF=\"%s.run.out\">execution output</A>\n" % (testName) )


        hf.write("<P>&nbsp;\n")
        hf.write("<PRE>\n")
    
        for line in diffLines:
            hf.write(line)

        hf.write("</PRE>\n")

        # show any visualizations
        doVis = getParam("%s.doVis" % (testName) )
        if (doVis):
            pngFile = getRecentFileName(fullWebDir, testName, ".png")
            hf.write("<P>&nbsp;\n")
            hf.write("<P><IMG SRC='%s' BORDER=0>" % (pngFile) )
    

    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")    

    hf.close()
    
    
    # switch back to the original directory
    os.chdir(currentDir)
	

#==============================================================================
# reportThisTestRun
#==============================================================================
def reportThisTestRun(suiteName, make_benchmarks, comment, note, cvsTime, tests, testDir, testFile, fullWebDir):
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

    newHead = newHead.replace("@TESTDIR@", suiteName)
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
    hf.write("<p>&nbsp;&nbsp;<b>cvs update on Parallel/:</b> <A HREF=\"%s\">%s</A>\n" %
             ("cvs.Parallel.out", "cvs.Parallel.out") )
    hf.write("<p>&nbsp;&nbsp;<b>cvs update on fParallel/:</b> <A HREF=\"%s\">%s</A>\n" %
             ("cvs.fParallel.out", "cvs.fParallel.out") )        
    hf.write("<p>&nbsp;\n")

    sourceTree = getParam("main.sourceTree")
    if (sourceTree == "Parallel"):
            hf.write("<p>&nbsp;&nbsp;<b>Parallel/ ChangeLog:</b> <A HREF=\"%s\">%s</A>\n" %
                     ("ChangeLog.Parallel", "ChangeLog.Parallel") )
            
    hf.write("<p>&nbsp;&nbsp;<b>fParallel/ ChangeLog:</b> <A HREF=\"%s\">%s</A>\n" %
             ("ChangeLog.fParallel", "ChangeLog.fParallel") )        
    hf.write("<p>&nbsp;\n")    

    hf.write("<P><TABLE BORDER=0 CELLPADDING=3>\n")
    
    # loop over the tests and add a line for each
    for test in tests:

        if (not make_benchmarks):
            
            # check if it passed or failed
            statusFile = "%s.status" % (test)

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
                     (test, test) )
        
            if (testPassed):
                hf.write("<TD><H3 class=\"passed\">PASSED</H3></TD></TR>\n")
            else:
                hf.write("<TD><H3 class=\"failed\">FAILED</H3></TD></TR>\n")

        
            hf.write("<TR><TD>&nbsp;</TD></TR>\n")


        else:

            if (getParam(test + ".restartTest")):
                continue

            if (getParam(test + ".compileTest")):
                continue

                
            # the benchmark was updated -- find the name of the new benchmark file
            benchStatusFile = "%s.status" % (test)

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
                          (test, benchFile) )
            else:
                 hf.write("<TR><TD>%s</TD><TD>&nbsp;</TD><TD><H3 class=\"failed\">BENCHMARK NOT UPDATED</H3></TD><TD>&nbsp;</TD><TD>(compilation or execution failed)</TD></TR>\n" %
                          (test) )
                 
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
def reportAllRuns(suiteName, webTopDir):

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

    title = "%s regression tests" % (suiteName)
    
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

    test(sys.argv)
