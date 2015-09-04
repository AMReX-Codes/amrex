#!/usr/bin/env python

"""
A simple regression test framework for a BoxLib-based code

There are several major sections to this source: the runtime parameter
routines, the test suite routines, and the report generation routines.
They are separated as such in this file.

This test framework understands source based out of the F_Src and
C_Src BoxLib frameworks.

"""

import ConfigParser
import datetime
import email
import argparse
import getpass
import os
import shlex
import shutil
import smtplib
import socket
import string
import subprocess
import sys
import tarfile
import time


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# T E S T   C L A S S E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class Test:

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

        self.compare_file_used = ""

        self.diffDir = ""
        self.diffOpts = ""

        self.addToCompileString = ""

        self.reClean = 0    # set automatically, not by users

        self.wallTime = 0   # set automatically, not by users

        self.nlevels = None  # set but running fboxinfo on the output

    def __cmp__(self, other):
        return cmp(self.value(), other.value())

    def value(self):
        return self.name


class Suite:

    def __init__ (self, args):

        self.args = args

        self.test_file_path = os.getcwd() + '/' + self.args.input_file[0]

        self.suiteName = "testDefault"
        self.sub_title = ""

        self.sourceTree = ""
        self.boxLibDir = ""
        self.sourceDir = ""
        self.testTopDir = ""
        self.webTopDir = ""

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

        # delete all plot/checkfiles but the plotfile used for comparison upon
        # completion
        self.purge_output = 0

    def get_bench_dir(self):
        bench_dir = self.testTopDir + self.suiteName + "-benchmarks/"
        if not os.path.isdir(bench_dir):
            if not args.make_benchmarks == None:
                os.mkdir(bench_dir)
            else:
                fail("ERROR: benchmark directory, %s, does not exist" % (bench_dir))
        return bench_dir

    def make_test_dirs(self):
        os.chdir(self.testTopDir)

        today_date = datetime.date.today()
        today = today_date.__str__()

        # figure out what the current output directory should be
        maxRuns = 100      # maximum number of tests in a given day

        test_dir = today + "/"

        # test output stored in a directory suiteName-tests/2007-XX-XX/
        # make sure that the suiteName-tests directory exists
        if not os.path.isdir(self.testTopDir + self.suiteName + "-tests/"):
            os.mkdir(self.testTopDir + self.suiteName + "-tests/")

        full_test_dir = self.testTopDir + self.suiteName + "-tests/" + test_dir

        if self.args.do_temp_run:
            test_dir = "TEMP_RUN/"
            full_test_dir = self.testTopDir + self.suiteName + "-tests/" + test_dir
            shutil.rmtree(full_test_dir)
        else:
            for i in range(1, maxRuns):
                if not os.path.isdir(full_test_dir): break
                test_dir = today + "-{:03d}/".format(i)
                full_test_dir = self.testTopDir + self.suiteName + "-tests/" + test_dir

        bold("testing directory is: " + test_dir, skip_before=1)
        os.mkdir(full_test_dir)

        # make the web directory -- this is where all the output and HTML will be
        # put, so it is easy to move the entire test website to a different disk
        full_web_dir = "%s/%s/"  % (self.webTopDir, test_dir)

        if self.args.do_temp_run:
            shutil.rmtree(full_web_dir)

        os.mkdir(full_web_dir)

        # copy the test file into the web output directory
        shutil.copy(self.test_file_path, full_web_dir)

        return test_dir, full_test_dir, full_web_dir

    def get_last_run(self):
        """ return the name of the directory corresponding to the previous
            run of the test suite """

        outdir = self.testTopDir + self.suiteName + "-tests/"

        # this will work through 2099
        dirs = [d for d in os.listdir(outdir) if (os.path.isdir(outdir + d) and
                                                  d.startswith("20"))]
        dirs.sort()

        return dirs[-1]

    def get_test_failures(self, test_dir):
        """ look at the test run in testDir and return the list of tests that
            failed """

        cwd = os.getcwd()

        outdir = self.testTopDir + self.suiteName + "-tests/"

        os.chdir(outdir + test_dir)

        failed = []

        for test in os.listdir("."):
            if not os.path.isdir(test): continue

            # the status files are in the web dir
            status_file = self.webTopDir + test_dir + "/%s.status" % (test)
            sf = open(status_file, "r")
            for line in sf:
                if line.find("FAILED") >= 0:
                    failed.append(test)

        os.chdir(cwd)
        return failed

    def make_realclean(self):
        run("{} BOXLIB_HOME={} {} {} realclean".format(
            self.MAKE, self.boxLibDir,
            self.extSrcCompString, self.extraBuildDirCompString))

    def build_f(self, opts="", target="", outfile=None):
        comp_string = "{} -j{} BOXLIB_HOME={} COMP={} {} {}".format(
            self.MAKE, self.numMakeJobs, self.boxLibDir, self.FCOMP, opts, target)
        print "  " + comp_string
        run(comp_string, outfile=outfile)
        return comp_string

    def run_test(self, test, base_command):
        if test.useMPI:
            testRunCommand = ""
            if test.useOMP:
	        testRunCommand = "OMP_NUM_THREADS={} ".format(test.numthreads)
            testRunCommand += self.MPIcommand
	    testRunCommand = testRunCommand.replace("@host@", self.MPIhost)
	    testRunCommand = testRunCommand.replace("@nprocs@", "{}".format(test.numprocs))
            testRunCommand = testRunCommand.replace("@command@", base_command)

        elif test.useOMP:
            testRunCommand = "OMP_NUM_THREADS={} ".format(test.numthreads)
            testRunCommand += base_command

        else:
            testRunCommand = base_command

        print "    " + testRunCommand
        systemCall(testRunCommand)
        return testRunCommand


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# R U N T I M E   P A R A M E T E R   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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

    if isInt(string): return int(string)
    elif isFloat(string): return float(string)
    else: return string.strip()


#==============================================================================
def load_params(args):
    """
    reads the parameter file and creates as list of test objects as well as
    the suite object
    """

    testList = []

    cp = ConfigParser.ConfigParser()
    cp.optionxform = str

    bold("loading " + args.input_file[0])

    try: cp.read(args.input_file[0])
    except:
        fail("ERROR: unable to read parameter file {}".format(file))

    # "main" is a special section containing the global suite parameters.
    mysuite = Suite(args)

    valid_options = mysuite.__dict__.keys()
    valid_options += ["extraBuildDir", "extraBuildDir2"]

    for opt in cp.options("main"):

        # get the value of the current option
        value = convertType(cp.get("main", opt))

        if opt in valid_options:

            if opt == "sourceTree":
                if not value in ["C_Src", "F_Src", "BoxLib"]:
                    fail("ERROR: invalid sourceTree")
                else:
                    mysuite.sourceTree = value

            elif opt == "sourceDir": mysuite.sourceDir = checkTestDir(value)
            elif opt == "boxLibDir": mysuite.boxLibDir = checkTestDir(value)
            elif opt == "testTopDir": mysuite.testTopDir = checkTestDir(value)
            elif opt == "webTopDir": mysuite.webTopDir = os.path.normpath(value) + "/"
            elif opt == "extSrcDir":
                mysuite.extSrcDir = checkTestDir(value)
                mysuite.useExtSrc = 1

            elif opt == "extraBuildDir":
                # we will keep the extra build directories in a list -- we
                # want them in the correct order in the list so we can simply
                # index later
                if len(mysuite.extraBuildDirs) == 0:
                    mysuite.extraBuildDirs.append(checkTestDir(value))
                else:
                    mysuite.extraBuildDirs.insert(0,checkTestDir(value))

            elif opt == "extraBuildDir2":
                mysuite.extraBuildDirs.append(checkTestDir(value))

            elif opt == "emailTo": mysuite.emailTo = value.split(",")

            else:
                # generic setting of the object attribute
                setattr(mysuite, opt, value)

        else:
            warning("WARNING: suite parameter %s not valid" % (opt))


    mysuite.useExtraBuild = len(mysuite.extraBuildDirs)
    mysuite.srcName = os.path.basename(os.path.normpath(mysuite.sourceDir))

    # update the additional compilation str for additional source dir
    if mysuite.useExtSrc:
        mysuite.extSrcName = os.path.basename(os.path.normpath(mysuite.extSrcDir))

        if mysuite.extSrcCompString != "":
            mysuite.extSrcCompString += "="+mysuite.extSrcDir

    # update additional compiled string for any extra build directory
    if mysuite.useExtraBuild > 0:
        for n in range(len(mysuite.extraBuildDirs)):
            mysuite.extraBuildNames.append(os.path.basename(os.path.normpath(mysuite.extraBuildDirs[n])))

        # since we are building in the extraBuildDir, we need to
        # tell make where the sourceDir is
        if mysuite.extraBuildDirCompString != "":
            mysuite.extraBuildDirCompString += "="+mysuite.sourceDir


    # BoxLib-only tests don't have a sourceDir
    if mysuite.sourceTree == "BoxLib": mysuite.sourceDir = mysuite.boxLibDir

    # checks
    if mysuite.sendEmailWhenFail:
        if mysuite.emailTo == [] or mysuite.emailBody == "":
            fail("ERROR: when sendEmailWhenFail = 1, you must specify emailTo and emailBody\n")

        if mysuite.emailFrom == "":
            mysuite.emailFrom = '@'.join((getpass.getuser(), socket.getfqdn()))

        if mysuite.emailSubject == "":
            mysuite.emailSubject = mysuite.suiteName+" Regression Test Failed"

    if (mysuite.sourceTree == "" or mysuite.boxLibDir == "" or
        mysuite.sourceDir == "" or mysuite.testTopDir == ""):
        fail("ERROR: required suite-wide directory not specified\n" + \
                 "(sourceTree, boxLibDir, sourceDir, testTopDir)")

    # Make sure the web dir is valid (or use the default is none specified)
    if mysuite.webTopDir == "":
        mysuite.webTopDir = "%s/%s-web/" % (mysuite.testTopDir, mysuite.suiteName)

    if not os.path.isdir(mysuite.webTopDir):
        try: os.mkdir(mysuite.webTopDir)
        except: fail("ERROR: unable to create the web directory: %s\n" %
                     (mysuite.webTopDir))

    # all other sections are tests
    bold("finding tests and checking parameters...", skip_before=1)

    for sec in cp.sections():

        if sec == "main": continue

        print "  %s" % (sec)

        # maximum test name length -- used for HTML formatting
        mysuite.lenTestName = max(mysuite.lenTestName, len(sec))

        # create the test object for this test
        mytest = Test(sec)

        invalid = 0

        # set the test object data by looking at all the options in
        # the current section of the parameter file
        valid_options = mytest.__dict__.keys()
        valid_options += ["aux1File", "aux2File", "aux3File"]
        valid_options += ["link1File", "link2File", "link3File"]

        for opt in cp.options(sec):

            # get the value of the current option
            value = convertType(cp.get(sec, opt))

            if opt in valid_options:

                if opt in ["aux1File", "aux2File", "aux3File"]:
                    mytest.auxFiles.append(value)

                elif opt in ["link1File", "link2File", "link3File"]:
                    mytest.linkFiles.append(value)

                else:
                    # generic setting of the object attribute
                    setattr(mytest, opt, value)

            else:
                warning("   WARNING: unrecognized parameter %s for test %s" % (opt, sec))


        # make sure that the build directory actually exists
        if mytest.useExtraBuildDir > 0:
            bDir = mysuite.extraBuildDirs[mytest.useExtraBuildDir-1] + mytest.buildDir
        else:
            bDir = mysuite.sourceDir + mytest.buildDir

        if not os.path.isdir(bDir):
            warning("   WARNING: invalid build directory: %s" % (bDir))
            invalid = 1


        # make sure all the require parameters are present
        if mytest.compileTest:
            if mytest.buildDir == "":
                warning("   WARNING: mandatory parameters for test %s not set" % (sec))
                invalid = 1

        else:
            if (mytest.buildDir == "" or mytest.inputFile == "" or
                (mysuite.sourceTree == "C_Src" and mytest.probinFile == "") or
                mytest.dim == -1):
                warning("   WARNING: mandatory parameters for test %s not set" % (sec))
                warning("            buildDir = %s" % (mytest.buildDir))
                warning("            inputFile = %s" % (mytest.inputFile))
                if (mysuite.sourceTree == "C_Src"):
                    warning("            probinFile = %s" % (mytest.probinFile))
                warning("            dim = %s" % (mytest.dim))

                invalid = 1

        # check the optional parameters
        if mytest.restartTest and mytest.restartFileNum == -1:
            warning("   WARNING: restart-test %s needs a restartFileNum" % (sec))
            invalid = 1

        if mytest.selfTest and mytest.stSuccessString == "":
            warning("   WARNING: self-test %s needs a stSuccessString" % (sec))
            invalid = 1

        if mytest.useMPI and mytest.numprocs == -1:
            warning("   WARNING: MPI parallel test %s needs numprocs" % (sec))
            invalid = 1

        if mytest.useOMP and mytest.numthreads == -1:
            warning("   WARNING: OpenMP parallel test %s needs numthreads" % (sec))
            invalid = 1

        if mytest.doVis and mytest.visVar == "":
            warning("   WARNING: test %s has visualization, needs visVar" % (sec))
            invalid = 1

        if mysuite.sourceTree == "BoxLib" and mytest.testSrcTree == "":
            warning("   WARNING: test %s is a BoxLib test but testSrcTree not set" % (sec))
            invalid = 1


        # add the current test object to the master list
        if not invalid:
            testList.append(mytest)
        else:
            warning("   WARNING: test %s will be skipped" % (sec))


    # if any runs are parallel, make sure that the MPIcommand is defined
    anyMPI = any([t.useMPI for t in testList])

    if anyMPI and mysuite.MPIcommand == "":
        fail("ERROR: some tests are MPI parallel, but MPIcommand not defined")

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
    print termColors.FAIL + str + termColors.ENDC
    sys.exit()

def testfail(str):
    print termColors.FAIL + str + termColors.ENDC

def warning(str):
    print termColors.WARNING + str + termColors.ENDC

def success(str):
    print termColors.SUCCESS + str + termColors.ENDC

def bold(str, skip_before=0):
    if skip_before == 1: print ""
    print termColors.BOLD + str + termColors.ENDC


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# S Y S T E M   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def systemCall(string):
    status = os.system('bash -c "' + string + '"')
    return status


def run(string, stdin=False, outfile=None):

    # shlex.split will preserve inner quotes
    prog = shlex.split(string)
    if stdin:
        p0 = subprocess.Popen(prog, stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
    else:
        p0 = subprocess.Popen(prog, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)

    stdout0, stderr0 = p0.communicate()
    if stdin: p0.stdin.close()
    rc = p0.returncode
    p0.stdout.close()

    if not outfile == None:
        try: cf = open(outfile, "w")
        except IOError:
            fail("  ERROR: unable to open file for writing")
        else:
            for line in stdout0:
                cf.write(line)
            cf.close()

    return stdout0, stderr0, rc


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# T E S T   S U I T E   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

def findBuildDirs(tests):
    """ given the list of test objects, find the set of UNIQUE build
        directories.  Note if we have the useExtraBuildDir flag set """

    build_dirs = []
    reClean = []

    for obj in tests:

        # keep track of the build directory and which source tree it is
        # in (e.g. the extra build dir)

        # first find the list of unique build directories
        dir_pair = (obj.buildDir, obj.useExtraBuildDir)
        if build_dirs.count(dir_pair) == 0:
            build_dirs.append(dir_pair)


        # re-make all problems that specify an extra compile argument,
        # just to make sure that any unique build commands are seen.
        if not obj.addToCompileString == "":
            reClean.append(dir_pair)

    for bdir, _ in reClean:
        for obj in tests:
            if obj.buildDir == bdir:
                obj.reClean = 1

    return build_dirs


#==============================================================================
def getLastPlotfile(outputDir, test):
    """given an output directory and the test name, find the last
       plotfile written.  Note: we give an error if the last
       plotfile is 0 """

    plts = [d for d in os.listdir(outputDir) if (os.path.isdir(d) and
                                                 d.startswith("{}_plt".format(test.name)))]
    if len(plts) == 0:
        warning("WARNING: test did not produce any output")
        return ""

    plts.sort()
    last_plot = plts.pop()

    if last_plot.endswith("00000"):
        warning("WARNING: only plotfile 0 was output -- skipping comparison")
        return ""

    return last_plot


#==============================================================================
def getRecentFileName(dir, base, extension):
    """ given the base and extension, find the most recent corresponding
    file """

    files = [f for f in os.listdir(dir) if (f.startswith(base) and
                                            f.endswith(extension))]

    files.sort(key=lambda x: os.path.getmtime(x))

    try: return files.pop()
    except: return None


#==============================================================================
def checkTestDir(dir_name):
   """ given a string representing a directory, check if it points to
       a valid directory.  If so, return the directory name """

   dir_name = os.path.normpath(dir_name) + "/"

   if not os.path.isdir(dir_name):
       fail("ERROR: {} is not a valid directory".format(dir_name))

   return dir_name


#==============================================================================
def doGITUpdate(topDir, root, outDir, gitbranch, githash):
   """ do a git update of the repository in topDir.  root is the name
       of the directory (used for labeling).  outDir is the full path
       to the directory where we will store the git output.  If githash
       is not empty, then we will check out that version instead of
       git-pulling."""

   os.chdir(topDir)

   # find out current branch so that we can go back later if we need.
   stdout0, stderr0, rc = run("git rev-parse --abbrev-ref HEAD")
   currentBranch = stdout0.rstrip('\n')

   if currentBranch != gitbranch:
       bold("git checkout %s in %s" % (gitbranch, topDir), skip_before=1)
       stdout, stderr, rc = run("git checkout {}".format(gitbranch), stdin=True)


   if githash == "" or githash == None:
       bold("'git pull' in %s" % (topDir), skip_before=1)

       # we need to be tricky here to make sure that the stdin is
       # presented to the user to get the password.
       stdout, stderr, rc = run("git pull", stdin=True,
                                outfile="git.{}.out".format(root))

   else:
       stdout, stderr, rc = run("git checkout {}".format(githash),
                                outfile="git.{}.out".format(root))

   # not sure if this is valid -- we are piping stderr into stdout
   # -- we should check the return code instead
   if stdout == "":
       fail("  ERROR: git update was unsuccessful")

   shutil.copy("git.{}.out".format(root),  outDir)

   return currentBranch


#==============================================================================
def saveGITHEAD(topDir, root, outDir):

   os.chdir(topDir)

   bold("saving git HEAD for %s/" % (root), skip_before=1)

   run("git rev-parse HEAD", outfile="git.{}.HEAD".format(root) )
   shutil.copy("git.{}.HEAD".format(root),  outDir)


#==============================================================================
def doGITback(topDir, root, gitbranch):
   """ do a git checkout of gitbranch in topDir.  root is the name
       of the directory (used for labeling). """

   os.chdir(topDir)

   bold("git checkout %s in %s" % (gitbranch, topDir), skip_before=1)

   stdout, stderr, rc = run("git checkout {}".format(gitbranch), stdin=True,
                            outfile="git.{}.out".format(root))

   if stdout == "":
       fail("  ERROR: git checkout was unsuccessful")


#==============================================================================
def makeGITChangeLog(gitDir, root, outDir):
    """ generate a ChangeLog git repository named root.  outDir is the
        full path to the directory where we will store the git output"""

    os.chdir(gitDir)

    bold("generating ChangeLog for %s/" % (root), skip_before=1)

    run("git log --name-only", outfile="ChangeLog.{}".format(root) )
    shutil.copy("ChangeLog.{}".format(root), outDir)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# test
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def testSuite(argv):

    usage = """
    testnew.py -h for options

    input file structure

      testfile.ini
          This is the input file that defines the tests that the
          suite knows about.  It has the format

            [main]
            boxLibDir      = < directory to the BoxLib/ directory >
            sourceDir      = < directory to the main Source directory >
            testTopDir     = < full path to test output directory >
            webTopDir      = < full path to test web directory >
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
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", type=int, default=-1,
                        help="restrict tests to a particular dimensionality")
    parser.add_argument("--make_benchmarks", type=str, default=None, metavar="comment",
                        help="make new benchmarks? (must provide a comment)")
    parser.add_argument("--no_update", type=str, default="None", metavar="name",
                        help="which codes to exclude from the git update? (None, All, or a comma-separated list of codes)")
    parser.add_argument("--single_test", type=str, default="", metavar="test-name",
                        help="name of a single test to run")
    parser.add_argument("--tests", type=str, default="", metavar="'test1 test2 test3'",
                        help="a space-separated list of tests to run")
    parser.add_argument("--do_temp_run", action="store_true",
                        help="is this a temporary run? (output not stored or logged)")
    parser.add_argument("--boxLibGitHash", type=str, default=None, metavar="hash",
                        help="git hash of a version of BoxLib.  If provided, this version will be used to run tests.")
    parser.add_argument("--sourceGitHash", type=str, default=None, metavar="hash",
                        help="git hash of a version of the source code.  For BoxLib tests, this will be ignored.")
    parser.add_argument("--extSrcGitHash", type=str, default=None, metavar="hash",
                        help="git hash of a version of the source code.  For BoxLib tests, this will be ignored.")
    parser.add_argument("--note", type=str, default="",
                        help="a note on the resulting test webpages")
    parser.add_argument("--complete_report_from_crash", type=str, default="", metavar="testdir",
                        help="complete the generation of the report from a crashed test suite run named testdir")
    parser.add_argument("--redo_failed", action="store_true",
                        help="only run the tests that failed last time")
    parser.add_argument("input_file", metavar="input-file", type=str, nargs=1,
                        help="the input file (INI format) containing the suite and test parameters")

    args=parser.parse_args()

    #--------------------------------------------------------------------------
    # read in the test information
    #--------------------------------------------------------------------------
    suite, testList = load_params(args)

    defined_tests = testList[:]

    # if we only want to run the tests that failed previously, remove the
    # others
    if args.redo_failed:
        last_run = suite.get_last_run()
        failed = suite.get_test_failures(last_run)

        testListold = testList[:]
        testList = [t for t in testListold if t.name in failed]

    # if we only want to run tests of a certain dimensionality, remove
    # the others
    if args.d in [1, 2, 3]:
        testListold = testList[:]
        testList = [t for t in testListold if t.dim == args.d]

    if len(testList) == 0:
        fail("No valid tests defined")

    activeTestList = [t.name for t in defined_tests]

    if not args.complete_report_from_crash == "":

        # make sure the web directory from the crash run exists
        full_web_dir = "%s/%s/"  % (suite.webTopDir, args.complete_report_from_crash)
        if not os.path.isdir(full_web_dir):
            fail("Crash directory does not exist")

        # find all the tests that completed in that web directory
        tests = []
        testFile = ""
        wasBenchmarkRun = 0
        for file in os.listdir(full_web_dir):
            if os.path.isfile(file) and file.endswith(".status"):
                index = string.rfind(file, ".status")
                tests.append(file[:index])

                f = open(full_web_dir + file, "r")
                for line in f:
                    if line.find("benchmarks updated") > 0:
                        wasBenchmarkRun = 1

            if os.path.isfile(file) and file.endswith(".ini"):
                testFile = file


        # create the report for this test run
        numFailed = reportThisTestRun(suite, wasBenchmarkRun,
                                      "recreated report after crash of suite",
                                      "",  0, 0, 0,
                                      tests, args.complete_report_from_crash, testFile, full_web_dir)


        # create the suite report
        bold("creating suite report...")
        tableHeight = min(max(suite.lenTestName, 4), 16)
        reportAllRuns(suite, activeTestList, suite.webTopDir, tableHeight=tableHeight)
        sys.exit("done")


    #--------------------------------------------------------------------------
    # figure out which git repos we will update
    #--------------------------------------------------------------------------
    no_update_low = args.no_update.lower()

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

        updateBoxLib = True
        updateSource = True
        updateExtSrc = True

        if "boxlib" in nouplist: updateBoxLib = False
        if suite.srcName.lower() in nouplist: updateSource = False
        if suite.extSrcName.lower() in nouplist: updateExtSrc = False

        # each extra build directory has its own update flag
        updateExtraBuild = []
        for e in suite.extraBuildNames:
            if e.lower() in nouplist:
                updateExtraBuild.append(False)
            else:
                updateExtraBuild.append(True)

    if not suite.useExtSrc: updateExtSrc = False
    if not suite.useExtraBuild: updateExtraBuild = [False]

    if suite.sourceTree == "BoxLib":
        updateSource = False # to avoid updating BoxLib twice.
                             # The update of BoxLib is controlled by updateBoxLib
        args.sourceGitHash = ""

    if args.boxLibGitHash: updateBoxLib = False
    if args.sourceGitHash: updateSource = False
    if args.extSrcGitHash: updateExtSrc = False

    #--------------------------------------------------------------------------
    # if we are doing a single test, remove all other tests
    # if we specified a list of tests, check each one
    # if we did both --single_test and --tests, complain
    #--------------------------------------------------------------------------
    if not args.single_test == "" and not args.tests == "":
        fail("ERROR: specify tests either by --single_test or --tests, not both")

    if not args.single_test == "":
        testsFind = [args.single_test]
    elif not args.tests == "":
        testsFind = args.tests.split()
    else:
        testsFind = []

    if len(testsFind) > 0:
        new_test_list = []
        for test in testsFind:
            _tmp = [o for o in testList if o.name == test]
            if len(_tmp) == 1:
                new_test_list += _tmp
            else:
                fail("ERROR: {} is not a valid test".format(test))

        testList = new_test_list


    #--------------------------------------------------------------------------
    # check bench dir and create output directories
    #--------------------------------------------------------------------------
    all_compile = all([t.compileTest == 1 for t in testList])

    if not all_compile:
        bench_dir = suite.get_bench_dir()

    testDir, full_test_dir, full_web_dir = suite.make_test_dirs()

    #--------------------------------------------------------------------------
    # do the source updates
    #--------------------------------------------------------------------------
    now = time.localtime(time.time())
    updateTime = time.strftime("%Y-%m-%d %H:%M:%S %Z", now)

    os.chdir(suite.testTopDir)

    if updateSource or args.sourceGitHash:

        # main suite
        sourceGitBranch_Orig = doGITUpdate(suite.sourceDir,
                                           suite.srcName, full_web_dir,
                                           suite.sourceGitBranch,
                                           args.sourceGitHash)

    if updateExtSrc or args.extSrcGitHash:

        # extra source
        if suite.useExtSrc:
            extSrcGitBranch_Orig = doGITUpdate(suite.extSrcDir,
                                               suite.extSrcName, full_web_dir,
                                               suite.extSrcGitBranch,
                                               args.extSrcGitHash)

    if any(updateExtraBuild):

        # extra build directory
        for n in range(suite.useExtraBuild):
            if updateExtraBuild[n]:
                extSrcGitBranch_Orig = doGITUpdate(suite.extraBuildDirs[n],
                                                   suite.extraBuildNames[n], full_web_dir,
                                                   "master",
                                                   "")

    if updateBoxLib or args.boxLibGitHash:

        # BoxLib
        boxLibGitBranch_Orig = doGITUpdate(suite.boxLibDir,
                                           "BoxLib", full_web_dir,
                                           suite.boxLibGitBranch,
                                           args.boxLibGitHash)

    #--------------------------------------------------------------------------
    # Save git HEADs
    #--------------------------------------------------------------------------
    saveGITHEAD(suite.boxLibDir, "BoxLib", full_web_dir)

    if suite.sourceTree != "BoxLib":
        saveGITHEAD(suite.sourceDir, suite.srcName, full_web_dir)

    if suite.useExtSrc:
        saveGITHEAD(suite.extSrcDir, suite.extSrcName, full_web_dir)

    for n in range(suite.useExtraBuild):
        saveGITHEAD(suite.extraBuildDirs[n], suite.extraBuildNames[n], full_web_dir)


    #--------------------------------------------------------------------------
    # generate the ChangeLogs
    #--------------------------------------------------------------------------
    if updateSource:

        # main suite
        makeGITChangeLog(suite.sourceDir, suite.srcName, full_web_dir)

    if updateExtSrc:

        # extra source
        if (suite.useExtSrc):
            makeGITChangeLog(suite.extSrcDir, suite.extSrcName, full_web_dir)


    # extra build directories
    for n in range(suite.useExtraBuild):
        if updateExtraBuild[n]:
            makeGITChangeLog(suite.extraBuildDirs[n], suite.extraBuildNames[n], full_web_dir)

    if updateBoxLib:

        # BoxLib
        makeGITChangeLog(suite.boxLibDir, "BoxLib", full_web_dir)


    #--------------------------------------------------------------------------
    # build the comparison and visualization tools
    #--------------------------------------------------------------------------
    bold("building the comparison tools...", skip_before=1)

    suite.compareToolDir = os.path.normpath(suite.boxLibDir) + "/Tools/Postprocessing/F_Src/"

    os.chdir(suite.compareToolDir)

    suite.make_realclean()

    suite.build_f(target="programs=fcompare", opts="NDEBUG=t MPI= ")
    compareExecutable = getRecentFileName(suite.compareToolDir,"fcompare",".exe")
    shutil.copy(compareExecutable, full_test_dir + "/fcompare.exe")

    suite.build_f(target="programs=fboxinfo", opts="NDEBUG=t MPI= ")
    compareExecutable = getRecentFileName(suite.compareToolDir,"fboxinfo",".exe")
    shutil.copy(compareExecutable, full_test_dir + "/fboxinfo.exe")

    if any([t for t in testList if t.dim == 2]):
        bold("building the 2-d visualization tools...", skip_before=1)
        suite.build_f(target="programs=fsnapshot2d", opts="NDEBUG=t MPI= ")
        vis2dExecutable = getRecentFileName(suite.compareToolDir,"fsnapshot2d",".exe")

    if any([t for t in testList if t.dim == 3]):
        bold("building the 3-d visualization tools...", skip_before=1)
        suite.build_f(opts="NDEBUG=t MPI= ", target="programs=fsnapshot3d")
        vis3dExecutable = getRecentFileName(suite.compareToolDir,"fsnapshot3d",".exe")


    #--------------------------------------------------------------------------
    # output test list
    #--------------------------------------------------------------------------
    bold("running tests: ", skip_before=1)
    for obj in testList:
        print "  %s " % obj.name


    #--------------------------------------------------------------------------
    # do a make clean, only once per build directory
    #--------------------------------------------------------------------------
    all_build_dirs = findBuildDirs(testList)

    bold("make clean in...", skip_before=1)

    for dir, source_tree in all_build_dirs:

        if source_tree > 0:
            print "  {} in {}".format(dir, suite.extraBuildNames[source_tree-1])
            os.chdir(suite.extraBuildDirs[source_tree-1] + dir)
        else:
            print "  {}".format(dir)
            os.chdir(suite.sourceDir + dir)

        suite.make_realclean()

    os.chdir(suite.testTopDir)


    #--------------------------------------------------------------------------
    # main loop over tests
    #--------------------------------------------------------------------------
    for test in testList:

        bold("working on test: {}".format(test.name), skip_before=1)

        if not args.make_benchmarks == None and (test.restartTest or test.compileTest or
                                test.selfTest):
            warning("  WARNING: test {} doesn't need benchmarks".format(test.name))
            warning("           skipping\n")
            continue


        #----------------------------------------------------------------------
        # make the run directory
        #----------------------------------------------------------------------
        outputDir = full_test_dir + test.name + '/'
        os.mkdir(outputDir)


        #----------------------------------------------------------------------
        # compile the code
        #----------------------------------------------------------------------
        if test.useExtraBuildDir > 0:
            bDir = suite.extraBuildDirs[test.useExtraBuildDir-1] + test.buildDir
        else:
            bDir = suite.sourceDir + test.buildDir

        os.chdir(bDir)

        if test.reClean == 1:
            # for one reason or another, multiple tests use different
            # build options, make clean again to be safe
            print "  re-making clean..."
            suite.make_realclean()


        print "  building..."

        if suite.sourceTree == "C_Src" or test.testSrcTree == "C_Src":

            buildOptions = ""
            exeSuffix = ""

            if test.debug:
                buildOptions += "DEBUG=TRUE "
                exeSuffix += ".DEBUG"
            else:
                buildOptions += "DEBUG=FALSE "

            if test.useMPI:
                buildOptions += "USE_MPI=TRUE "
                exeSuffix += ".MPI"
            else:
                buildOptions += "USE_MPI=FALSE "

            if test.useOMP:
                buildOptions += "USE_OMP=TRUE "
                exeSuffix += ".OMP"
            else:
                buildOptions += "USE_OMP=FALSE "

            if test.useExtraBuildDir > 0:
                buildOptions += suite.extraBuildDirCompString + " "

            executable = "%s%dd" % (suite.suiteName, test.dim) + exeSuffix + ".ex"

            compString = "%s -j%s BOXLIB_HOME=%s %s %s DIM=%d %s COMP=%s FCOMP=%s executable=%s" % \
                (suite.MAKE, suite.numMakeJobs, suite.boxLibDir,
                 suite.extSrcCompString, test.addToCompileString,
                 test.dim, buildOptions, suite.COMP, suite.FCOMP,
                 executable)

            print "    " + compString
            so, se, r = run(compString,
                            outfile="{}/{}.make.out".format(outputDir, test.name))

        elif suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src":

            buildOptions = ""

            if test.debug:
                buildOptions += "NDEBUG= "
            else:
                buildOptions += "NDEBUG=t "

            if test.useMPI:
                buildOptions += "MPI=t "
            else:
                buildOptions += "MPI= "

            if test.useOMP:
                buildOptions += "OMP=t "
            else:
                buildOptions += "OMP= "

            if test.useExtraBuildDir > 0:
                buildOptions += suite.extraBuildDirCompString + " "

            compString = suite.build_f(opts="{} {} {}".format(
                suite.extSrcCompString, test.addToCompileString, buildOptions),
                          outfile="{}/{}.make.out".format(outputDir, test.name))

            # we need a better way to get the executable name here
            executable = getRecentFileName(bDir,"main",".exe")


        if test.compileTest:

            # compilation tests are done now -- just make the report and ...
            shutil.copy("%s/%s.make.out"    % (outputDir, test.name), full_web_dir)

            print "  creating problem test report ..."
            reportSingleTest(suite, test, compString, "", testDir, full_web_dir)

            # ... skip to the next test in the loop
            continue


        #----------------------------------------------------------------------
        # copy the necessary files over to the run directory
        #----------------------------------------------------------------------
        print "  copying files to run directory..."

        try: shutil.copy(executable, outputDir)
        except (IOError, AttributeError):

            # compilation failed.  First copy the make.out into the
            # web directory and then report
            shutil.copy("%s/%s.make.out" % (outputDir, test.name), full_web_dir)

            errorMsg = "    ERROR: compilation failed"
            reportTestFailure(errorMsg, test, testDir, full_web_dir, compString=compString)
            continue

        try: shutil.copy(test.inputFile, outputDir)
        except IOError:
            errorMsg = "    ERROR: unable to copy input file: %s" % test.inputFile
            reportTestFailure(errorMsg, test, testDir, full_web_dir)
            continue

	# sometimes the input file was in a subdirectory under the
	# build directory.  Keep only the input file for latter
	index = string.rfind(test.inputFile, "/")
	if index > 0:
	   test.inputFile = test.inputFile[index+1:]


        # if we are a "C_Src" build, we need the probin file
        if (suite.sourceTree == "C_Src" or \
                (test.testSrcTree == "C_Src" and test.probinFile != "")):
            try: shutil.copy(test.probinFile, outputDir)
            except IOError:
                errorMsg = "    ERROR: unable to copy probin file: %s" % test.probinFile
                reportTestFailure(errorMsg, test, testDir, full_web_dir)
                continue

            # sometimes the probin file was in a subdirectory under the
            # build directory.  Keep only the probin file for latter
            index = string.rfind(test.probinFile, "/")
            if index > 0:
               test.probinFile = test.probinFile[index+1:]


        # python doesn't allow labelled continue statements, so we
        # use skip_to_next_test to decide if we need to skip to
        # the next test
        skip_to_next_test = 0
        for file in test.auxFiles:
            try: shutil.copy(file, outputDir)
            except IOError:
                errorMsg = "    ERROR: unable to copy aux file: %s" % file
                reportTestFailure(errorMsg, test, testDir, full_web_dir)
                skip_to_next_test = 1
                break

        if skip_to_next_test: continue

        # python doesn't allow labelled continue statements, so we
        # use skip_to_next_test to decide if we need to skip to
        # the next test
        skip_to_next_test = 0
        for file in test.linkFiles:
            if not os.path.exists(file):
                errorMsg = "    ERROR: link file %s does not exist" % file
                reportTestFailure(errorMsg, test, testDir, full_web_dir)
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
                    reportTestFailure(errorMsg, test, testDir, full_web_dir)
                    skip_to_next_test = 1
                    break

        if skip_to_next_test: continue


        #----------------------------------------------------------------------
        # run the test
        #----------------------------------------------------------------------
        print "  running the test..."

        os.chdir(outputDir)

        test.wallTime = time.time()

        if suite.sourceTree == "C_Src" or test.testSrcTree == "C_Src":

            base_command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk" % \
                           (executable, test.inputFile, test.name, test.name)

            # keep around the checkpoint files only for the restart runs
            if test.restartTest:
                base_command += " amr.checkpoint_files_output=1 amr.check_int=%d" % \
                                (test.restartFileNum)
            else:
                base_command += " amr.checkpoint_files_output=0"

            base_command += " >& %s.run.out < /dev/null" % (test.name)

        elif suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src":

            base_command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk " % \
                           (executable, test.inputFile, test.name, test.name)

            # keep around the checkpoint files only for the restart runs
            if not test.restartTest: base_command += " --chk_int 0 "

            base_command += "%s >& %s.run.out" % \
                            (suite.globalAddToExecString, test.name)

        testRunCommand = suite.run_test(test, base_command)


        # if it is a restart test, then rename the final output file and
        # restart the test
        if test.restartTest:
            lastFile = getLastPlotfile(outputDir, test)

            if lastFile == "":
                errorMsg = "ERROR: test did not produce output.  Restart test not possible"
                reportTestFailure(errorMsg, test, testDir, full_web_dir)
                continue

            origLastFile = "orig_%s" % (lastFile)
            shutil.move(lastFile, origLastFile)

            if test.diffDir:
                origDiffDir = "orig_%s" % (test.diffDir)
                shutil.move(test.diffDir, origDiffDir)

            # get the file number to restart from
            restartFile = "%s_chk%5.5d" % (test.name, test.restartFileNum)

            print "  restarting from %s ... " % (restartFile)

            if suite.sourceTree == "C_Src" or test.testSrcTree == "C_Src":

                base_command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 amr.restart=%s >> %s.run.out 2>&1" % \
                        (executable, test.inputFile, test.name, test.name, restartFile, test.name)

            elif suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src":

                base_command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 --restart %d %s >> %s.run.out 2>&1" % \
                        (executable, test.inputFile, test.name, test.name, test.restartFileNum, suite.globalAddToExecString, test.name)

            testRunCommand = suite.run_test(test, base_command)

        test.wallTime = time.time() - test.wallTime


        #----------------------------------------------------------------------
        # do the comparison
        #----------------------------------------------------------------------
        if not test.selfTest:

            if test.outputFile == "":
                if test.compareFile == "":
                    compareFile = getLastPlotfile(outputDir, test)
                else:
                    compareFile = test.compareFile
                outputFile = compareFile
            else:
                outputFile = test.outputFile
                compareFile = test.name+'_'+outputFile


            # get the number of levels for reporting
            prog = "../fboxinfo.exe -l {}".format(outputFile)
            stdout0, stderr0, rc = run(prog)
            test.nlevels = stdout0.rstrip('\n')
            if not isInt(test.nlevels):
                test.nlevels = ""

            if args.make_benchmarks == None:

                print "  doing the comparison..."
                print "    comparison file: ", outputFile

                test.compare_file_used = outputFile

                if not test.restartTest:
                    benchFile = bench_dir + compareFile
                else:
                    benchFile = origLastFile

                # see if it exists
                # note, with BoxLib, the plotfiles are actually directories

                if not os.path.isdir(benchFile):
                    warning("    WARNING: no corresponding benchmark found")
                    benchFile = ""

                    cf = open("%s.compare.out" % (test.name), 'w')
                    cf.write("WARNING: no corresponding benchmark found\n")
                    cf.write("         unable to do a comparison\n")
                    cf.close()

                else:
                    if not compareFile == "":

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

                if not test.diffDir == "":
                    if not test.restartTest:
                        diffDirBench = bench_dir + '/' + test.name + '_' + test.diffDir
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

                    if diffstatus == 0:
                        cf = open("%s.compare.out" % (test.name), 'a')
                        cf.write("diff was SUCCESSFUL\n")
                        cf.close()


            else:   # make_benchmarks

                print "  storing output of %s as the new benchmark..." % (test.name)
                print "     new benchmark file: ", compareFile

                if not compareFile == "":
                    if not outputFile == compareFile:
                        source_file = outputFile
                    else:
                        source_file = compareFile

                    try: shutil.rmtree("{}/{}".format(bench_dir, compareFile))
                    except: pass
                    shutil.copytree(source_file, "{}/{}".format(bench_dir, compareFile))

                    cf = open("%s.status" % (test.name), 'w')
                    cf.write("benchmarks updated.  New file:  %s\n" % (compareFile) )
                    cf.close()
                else:
                    cf = open("%s.status" % (test.name), 'w')
                    cf.write("benchmarks failed")
                    cf.close()

                if not test.diffDir == "":
                    diffDirBench = "{}/{}_{}".format(bench_dir, test.name, test.diffDir)
                    if os.path.isdir(diffDirBench):
                        shutil.rmtree(diffDirBench)
                        shutil.copytree(test.diffDir, diffDirBench)
                    else:
                        shutil.copy(test.diffDir, diffDirBench)
                    print "     new diffDir: {}_{}".format(test.name, test.diffDir)

        else:   # selfTest

            if args.make_benchmarks == None:

                print "  looking for selfTest success string: %s ..." % test.stSuccessString

                try: of = open("%s.run.out" % (test.name), 'r')
                except IOError:
                    warning("WARNING: no output file found")
                    compareSuccessful = 0
                    outLines = ['']
                else:
                    outLines = of.readlines()

                    # successful comparison is indicated by PLOTFILES AGREE
                    compareSuccessful = 0

                    for line in outLines:
                        if line.find(test.stSuccessString) >= 0:
                            compareSuccessful = 1
                            break

                    of.close()

                cf = open("%s.compare.out" % (test.name), 'w')

                if compareSuccessful:
                    cf.write("SELF TEST SUCCESSFUL\n")
                else:
                    cf.write("SELF TEST FAILED\n")

                cf.close()



        #----------------------------------------------------------------------
        # do any requested visualization (2- and 3-d only)
        #----------------------------------------------------------------------
        if test.doVis and args.make_benchmarks == None:

            if not outputFile == "":

                print "  doing the visualization..."

                if test.dim == 2:
                    systemCall('%s/%s --palette %s/Palette -cname "%s" -p "%s" >& /dev/null' %
                              (suite.compareToolDir, vis2dExecutable, suite.compareToolDir,
                               test.visVar, outputFile) )
                elif test.dim == 3:
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
        if not test.analysisRoutine == "" and args.make_benchmarks == None:

            if not outputFile == "":

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
        if args.make_benchmarks == None:
            shutil.copy("%s.run.out"     % (test.name), full_web_dir)
            shutil.copy("%s.make.out"    % (test.name), full_web_dir)
            shutil.copy("%s.compare.out" % (test.name), full_web_dir)

            shutil.copy(test.inputFile, "%s/%s.%s" % (full_web_dir, test.name, test.inputFile) )

            if suite.sourceTree == "C_Src":
                shutil.copy(test.probinFile, "%s/%s.%s" % (full_web_dir, test.name, test.probinFile) )

            for file in test.auxFiles:

                # sometimes the auxFile was in a subdirectory under the
                # build directory.
                index = string.rfind(file, "/")
                if index > 0:
                    file = file[index+1:]

                shutil.copy(file, "%s/%s.%s" % (full_web_dir, test.name, file) )

            if test.doVis:
               png_file = getRecentFileName(outputDir, "", ".png")
               if not png_file is None:
                   try: shutil.copy(png_file, full_web_dir)
                   except IOError:
                       # visualization was not successful.  Reset doVis
                       test.doVis = 0

            if not test.analysisRoutine == "":
                try: shutil.copy(test.analysisOutputImage, full_web_dir)
                except IOError:
                    # analysis was not successful.  Reset the output image
                    test.analysisOutputImage = ""


        else:
            shutil.copy("%s.status" % (test.name), full_web_dir)


        #----------------------------------------------------------------------
        # archive (or delete) the output
        #----------------------------------------------------------------------
        print "  archiving the output..."
        for file in os.listdir(outputDir):
            if (os.path.isdir(file) and
                (file.startswith("%s_plt" % (test.name)) or
                 file.startswith("%s_chk" % (test.name)) ) ):

                if suite.purge_output == 1 and not file == outputFile:
                    # delete the plt/chk file
                    if os.path.isdir(file):
                        try: shutil.rmtree(file)
                        except:
                            warning("    WARNING: unable to remove {}".format(file))

                else:
                    # tar it up
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
        if args.make_benchmarks == None:
            print "  creating problem test report ..."
            reportSingleTest(suite, test, compString, testRunCommand, testDir, full_web_dir)


    #--------------------------------------------------------------------------
    # write the report for this instance of the test suite
    #--------------------------------------------------------------------------
    bold("creating new test report...", skip_before=1)
    numFailed = reportThisTestRun(suite, args.make_benchmarks, args.note,
                                  updateTime,  updateBoxLib,
                                  updateSource, updateExtSrc,
                                  testList, testDir, args.input_file[0], full_web_dir)


    # make sure that all of the files in the web directory are world readable
    for file in os.listdir(full_web_dir):
       currentFile = full_web_dir + file

       if os.path.isfile(currentFile):
          os.chmod(currentFile, 0644)

    if updateBoxLib or args.boxLibGitHash:
        doGITback(suite.boxLibDir, "BoxLib", boxLibGitBranch_Orig)

    if updateSource or args.sourceGitHash:
        doGITback(suite.sourceDir, suite.srcName, sourceGitBranch_Orig)

    if updateExtSrc or args.extSrcGitHash:
        doGITback(suite.extSrcDir, suite.extSrcName, extSrcGitBranch_Orig)

    #--------------------------------------------------------------------------
    # For temporary run, return now without creating suote report.
    #--------------------------------------------------------------------------
    if args.do_temp_run:
        return numFailed


    #--------------------------------------------------------------------------
    # generate the master report for all test instances
    #--------------------------------------------------------------------------
    bold("creating suite report...", skip_before=1)
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

    if numFailed > 0 and suite.sendEmailWhenFail:
        bold("sending email...", skip_before=1)
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

#summary th {background-color: grey;
    color: yellow;
    text-align: center;
    height: 2em;
    padding-bottom: 3px;
    padding-left: 5px;
    padding-right: 5px;}


#summary td {background: transparent;}

#summary tr:nth-child(even) {background: #dddddd;}
#summary tr:nth-child(odd) {background: #eeeeee;}

#summary tr.special {background: #ccccff;}
#summary td.highlight {color: red;}

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
<CENTER><H2>@SUBTITLE@</H2></CENTER>
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


class HTMLTable(object):
    def __init__(self, out_file, columns=1):
        self.hf = out_file
        self.columns = columns

    def start_table(self):
        self.hf.write("<div id='summary'>\n")
        self.hf.write("<p><table>\n")

    def header(self, header_list):
        n = len(header_list)
        line = "<tr>"+n*"<th>{}</th>"+"</tr>\n"
        self.hf.write(line.format(*header_list))

    def print_single_row(self, row):
        self.hf.write("<tr class='special'><td colspan={}>".format(self.columns)+row+"</td></tr>\n")

    def print_row(self, row_list, highlight=False):
        n = len(row_list)
        if highlight:
            line = "<tr>"+n*"<td class='highlight'>{}</td>"+"</tr>\n"
        else:
            line = "<tr>"+n*"<td>{}</td>"+"</tr>\n"
        self.hf.write(line.format(*row_list))

    def end_table(self):
        self.hf.write("</table>\n")
        self.hf.write("</div>\n")



#==============================================================================
# reportSingleTest
#==============================================================================
def reportSingleTest(suite, test, compileCommand, runCommand, testDir, full_web_dir):
    """ generate a single problem's test result page """

    # get the current directory
    currentDir = os.getcwd()

    # switch to the web directory and open the report file
    os.chdir(full_web_dir)

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
        # successful compilation be indicated by SUCCESS or
        # Nothing to be done for `all'.  Look for both
        compileSuccessful = 0

        for line in cf:
            if (line.find("SUCCESS") >= 0 or
                line.find("is up to date.") >= 0 or
                line.find("Nothing to be done") >= 0):
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
                if (line.find("PLOTFILES AGREE") >= 0 or
                    line.find("SELF TEST SUCCESSFUL") >= 0):
                    compareSuccessful = 1
                    break

            if compareSuccessful:
                if not test.diffDir == "":
                    compareSuccessful = 0
                    for line in diffLines:
                        if line.find("diff was SUCCESSFUL") >= 0:
                            compareSuccessful = 1
                            break

            cf.close()


    #--------------------------------------------------------------------------
    # write out the status file for this problem, with either
    # PASSED or FAILED
    #--------------------------------------------------------------------------
    status_file = "%s.status" % (test.name)
    sf = open(status_file, 'w')

    if (compileSuccessful and
        (test.compileTest or (not test.compileTest and compareSuccessful))):
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

    if not test.compileTest:

        if test.debug:
            hf.write("<p><b>Debug test</b>\n")
            hf.write("<p>&nbsp;\n")

        if test.useMPI:
            hf.write("<P><b>Parallel (MPI) Run</b><br>numprocs = %d\n" % (test.numprocs) )
            hf.write("<P>&nbsp;\n")


        if test.useOMP:
            hf.write("<P><b>OpenMP Run</b><br>numthreads = %d\n" % (test.numthreads) )
            hf.write("<P>&nbsp;\n")

        hf.write("<p><b>Execution Time</b> (seconds) = %f\n" % (test.wallTime))


        # is this a restart test?
        if test.restartTest:

            hf.write("<P><b>Restart Test</b><br>Job was run as normal and then restarted from checkpoint # %d, and the two final outputs were compared\n" % (test.restartFileNum) )

        hf.write("<P>&nbsp;\n")

        # write out the information about the test
        hf.write("<P><b>input file:</b> <A HREF=\"%s.%s\">%s</A>\n" %
                 (test.name, test.inputFile, test.inputFile) )


        if (suite.sourceTree == "C_Src"):
            hf.write("<P><b>probin file:</b> <A HREF=\"%s.%s\">%s</A>\n" %
                     (test.name, test.probinFile, test.probinFile) )


        for i, file in enumerate(test.auxFiles):
            # sometimes the auxFile was in a subdirectory under the
            # build directory.
            index = string.rfind(file, "/")
            if (index > 0):
                root_file = file[index+1:]
            else:
                root_file = file

            hf.write("<P><b>aux%dFile:</b> <A HREF=\"%s.%s\">%s</A>\n" %
                     (i+1, test.name, root_file, file) )


        hf.write("<P><b>dimensionality:</b> %s\n" % (test.dim) )

        hf.write("<P>&nbsp;\n")


    # write out the compilation report
    if compileSuccessful:
        hf.write("<P><H3 CLASS=\"passed\">Compilation Successful</H3></P>\n")
    else:
        hf.write("<P><H3 CLASS=\"failed\">Compilation Failed</H3></P>\n")

    hf.write("<P>Compliation command:<br>\n&nbsp; %s\n" % (compileCommand) )
    hf.write("<P><A HREF=\"%s.make.out\">make output</A>\n" % (test.name) )

    hf.write("<P>&nbsp;\n")


    if not test.compileTest:

        # write out the comparison report
        if compareSuccessful:
            hf.write("<P><H3 CLASS=\"passed\">Comparison Successful</H3></P>\n")
        else:
            hf.write("<P><H3 CLASS=\"failed\">Comparison Failed</H3></P>\n")

        hf.write("<P>Execution command:<br>\n&nbsp; %s\n" % (runCommand) )
        hf.write("<P><A HREF=\"%s.run.out\">execution output</A>\n" % (test.name) )

        hf.write("<P>&nbsp;\n")

        # parse the compare output and make an HTML table
        ht = HTMLTable(hf, columns=3)
        in_diff_region = False

        box_error = False
        for line in diffLines:

            if "number of boxes do not match" in line:
                box_error = True
                break
            
            if not in_diff_region:
                if line.find("fcompare") > 1:
                    hf.write("<pre>"+line+"</pre>\n")

                    ht.start_table()
                    continue

                if line.strip().startswith("diff"):
                    ht.end_table()
                    hf.write("<pre>\n")

                    hf.write(line.strip())
                    in_diff_region = True

                if line.strip().startswith("level"):
                    ht.print_single_row(line.strip())
                    continue

                if line.strip().startswith("-----"):
                    continue

                if line.strip().startswith("<<<"):
                    ht.print_single_row(line.strip().replace('<','&lt;').replace('>','&gt;'))
                    continue

                fields = [q.strip() for q in line.split("  ") if not q == ""]
    
                if fields[0].startswith("variable"):
                    ht.header(fields)
                    continue

                if len(fields) == 2:
                    ht.header([" "] + fields)
                    continue

                if len(fields) == 1:
                    continue

                else:
                    abs_err = float(fields[1])
                    rel_err = float(fields[2])
                    if abs(rel_err) > 1.e-6:
                        ht.print_row([fields[0], abs_err, rel_err], highlight=True)
                    else:
                        ht.print_row([fields[0], abs_err, rel_err])

            else:
                # diff region
                hf.write(line.strip())

        if in_diff_region:
            hf.write("</pre>\n")
        else:
            ht.end_table()

        if box_error:
            hf.write("<p>number of boxes do not match</p>\n")
        
        # show any visualizations
        if test.doVis:
            png_file = getRecentFileName(full_web_dir, test.name, ".png")
            if not png_file is None:
                hf.write("<P>&nbsp;\n")
                hf.write("<P><IMG SRC='%s' BORDER=0>" % (png_file) )

        # show any analysis
        if not test.analysisOutputImage == "":
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
def reportTestFailure(message, test, testDir, full_web_dir, compString=None):
    """ generate a simple report for an error encountered while performing
        the test """

    testfail("    aborting test")
    testfail(message)

    # get the current directory
    currentDir = os.getcwd()

    # switch to the web directory and open the report file
    os.chdir(full_web_dir)

    #--------------------------------------------------------------------------
    # write out the status file for this problem -- FAILED
    #--------------------------------------------------------------------------
    status_file = "%s.status" % (test.name)
    sf = open(status_file, 'w')
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
def reportThisTestRun(suite, make_benchmarks, note, updateTime,
                      updateBoxLib, updateSource, updateExtSrc,
                      testList, testDir, testFile, full_web_dir):
    """ generate the master page for a single run of the test suite """

    # get the current directory
    currentDir = os.getcwd()

    # switch to the web directory and open the report file
    os.chdir(full_web_dir)


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

    if not make_benchmarks == None:
       hf.write("<p><b>Benchmarks updated</b><br>comment: <font color=\"gray\">{}</font>\n".format(make_benchmarks) )
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

    if make_benchmarks == None:
        ht = HTMLTable(hf, columns=11)
        ht.start_table()
        ht.header(["test name", "dim", "compare plotfile",
                   "# levels", "MPI (# procs)", "OMP (# threads)", "debug?", 
                   "compile?", "restart?", "wall time", "result"])

    else:
        ht = HTMLTable(hf, columns=3)
        ht.start_table()
        ht.header(["test name", "result", "comment"])

    # loop over the tests and add a line for each
    for test in testList:

        if make_benchmarks == None:

            # check if it passed or failed
            status_file = "%s.status" % (test.name)

            sf = open(status_file, 'r')

            testPassed = 0

            for line in sf:
                if line.find("PASSED") >= 0:
                    testPassed = 1
                    numPassed += 1
                    break

            if not testPassed:
                numFailed += 1


            sf.close()

            row_info = []
            row_info.append("<a href=\"%s.html\">{}</a>".format(test.name, test.name))
            row_info.append(test.dim)
            row_info.append(test.compare_file_used)

            if not test.nlevels == None:
                row_info.append(test.nlevels)
            else:
                row_info.append("")

            if test.useMPI:
                row_info.append("&check; ({})".format(test.numprocs))
            else:
                row_info.append("")

            # OMP ?
            if test.useOMP:
                row_info.append("&check; ({})".format(test.numthreads))
            else:
                row_info.append("")

            # debug ?
            if test.debug:
                row_info.append("&check;")
            else:
                row_info.append("")

            # compile ?
            if test.compileTest:
                row_info.append("&check;")
            else:
                row_info.append("")

            # restart ?
            if test.restartTest:
                row_info.append("&check;")
            else:
                row_info.append("")

            # wallclock time
            row_info.append("{:.3f} s".format(test.wallTime))

            if testPassed:
                row_info.append("<h3 class=\"passed\">PASSED</h3>")
            else:
                row_info.append("<h3 class=\"failed\">FAILED</h3>")

            ht.print_row(row_info)

        else:
            if test.restartTest: continue
            if test.compileTest: continue
            if test.selfTest: continue

            # the benchmark was updated -- find the name of the new benchmark file
            benchStatusFile = "%s.status" % (test.name)

            bf = open(benchStatusFile, 'r')

            benchFile = "none"

            for line in bf:
                index = line.find("file:")
                if index >= 0:
                    benchFile = line[index+5:]
                    break

            row_info = []
            row_info.append("{}".format(test.name))
            if not benchFile == "none":
                row_info.append("<h3 class=\"benchmade\">BENCHMARK UPDATED</h3>")
                row_info.append("(new benchmark file is {})".format(benchFile))
            else:
                row_info.append("<h3 class=\"failed\">BENCHMARK NOT UPDATED</h3>")
                row_info.append("(compilation or execution failed)")

            ht.print_row(row_info)

    ht.end_table()

    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")
    hf.close()


    #--------------------------------------------------------------------------
    # write out a status file for all the tests
    #--------------------------------------------------------------------------

    index = testDir.find("/")
    if index > 0:
        status_file = testDir[0:index] + ".status"
    else:
        status_file = testDir + ".status"

    sf = open(status_file, 'w')

    if make_benchmarks == None:
        if numFailed == 0:
            sf.write("ALL PASSED\n")
        elif numFailed > 0 and numPassed > 0:
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
        if file.startswith("20") and os.path.isdir(file):

            # look for the status file
            status_file = file + '/' + file + '.status'

            if os.path.isfile(status_file):
                validDirs.append(file)


    validDirs.sort()
    validDirs.reverse()


    #--------------------------------------------------------------------------
    # now find all of the unique problems in the test directories
    #--------------------------------------------------------------------------
    for dir in validDirs:

        for file in os.listdir(webTopDir + dir):

            if file.endswith(".status") and not file.startswith("20"):

                index = string.rfind(file, ".status")
                testName = file[0:index]

                if allTests.count(testName) == 0:
                    if (not suite.reportActiveTestsOnly) or (testName in activeTestList):
                        allTests.append(testName)


    allTests.sort()


    #--------------------------------------------------------------------------
    # generate the HTML
    #--------------------------------------------------------------------------
    htmlFile = "index.html"

    title = "%s regression tests" % (suite.suiteName)

    hf = open(htmlFile, 'w')

    header = MainHeader.replace("@TITLE@", title).replace("@SUBTITLE@", suite.sub_title)
    
    if suite.goUpLink:
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
            status_file = "%s/%s/%s.status" % (webTopDir, dir, test)
            if os.path.isfile(status_file):
                valid = 1
                break

        if not valid: continue

        # write out the directory (date)
        hf.write("<TR><TD class='date'><SPAN CLASS='nobreak'><A class='main' HREF=\"%s/index.html\">%s&nbsp;</A></SPAN></TD>\n" %
                 (dir, dir) )

        for test in allTests:

            # look to see if the current test was part of this suite run
            status_file = "%s/%s/%s.status" % (webTopDir, dir, test)
            status = 0

            if os.path.isfile(status_file):

                sf = open(status_file, 'r')

                # status = -1 (failed); 1 (passed); 10 (benchmark update)
                status = -1
                for line in sf:
                    if line.find("PASSED") >= 0:
                        status = 1
                        break
                    elif line.find("FAILED") >= 0:
                        status = -1
                        break
                    elif line.find("benchmarks updated") >= 0:
                        status = 10
                        break

                sf.close()


            # write out this test's status
            if status == 1:
                hf.write("<TD ALIGN=CENTER title=\"%s\" class=\"passed\"><H3><a href=\"%s/%s.html\" class=\"passed\">:)</a></H3></TD>\n" % (test, dir, test))
            elif status == -1:
                hf.write("<TD ALIGN=CENTER title=\"%s\" class=\"failed\"><H3><a href=\"%s/%s.html\" class=\"failed\">&nbsp;!&nbsp;</a></H3></TD>\n" % (test, dir, test))
            elif status == 10:
                hf.write("<TD ALIGN=CENTER title=\"%s\" class=\"benchmade\"><H3>U</H3></TD>\n" % (test))
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
