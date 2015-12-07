#!/usr/bin/env python

"""
A simple regression test framework for a BoxLib-based code

There are several major sections to this source: the runtime parameter
routines, the test suite routines, and the report generation routines.
They are separated as such in this file.

This test framework understands source based out of the F_Src and
C_Src BoxLib frameworks.

"""

from __future__ import print_function

try: import ConfigParser as configparser
except ImportError:
    import configparser   # python 3

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
class Test(object):

    def __init__ (self, name):

        self.name = name

        self.log = None

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

        self.comp_string = None  # set automatically
        self.run_command = None  # set automatically
        
    def __lt__(self, other):
        return self.value() < other.value()

    def value(self):
        return self.name

    def get_last_plotfile(self, output_dir=None):
        """ Find the last plotfile written.  Note: we give an error if the
            last plotfile is 0.  If output_dir is specified, then we use
            that instead of the default
        """

        if output_dir is None:
            output_dir = self.output_dir   # not yet implemented

        plts = [d for d in os.listdir(output_dir) if \
                (os.path.isdir(d) and d.startswith("{}_plt".format(self.name))) or \
                (os.path.isfile(d) and d.startswith("{}_plt".format(self.name)) and d.endswith(".tgz"))]
    
        if len(plts) == 0:
            self.log.warn("WARNING: test did not produce any output")
            return ""

        plts.sort()
        last_plot = plts.pop()

        if last_plot.endswith("00000"):
            self.log.warn("WARNING: only plotfile 0 was output -- skipping comparison")
            return ""

        return last_plot


class Suite(object):

    def __init__ (self, args):

        self.args = args
        
        self.repos = {}

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

        self.boxLibGitBranch = "development"
        self.sourceGitBranch = "development"
        self.extSrcGitBranch = "development"

        # this will hold the # of extra build directories we need to worry
        # about.  It is set automatically, not by users
        self.useExtraBuild = 0

        self.extra_build_dirs = []
        self.extra_build_branches = []

        # this should be the environment variable name that should be
        # set so the builds in the extraBuildDir can see the main
        # source.  This environment variable will be set to the
        # sourceTree path and included on the make lines
        self.extraBuildDirCompString = ""

        # these are set automatically by the script -- they hold the
        # basename of the various source directories
        self.srcName = ""
        self.extSrcName = ""
        self.extra_build_names = []

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

        self.log = None

    def check_test_dir(self, dir_name):
        """ given a string representing a directory, check if it points to
            a valid directory.  If so, return the directory name """

        dir_name = os.path.normpath(dir_name) + "/"

        if not os.path.isdir(dir_name):
            self.log.fail("ERROR: {} is not a valid directory".format(dir_name))

        return dir_name

    def get_tests_to_run(self, test_list_old):
        """ perform various tests based on the runtime options to determine
            which of the tests in the input file we run """

        # if we only want to run the tests that failed previously,
        # remove the others
        if self.args.redo_failed or not self.args.copy_benchmarks is None:
            last_run = self.get_last_run()
            failed = self.get_test_failures(last_run)

            test_list = [t for t in test_list_old if t.name in failed]        
        else:
            test_list = test_list_old[:]

        # if we only want to run tests of a certain dimensionality, remove
        # the others
        if self.args.d in [1, 2, 3]:
            test_list = [t for t in test_list_old if t.dim == self.args.d]

        # if we are doing a single test, remove all other tests; if we
        # specified a list of tests, check each one; if we did both
        # --single_test and --tests, complain
        if not self.args.single_test == "" and not self.args.tests == "":
            self.log.fail("ERROR: specify tests either by --single_test or --tests, not both")

        if not self.args.single_test == "":
            tests_find = [self.args.single_test]
        elif not self.args.tests == "":
            tests_find = self.args.tests.split()
        else:
            tests_find = []

        if len(tests_find) > 0:
            new_test_list = []
            for test in tests_find:
                _tmp = [o for o in test_list if o.name == test]
                if len(_tmp) == 1:
                    new_test_list += _tmp
                else:
                    self.log.fail("ERROR: {} is not a valid test".format(test))

            test_list = new_test_list

        if len(test_list) == 0:
            self.log.fail("No valid tests defined")

        return test_list

    def get_bench_dir(self):
        bench_dir = self.testTopDir + self.suiteName + "-benchmarks/"
        if not os.path.isdir(bench_dir):
            if not self.args.make_benchmarks == None:
                os.mkdir(bench_dir)
            else:
                self.log.fail("ERROR: benchmark directory, %s, does not exist" % (bench_dir))
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
            if os.path.isdir(full_test_dir):
                shutil.rmtree(full_test_dir)
        else:
            for i in range(1, maxRuns):
                if not os.path.isdir(full_test_dir): break
                test_dir = today + "-{:03d}/".format(i)
                full_test_dir = self.testTopDir + self.suiteName + "-tests/" + test_dir

        self.log.skip()
        self.log.bold("testing directory is: " + test_dir)
        os.mkdir(full_test_dir)

        # make the web directory -- this is where all the output and HTML will be
        # put, so it is easy to move the entire test website to a different disk
        full_web_dir = "%s/%s/"  % (self.webTopDir, test_dir)

        if self.args.do_temp_run:
            if os.path.isdir(full_web_dir):
                shutil.rmtree(full_web_dir)

        os.mkdir(full_web_dir)

        # copy the test file into the web output directory
        shutil.copy(self.test_file_path, full_web_dir)

        self.test_dir = test_dir
        self.full_test_dir = full_test_dir
        self.full_web_dir = full_web_dir

    def get_last_run(self):
        """ return the name of the directory corresponding to the previous
            run of the test suite """

        outdir = self.testTopDir + self.suiteName + "-tests/"

        # this will work through 2099
        if os.path.isdir(outdir):
            dirs = [d for d in os.listdir(outdir) if (os.path.isdir(outdir + d) and
                                                      d.startswith("20"))]
            dirs.sort()

            return dirs[-1]
        else:
            return None

    def get_test_failures(self, test_dir):
        """ look at the test run in test_dir and return the list of tests that
            failed """

        cwd = os.getcwd()

        outdir = self.testTopDir + self.suiteName + "-tests/"

        os.chdir(outdir + test_dir)

        failed = []

        for test in os.listdir("."):
            if not os.path.isdir(test): continue

            # the status files are in the web dir
            status_file = "{}/{}/{}.status".format(self.webTopDir, test_dir, test)
            with open(status_file, "r") as sf:
                for line in sf:
                    if line.find("FAILED") >= 0:
                        failed.append(test)

        os.chdir(cwd)
        return failed

    def make_realclean(self):
        cmd = "{} BOXLIB_HOME={} {} {} realclean".format(
            self.MAKE, self.boxLibDir,
            self.extSrcCompString, self.extraBuildDirCompString)
        run(cmd)

    def build_f(self, opts="", target="", outfile=None):
        comp_string = "{} -j{} BOXLIB_HOME={} COMP={} {} {}".format(
            self.MAKE, self.numMakeJobs, self.boxLibDir, self.FCOMP, opts, target)
        self.log.log(comp_string)
        run(comp_string, outfile=outfile)
        return comp_string

    def run_test(self, test, base_command):
        test_env = None
        if test.useOMP:
            test_env = dict(os.environ, OMP_NUM_THREADS="{}".format(test.numthreads))

        if test.useMPI:
            test_run_command = self.MPIcommand
            test_run_command = test_run_command.replace("@host@", self.MPIhost)
            test_run_command = test_run_command.replace("@nprocs@", "{}".format(test.numprocs))
            test_run_command = test_run_command.replace("@command@", base_command)
        else:
            test_run_command = base_command

        self.log.log(test_run_command)
        sout, serr, ierr = run(test_run_command, stdin=True, outfile="{}.run.out".format(test.name), env=test_env)
        test.run_command = test_run_command

    def build_tools(self, test_list):

        self.compare_tool_dir = "{}/Tools/Postprocessing/F_Src/".format(
            os.path.normpath(self.boxLibDir))

        os.chdir(self.compare_tool_dir)

        self.make_realclean()

        tools = ["fcompare", "fboxinfo"]
        if any([t for t in test_list if t.dim == 2]): tools.append("fsnapshot2d")
        if any([t for t in test_list if t.dim == 3]): tools.append("fsnapshot3d")

        self.tools = {}

        self.log.skip()
        self.log.bold("building tools...")
        self.log.indent()

        for t in tools:
            self.log.log("building {}...".format(t))
            self.build_f(target="programs={}".format(t), opts="NDEBUG=t MPI= ")
            exe = get_recent_filename(self.compare_tool_dir, t, ".exe")
            self.tools[t] = "{}/{}".format(self.compare_tool_dir, exe)

        self.log.outdent()


class Repo(object):
    """ a simple class to manage our git operations """
    def __init__(self, suite, directory, name, 
                 branch_wanted=None, hash_wanted=None, update=True):
        self.suite = suite
        self.dir = directory
        self.name = name
        self.branch_wanted = branch_wanted
        self.branch_orig = None
        self.hash_wanted = hash_wanted
        self.hash_current = None

        self.update = update
        if hash_wanted:
            self.update = False

    def git_update(self):
        """ Do a git update of the repository.  If githash is not empty, then
            we will check out that version instead of git-pulling. """

        os.chdir(self.dir)

        # find out current branch so that we can go back later if we need.
        stdout0, stderr0, rc = run("git rev-parse --abbrev-ref HEAD")
        self.branch_orig = stdout0.rstrip('\n')

        if self.branch_orig != self.branch_wanted:
            self.suite.log.log("git checkout {} in {}".format(self.branch_wanted, self.dir))
            stdout, stderr, rc = run("git checkout {}".format(self.branch_wanted), 
                                     stdin=True)
        else:
            self.branch_wanted = self.branch_orig
            
        if self.hash_wanted == "" or self.hash_wanted == None:
            self.suite.log.log("'git pull' in {}".format(self.dir))

            # we need to be tricky here to make sure that the stdin is
            # presented to the user to get the password.
            stdout, stderr, rc = run("git pull", stdin=True,
                                     outfile="git.{}.out".format(self.name))

        else:
            stdout, stderr, rc = run("git checkout {}".format(self.hash_wanted),
                                     outfile="git.{}.out".format(self.name))

        if not rc == 0:
            self.suite.log.fail("ERROR: git update was unsuccessful")

        shutil.copy("git.{}.out".format(self.name), self.suite.full_web_dir)

    def save_head(self):

        os.chdir(self.dir)

        self.suite.log.log("saving git HEAD for {}/".format(self.name))

        stdout, stderr, rc = run("git rev-parse HEAD", 
                                 outfile="git.{}.HEAD".format(self.name) )

        self.hash_current = stdout
        shutil.copy("git.{}.HEAD".format(self.name), self.suite.full_web_dir)

    def make_changelog(self):
        """ generate a ChangeLog git repository, and copy it to the
            web directory"""

        os.chdir(self.dir)

        self.suite.log.log("generating ChangeLog for {}/".format(self.name))

        run("git log --name-only", outfile="ChangeLog.{}".format(self.name), outfile_mode="w")
        shutil.copy("ChangeLog.{}".format(self.name), self.suite.full_web_dir)

    def git_back(self):
        """ switch the repo back to its original branch """

        os.chdir(self.dir)
        self.suite.log.log("git checkout {} in {}".format(self.branch_orig, self.dir))

        stdout, stderr, rc = run("git checkout {}".format(self.branch_orig), 
                                 stdin=True,
                                 outfile="git.{}.out".format(self.name))

        if not rc == 0:
            self.suite.log.fail("ERROR: git checkout was unsuccessful")


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# R U N T I M E   P A R A M E T E R   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

def convert_type(string):
    """ return an integer, float, or string from the input string """
    try: int(string)
    except: pass
    else: return int(string)

    try: float(string)
    except: pass
    else: return float(string)

    return string.strip()    

def load_params(args):
    """
    reads the parameter file and creates as list of test objects as well as
    the suite object
    """

    testList = []

    cp = configparser.ConfigParser()    # note, need strict=False for Python3
    cp.optionxform = str

    log = Log()

    log.bold("loading " + args.input_file[0])
    
    try: cp.read(args.input_file[0])
    except:
        log.fail("ERROR: unable to read parameter file {}".format(file))

    # "main" is a special section containing the global suite parameters.
    mysuite = Suite(args)

    mysuite.log = log

    valid_options = list(mysuite.__dict__.keys())
    valid_options += ["extraBuildDir", "extraBuildDir2"]

    for opt in cp.options("main"):

        # get the value of the current option
        value = convert_type(cp.get("main", opt))

        if opt in valid_options:

            if opt == "sourceTree":
                if not value in ["C_Src", "F_Src", "BoxLib"]:
                    mysuite.log.fail("ERROR: invalid sourceTree")
                else:
                    mysuite.sourceTree = value

            elif opt == "sourceDir": mysuite.sourceDir = mysuite.check_test_dir(value)
            elif opt == "boxLibDir": mysuite.boxLibDir = mysuite.check_test_dir(value)
            elif opt == "testTopDir": mysuite.testTopDir = mysuite.check_test_dir(value)
            elif opt == "webTopDir": mysuite.webTopDir = os.path.normpath(value) + "/"
            elif opt == "extSrcDir":
                mysuite.extSrcDir = mysuite.check_test_dir(value)
                mysuite.useExtSrc = 1

            elif opt == "extraBuildDir":
                # we will keep the extra build directories in a list -- we
                # want them in the correct order in the list so we can simply
                # index later

                try: branch = convert_type(cp.get("main", "extraBuildBranch"))
                except: branch = "master"
                
                if len(mysuite.extra_build_dirs) == 0:
                    mysuite.extra_build_dirs.append(mysuite.check_test_dir(value))
                    mysuite.extra_build_branches.append(branch)
                else:
                    mysuite.extra_build_dirs.insert(0, mysuite.check_test_dir(value))
                    mysuite.extra_build_branches.insert(0, branch)

            elif opt == "extraBuildDir2":
                try: branch = convert_type(cp.get("main", "extraBuildBranch2"))
                except: branch = "master"

                mysuite.extra_build_dirs.append(mysuite.check_test_dir(value))
                mysuite.extra_build_branches.append(branch)
                
            elif opt == "emailTo": mysuite.emailTo = value.split(",")

            else:
                # generic setting of the object attribute
                setattr(mysuite, opt, value)

        else:
            mysuite.log.warn("WARNING: suite parameter %s not valid" % (opt))


    mysuite.useExtraBuild = len(mysuite.extra_build_dirs)

    # create the repo objects
    mysuite.repos["BoxLib"] = Repo(mysuite, mysuite.boxLibDir, "BoxLib",
                                   branch_wanted=mysuite.boxLibGitBranch, 
                                   hash_wanted=mysuite.args.boxLibGitHash)

    mysuite.srcName = os.path.basename(os.path.normpath(mysuite.sourceDir))
    if not mysuite.sourceTree == "BoxLib":
        mysuite.repos["source"] = Repo(mysuite, mysuite.sourceDir, mysuite.srcName,
                                       branch_wanted=mysuite.sourceGitBranch, 
                                       hash_wanted=mysuite.args.sourceGitHash)

    if mysuite.useExtSrc:
        mysuite.extSrcName = os.path.basename(os.path.normpath(mysuite.extSrcDir))

        # update the additional compilation str for additional source dir
        if mysuite.extSrcCompString != "":
            mysuite.extSrcCompString += "="+mysuite.extSrcDir

        mysuite.repos["extra_source"] = Repo(mysuite, mysuite.extSrcDir, 
                                             mysuite.extSrcName,
                                             branch_wanted=mysuite.extSrcGitBranch,
                                             hash_wanted=mysuite.args.extSrcGitHash)

    # update additional compiled string for any extra build directory
    if mysuite.useExtraBuild > 0:
        for n in range(len(mysuite.extra_build_dirs)):
            extra_build_name = os.path.basename(os.path.normpath(mysuite.extra_build_dirs[n]))
            mysuite.extra_build_names.append(extra_build_name)

            mysuite.repos["extra_build-{}".format(n)] = \
                Repo(mysuite, mysuite.extra_build_dirs[n],
                     extra_build_name, 
                     branch_wanted=mysuite.extra_build_branches[n])

        # since we are building in the extraBuildDir, we need to
        # tell make where the sourceDir is
        if mysuite.extraBuildDirCompString != "":
            mysuite.extraBuildDirCompString += "="+mysuite.sourceDir

    # BoxLib-only tests don't have a sourceDir
    if mysuite.sourceTree == "BoxLib": mysuite.sourceDir = mysuite.boxLibDir

    # checks
    if mysuite.sendEmailWhenFail and not args.send_no_email:
        if mysuite.emailTo == [] or mysuite.emailBody == "":
            mysuite.log.fail("ERROR: when sendEmailWhenFail = 1, you must specify emailTo and emailBody\n")

        if mysuite.emailFrom == "":
            mysuite.emailFrom = '@'.join((getpass.getuser(), socket.getfqdn()))

        if mysuite.emailSubject == "":
            mysuite.emailSubject = mysuite.suiteName+" Regression Test Failed"

    if (mysuite.sourceTree == "" or mysuite.boxLibDir == "" or
        mysuite.sourceDir == "" or mysuite.testTopDir == ""):
        mysuite.log.fail("ERROR: required suite-wide directory not specified\n" + \
                         "(sourceTree, boxLibDir, sourceDir, testTopDir)")

    # Make sure the web dir is valid (or use the default is none specified)
    if mysuite.webTopDir == "":
        mysuite.webTopDir = "{}/{}-web/".format(mysuite.testTopDir, mysuite.suiteName)

    if not os.path.isdir(mysuite.webTopDir):
        try: os.mkdir(mysuite.webTopDir)
        except: 
            mysuite.log.fail("ERROR: unable to create the web directory: {}\n".format(
                mysuite.webTopDir))

    # all other sections are tests
    mysuite.log.skip()
    mysuite.log.bold("finding tests and checking parameters...")

    for sec in cp.sections():

        if sec == "main": continue

        # maximum test name length -- used for HTML formatting
        mysuite.lenTestName = max(mysuite.lenTestName, len(sec))

        # create the test object for this test
        mytest = Test(sec)
        mytest.log = log
        invalid = 0

        # set the test object data by looking at all the options in
        # the current section of the parameter file
        valid_options = list(mytest.__dict__.keys())
        valid_options += ["aux1File", "aux2File", "aux3File"]
        valid_options += ["link1File", "link2File", "link3File"]

        for opt in cp.options(sec):

            # get the value of the current option
            value = convert_type(cp.get(sec, opt))

            if opt in valid_options:

                if opt in ["aux1File", "aux2File", "aux3File"]:
                    mytest.auxFiles.append(value)

                elif opt in ["link1File", "link2File", "link3File"]:
                    mytest.linkFiles.append(value)

                else:
                    # generic setting of the object attribute
                    setattr(mytest, opt, value)

            else:
                mysuite.log.warn("WARNING: unrecognized parameter {} for test {}".format(opt, sec))


        # make sure that the build directory actually exists
        if mytest.useExtraBuildDir > 0:
            bDir = mysuite.extra_build_dirs[mytest.useExtraBuildDir-1] + mytest.buildDir
        else:
            bDir = mysuite.sourceDir + mytest.buildDir

        if not os.path.isdir(bDir):
            mysuite.log.warn("WARNING: invalid build directory: {}".format(bDir))
            invalid = 1


        # make sure all the require parameters are present
        if mytest.compileTest:
            if mytest.buildDir == "":
                mysuite.log.warn("WARNING: mandatory parameters for test {} not set".format(sec))
                invalid = 1

        else:
            if (mytest.buildDir == "" or mytest.inputFile == "" or
                (mysuite.sourceTree == "C_Src" and mytest.probinFile == "") or
                mytest.dim == -1):
                mysuite.log.warn("WARNING: mandatory parameters for test {} not set".format(sec))
                mysuite.log.warn("         buildDir = {}".format(mytest.buildDir))
                mysuite.log.warn("         inputFile = {}".format(mytest.inputFile))
                if mysuite.sourceTree == "C_Src":
                    mysuite.log.warn("         probinFile = {}".format(mytest.probinFile))
                mysuite.log.warn("            dim = {}".format(mytest.dim))

                invalid = 1

        # check the optional parameters
        if mytest.restartTest and mytest.restartFileNum == -1:
            mysuite.log.warn("WARNING: restart-test {} needs a restartFileNum".format(sec))
            invalid = 1

        if mytest.selfTest and mytest.stSuccessString == "":
            mysuite.log.warn("WARNING: self-test {} needs a stSuccessString".format(sec))
            invalid = 1

        if mytest.useMPI and mytest.numprocs == -1:
            mysuite.log.warn("WARNING: MPI parallel test {} needs numprocs".format(sec))
            invalid = 1

        if mytest.useOMP and mytest.numthreads == -1:
            mysuite.log.warn("WARNING: OpenMP parallel test {} needs numthreads".format(sec))
            invalid = 1

        if mytest.doVis and mytest.visVar == "":
            mysuite.log.warn("WARNING: test {} has visualization, needs visVar".format(sec))
            invalid = 1

        if mysuite.sourceTree == "BoxLib" and mytest.testSrcTree == "":
            mysuite.log.warn("WARNING: test {} is a BoxLib test but testSrcTree not set".format(sec))
            invalid = 1


        # add the current test object to the master list
        if not invalid:
            testList.append(mytest)
        else:
            mysuite.log.warn("WARNING: test {} will be skipped".format(sec))


    # if any runs are parallel, make sure that the MPIcommand is defined
    anyMPI = any([t.useMPI for t in testList])

    if anyMPI and mysuite.MPIcommand == "":
        mysuite.log.fail("ERROR: some tests are MPI parallel, but MPIcommand not defined")

    testList.sort()

    return mysuite, testList


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# O U T P U T   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class Log(object):
    def __init__(self, output_file=None):

        # http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python
        # which in-turn cites the blender build scripts
        self.warn_color = '\033[33m'
        self.success_color = '\033[32m'
        self.fail_color = '\033[31m'
        self.bold_color = '\033[1m'
        self.end_color = '\033[0m'

        self.current_indent = 0
        self.indent_str = ""

        if not output_file is None:
            self.of = output_file
        else:
            self.of = None
            
    def indent(self):
        self.current_indent += 1
        self.indent_str = self.current_indent*"   "

    def outdent(self):
        self.current_indent -= 1
        self.current_indent = max(0, self.current_indent)
        self.indent_str = self.current_indent*"   "

    def fail(self, string):
        nstr = self.fail_color + string + self.end_color
        print("{}{}".format(self.indent_str, nstr))
        self.close_log()
        sys.exit()

    def testfail(self, string):
        nstr = self.fail_color + string + self.end_color
        print("{}{}".format(self.indent_str, nstr))

    def warn(self, string):
        nstr = self.warn_color + string + self.end_color
        print("{}{}".format(self.indent_str, nstr))

    def success(self, string):
        nstr = self.success_color + string + self.end_color
        print("{}{}".format(self.indent_str, nstr))

    def log(self, string):
        print("{}{}".format(self.indent_str, string))

    def skip(self):
        print("")

    def bold(self, string):
        nstr = self.bold_color + string + self.end_color
        print("{}{}".format(self.indent_str, nstr))

    def close_log(self):
        if not self.of is None: self.of.close()
    
        
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# S Y S T E M   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def get_args(arg_string=None):
    """ parse the commandline arguments.  If arg_string is present, we
        parse from there, otherwise we use the default (sys.argv) """

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", type=int, default=-1,
                        help="restrict tests to a particular dimensionality")
    parser.add_argument("--make_benchmarks", type=str, default=None, metavar="comment",
                        help="make new benchmarks? (must provide a comment)")
    parser.add_argument("--copy_benchmarks", type=str, default=None, metavar="comment",
                        help="copy the last plotfiles from the failed tests of the most recent run as the new benchmarks.  No git pull is done and no new runs are performed (must provide a comment)")
    parser.add_argument("--no_update", type=str, default="None", metavar="name",
                        help="which codes to exclude from the git update? (None, All, or a comma-separated list of codes)")
    parser.add_argument("--single_test", type=str, default="", metavar="test-name",
                        help="name of a single test to run")
    parser.add_argument("--tests", type=str, default="", metavar="'test1 test2 test3'",
                        help="a space-separated list of tests to run")
    parser.add_argument("--do_temp_run", action="store_true",
                        help="is this a temporary run? (output not stored or logged)")
    parser.add_argument("--send_no_email", action="store_true",
                        help="do not send emails when tests fail")
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

    if not arg_string is None:
        args = parser.parse_args(arg_string)
    else:
        args = parser.parse_args()

    return args


def run(string, stdin=False, outfile=None, store_command=False, env=None, outfile_mode="a", log=None):

    # shlex.split will preserve inner quotes
    prog = shlex.split(string)
    if stdin:
        p0 = subprocess.Popen(prog, stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, env=env)
    else:
        p0 = subprocess.Popen(prog, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, env=env)

    stdout0, stderr0 = p0.communicate()
    if stdin: p0.stdin.close()
    rc = p0.returncode
    p0.stdout.close()

    if not outfile == None:
        try: cf = open(outfile, outfile_mode)
        except IOError:
            log.fail("  ERROR: unable to open file for writing")
        else:
            if store_command:
                cf.write(string)
            for line in stdout0:
                cf.write(line)
            cf.close()

    return stdout0, stderr0, rc


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# T E S T   S U I T E   R O U T I N E S
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

def find_build_dirs(tests):
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

def copy_benchmarks(old_full_test_dir, full_web_dir, test_list, bench_dir, log):
    """ copy the last plotfile output from each test in testList 
        into the benchmark directory.  Also copy the diffDir, if
        it exists """
    td = os.getcwd()
    
    for t in test_list:
        wd = "{}/{}".format(old_full_test_dir, t.name)
        os.chdir(wd)

        if t.compareFile == "":
            p = t.get_last_plotfile(output_dir=wd)
        else:
            if not os.path.isdir(t.compareFile):
                p = get_recent_filename(wd, t.compareFile, ".tgz")
            else:
                p = t.compareFile
            
        if not p == "": 
            if p.endswith(".tgz"):
                try:
                    tg = tarfile.open(name=p, mode="r:gz")
                    tg.extractall()
                except:
                    log.fail("ERROR extracting tarfile")
                idx = p.rfind(".tgz")
                p = p[:idx]

            try: shutil.rmtree("{}/{}".format(bench_dir, p))
            except: pass
            shutil.copytree(p, "{}/{}".format(bench_dir, p))

            with open("{}/{}.status".format(full_web_dir, t.name), 'w') as cf:
                cf.write("benchmarks updated.  New file:  {}\n".format(p) )

        else:   # no benchmark exists
            with open("{}/{}.status".format(full_web_dir, t.name), 'w') as cf:
                cf.write("benchmarks update failed")
            
        # is there a diffDir to copy too?
        if not t.diffDir == "":
            diff_dir_bench = "{}/{}_{}".format(bench_dir, t.name, t.diffDir)
            if os.path.isdir(diff_dir_bench):
                shutil.rmtree(diff_dir_bench)
                shutil.copytree(t.diffDir, diff_dir_bench)
            else:
                shutil.copy(t.diffDir, diff_dir_bench)
            log.log("new diffDir: {}_{}".format(t.name, t.diffDir))
            
        os.chdir(td)
        
def get_recent_filename(dir, base, extension):
    """ find the most recent file matching the base and extension """

    files = [f for f in os.listdir(dir) if (f.startswith(base) and
                                            f.endswith(extension))]

    files.sort(key=lambda x: os.path.getmtime(x))

    try: return files.pop()
    except: return None

def convert_to_f_make_flag(opt, test_not=False):
    if test_not:
        if opt: return " "
        else: return "t"
    else:
        if opt: return "t"
        else: return " "


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# test
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def test_suite(argv):

    usage = """
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
    args=get_args(arg_string=argv)

    #--------------------------------------------------------------------------
    # read in the test information
    #--------------------------------------------------------------------------
    suite, testList = load_params(args)

    activeTestList = [t.name for t in testList]

    testList = suite.get_tests_to_run(testList)

    suite.log.skip()
    suite.log.bold("running tests: ")
    suite.log.indent()
    for obj in testList:
        suite.log.log(obj.name)
    suite.log.outdent()

    if not args.complete_report_from_crash == "":

        # make sure the web directory from the crash run exists
        suite.full_web_dir = "%s/%s/"  % (suite.webTopDir, args.complete_report_from_crash)
        if not os.path.isdir(suite.full_web_dir):
            suite.log.fail("Crash directory does not exist")

        # find all the tests that completed in that web directory
        tests = []
        testFile = ""
        was_benchmark_run = 0
        for file in os.listdir(suite.full_web_dir):
            if os.path.isfile(file) and file.endswith(".status"):
                index = string.rfind(file, ".status")
                tests.append(file[:index])

                with open(suite.full_web_dir + file, "r") as f:
                    for line in f:
                        if line.find("benchmarks updated") > 0:
                            was_benchmark_run = 1

            if os.path.isfile(file) and file.endswith(".ini"):
                testFile = file


        # create the report for this test run
        num_failed = report_this_test_run(suite, was_benchmark_run,
                                          "recreated report after crash of suite",
                                          "",  
                                          tests, args.complete_report_from_crash, testFile)

        # create the suite report
        suite.log.bold("creating suite report...")
        report_all_runs(suite, activeTestList)
        suite.log.close_log()
        sys.exit("done")


    #--------------------------------------------------------------------------
    # check bench dir and create output directories
    #--------------------------------------------------------------------------
    all_compile = all([t.compileTest == 1 for t in testList])

    if not all_compile:
        bench_dir = suite.get_bench_dir()

    last_run = suite.get_last_run()        
    suite.make_test_dirs()

    if not args.copy_benchmarks is None:
        old_full_test_dir = suite.testTopDir + suite.suiteName + "-tests/" + last_run
        copy_benchmarks(old_full_test_dir, suite.full_web_dir, testList, bench_dir, suite.log)

        num_failed = report_this_test_run(suite, args.copy_benchmarks,   # plays the role of make_benchmarks here
                                          "copy_benchmarks used -- no new tests run",
                                          "",  
                                          testList, args.input_file[0])
        report_all_runs(suite, activeTestList)        
        sys.exit("done")

    
    #--------------------------------------------------------------------------
    # figure out what needs updating and do the git updates, save the
    # current hash / HEAD, and make a ChangeLog
    # --------------------------------------------------------------------------
    now = time.localtime(time.time())
    updateTime = time.strftime("%Y-%m-%d %H:%M:%S %Z", now)

    no_update = args.no_update.lower()
    if not args.copy_benchmarks is None:
        no_update = "all"

    # the default is to update everything, unless we specified a hash
    # when constructing the Repo object
    if no_update == "none":
        pass

    elif no_update == "all":
        for k in suite.repos:
            suite.repos[k].update = False

    else:
        nouplist = no_update.split(",")

        if "boxlib" in nouplist: suite.repos["BoxLib"].update = False
        if suite.srcName.lower() in nouplist: suite.repos["source"].update = False
        if suite.extSrcName.lower() in nouplist: suite.repos["extra_source"].update = False

        # each extra build directory has its own update flag
        for n, e in enumerate(suite.extra_build_names):
            if e.lower() in nouplist:
                suite.repos["extra_build-{}".format(n)].update = False

    os.chdir(suite.testTopDir)

    for k in suite.repos:
        suite.log.skip()
        suite.log.bold("repo: {}".format(suite.repos[k].name))
        suite.log.indent()

        if suite.repos[k].update or suite.repos[k].hash_wanted:
            suite.repos[k].git_update()

        suite.repos[k].save_head()

        if suite.repos[k].update:
            suite.repos[k].make_changelog()

        suite.log.outdent()

    #--------------------------------------------------------------------------
    # build the tools and do a make clean, only once per build directory
    #--------------------------------------------------------------------------
    suite.build_tools(testList)

    all_build_dirs = find_build_dirs(testList)

    suite.log.skip()
    suite.log.bold("make clean in...")

    for dir, source_tree in all_build_dirs:

        if source_tree > 0:
            suite.log.log("{} in {}".format(dir, suite.extra_build_names[source_tree-1]))
            os.chdir(suite.extra_build_dirs[source_tree-1] + dir)
        else:
            suite.log.log("{}".format(dir))
            os.chdir(suite.sourceDir + dir)

        suite.make_realclean()

    os.chdir(suite.testTopDir)


    #--------------------------------------------------------------------------
    # main loop over tests
    #--------------------------------------------------------------------------
    for test in testList:

        suite.log.outdent()  # just to make sure we have no indentation
        suite.log.skip()
        suite.log.bold("working on test: {}".format(test.name))
        suite.log.indent()

        if not args.make_benchmarks == None and (test.restartTest or test.compileTest or
                                                 test.selfTest):
            suite.log.warn("  WARNING: test {} doesn't need benchmarks... skipping".format(test.name))
            continue

        outputDir = suite.full_test_dir + test.name + '/'
        os.mkdir(outputDir)


        #----------------------------------------------------------------------
        # compile the code
        #----------------------------------------------------------------------
        if test.useExtraBuildDir > 0:
            bDir = suite.extra_build_dirs[test.useExtraBuildDir-1] + test.buildDir
        else:
            bDir = suite.sourceDir + test.buildDir

        os.chdir(bDir)

        if test.reClean == 1:
            # for one reason or another, multiple tests use different
            # build options, make clean again to be safe
            suite.log.log("re-making clean...")
            suite.make_realclean()

        suite.log.log("building...")

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

            comp_string = "%s -j%s BOXLIB_HOME=%s %s %s DIM=%d %s COMP=%s FCOMP=%s executable=%s" % \
                (suite.MAKE, suite.numMakeJobs, suite.boxLibDir,
                 suite.extSrcCompString, test.addToCompileString,
                 test.dim, buildOptions, suite.COMP, suite.FCOMP,
                 executable)

            suite.log.log(comp_string)
            so, se, r = run(comp_string,
                            outfile="{}/{}.make.out".format(outputDir, test.name))

        elif suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src":

            build_options = ""
            build_options += "NDEBUG={} ".format(convert_to_f_make_flag(test.debug, test_not=True))
            build_options += "MPI={} ".format(convert_to_f_make_flag(test.useMPI))
            build_options += "OMP={} ".format(convert_to_f_make_flag(test.useOMP))

            if test.useExtraBuildDir > 0:
                build_options += suite.extraBuildDirCompString + " "

            comp_string = suite.build_f(opts="{} {} {}".format(
                suite.extSrcCompString, test.addToCompileString, build_options),
                          outfile="{}/{}.make.out".format(outputDir, test.name))

            # we need a better way to get the executable name here
            executable = get_recent_filename(bDir,"main",".exe")

        test.comp_string = comp_string

        if test.compileTest:

            # compilation tests are done now -- just make the report and ...
            shutil.copy("%s/%s.make.out"    % (outputDir, test.name), suite.full_web_dir)

            suite.log.log("creating problem test report ...")
            report_single_test(suite, test)

            # ... skip to the next test in the loop
            continue


        #----------------------------------------------------------------------
        # copy the necessary files over to the run directory
        #----------------------------------------------------------------------
        suite.log.log("copying files to run directory...")

        try: shutil.copy(executable, outputDir)
        except (IOError, AttributeError):

            # compilation failed.  First copy the make.out into the
            # web directory and then report
            shutil.copy("%s/%s.make.out" % (outputDir, test.name), suite.full_web_dir)

            errorMsg = "    ERROR: compilation failed"
            report_test_failure(suite, errorMsg, test)
            continue

        try: shutil.copy(test.inputFile, outputDir)
        except IOError:
            errorMsg = "    ERROR: unable to copy input file: %s" % test.inputFile
            report_test_failure(suite, errorMsg, test)
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
                report_test_failure(suite, errorMsg, test)
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
                report_test_failure(suite, errorMsg, test)
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
                report_test_failure(suite, errorMsg, test)
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
                    report_test_failure(suite, errorMsg, test)
                    skip_to_next_test = 1
                    break

        if skip_to_next_test: continue


        #----------------------------------------------------------------------
        # run the test
        #----------------------------------------------------------------------
        suite.log.log("running the test...")

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

        elif suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src":

            base_command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk " % \
                           (executable, test.inputFile, test.name, test.name)

            # keep around the checkpoint files only for the restart runs
            if not test.restartTest: base_command += " --chk_int 0 "

            base_command += "{}".format(suite.globalAddToExecString)

        suite.run_test(test, base_command)


        # if it is a restart test, then rename the final output file and
        # restart the test
        if test.restartTest:
            lastFile = test.get_last_plotfile(output_dir=outputDir)

            if lastFile == "":
                errorMsg = "ERROR: test did not produce output.  Restart test not possible"
                report_test_failure(suite, errorMsg, test)
                continue

            origLastFile = "orig_%s" % (lastFile)
            shutil.move(lastFile, origLastFile)

            if test.diffDir:
                origDiffDir = "orig_%s" % (test.diffDir)
                shutil.move(test.diffDir, origDiffDir)

            # get the file number to restart from
            restartFile = "%s_chk%5.5d" % (test.name, test.restartFileNum)

            suite.log.log("restarting from {} ... ".format(restartFile))

            if suite.sourceTree == "C_Src" or test.testSrcTree == "C_Src":

                base_command = "./%s %s amr.plot_file=%s_plt amr.check_file=%s_chk amr.checkpoint_files_output=0 amr.restart=%s" % \
                        (executable, test.inputFile, test.name, test.name, restartFile)

            elif suite.sourceTree == "F_Src" or test.testSrcTree == "F_Src":

                base_command = "./%s %s --plot_base_name %s_plt --check_base_name %s_chk --chk_int 0 --restart %d %s" % \
                        (executable, test.inputFile, test.name, test.name, test.restartFileNum, suite.globalAddToExecString)

            suite.run_test(test, base_command)

        test.wallTime = time.time() - test.wallTime


        #----------------------------------------------------------------------
        # do the comparison
        #----------------------------------------------------------------------
        if not test.selfTest:

            if test.outputFile == "":
                if test.compareFile == "":
                    compareFile = test.get_last_plotfile(output_dir=outputDir)
                else:
                    compareFile = test.compareFile
                outputFile = compareFile
            else:
                outputFile = test.outputFile
                compareFile = test.name+'_'+outputFile


            # get the number of levels for reporting
            prog = "{} -l {}".format(suite.tools["fboxinfo"], outputFile)
            stdout0, stderr0, rc = run(prog)
            test.nlevels = stdout0.rstrip('\n')
            if not type(convert_type(test.nlevels)) is int:
                test.nlevels = ""

            if args.make_benchmarks == None:

                suite.log.log("doing the comparison...")
                suite.log.indent()
                suite.log.log("comparison file: {}".format(outputFile))

                test.compare_file_used = outputFile

                if not test.restartTest:
                    benchFile = bench_dir + compareFile
                else:
                    benchFile = origLastFile

                # see if it exists
                # note, with BoxLib, the plotfiles are actually directories

                if not os.path.isdir(benchFile):
                    suite.log.warn("WARNING: no corresponding benchmark found")
                    benchFile = ""

                    cf = open("%s.compare.out" % (test.name), 'w')
                    cf.write("WARNING: no corresponding benchmark found\n")
                    cf.write("         unable to do a comparison\n")
                    cf.close()

                else:
                    if not compareFile == "":

                        suite.log.log("benchmark file: {}".format(benchFile))

                        command = "{} -n 0 --infile1 {} --infile2 {}".format(
                            suite.tools["fcompare"], benchFile, outputFile)
                        sout, serr, ierr = run(command, outfile="{}.compare.out".format(test.name), store_command=True)

                    else:
                        suite.log.warn("WARNING: unable to do a comparison")

                        cf = open("%s.compare.out" % (test.name), 'w')
                        cf.write("WARNING: run did not produce any output\n")
                        cf.write("         unable to do a comparison\n")
                        cf.close()

                suite.log.outdent()

                if not test.diffDir == "":
                    if not test.restartTest:
                        diffDirBench = bench_dir + '/' + test.name + '_' + test.diffDir
                    else:
                        diffDirBench = origDiffDir

                    suite.log.log("doing the diff...")
                    suite.log.log("diff dir: {}".format(test.diffDir))

                    command = "diff %s -r %s %s" \
                        % (test.diffOpts, diffDirBench, test.diffDir)

                    outfile = "{}.compare.out".format(test.name)
                    sout, serr, diff_status = run(command, outfile=outfile, store_command=True)

                    if diff_status == 0:
                        with open("{}.compare.out".format(test.name), 'a') as cf:
                            cf.write("\ndiff was SUCCESSFUL\n")

            else:   # make_benchmarks

                suite.log.log("storing output of {} as the new benchmark...".format(test.name))
                suite.log.indent()
                suite.log.log("new benchmark file: {}".format(compareFile))
                suite.log.outdent()

                if not compareFile == "":
                    if not outputFile == compareFile:
                        source_file = outputFile
                    else:
                        source_file = compareFile

                    try: shutil.rmtree("{}/{}".format(bench_dir, compareFile))
                    except: pass
                    shutil.copytree(source_file, "{}/{}".format(bench_dir, compareFile))

                    with open("%s.status" % (test.name), 'w') as cf:
                        cf.write("benchmarks updated.  New file:  %s\n" % (compareFile) )

                else:
                    with open("%s.status" % (test.name), 'w') as cf:
                        cf.write("benchmarks failed")

                if not test.diffDir == "":
                    diffDirBench = "{}/{}_{}".format(bench_dir, test.name, test.diffDir)
                    if os.path.isdir(diffDirBench):
                        shutil.rmtree(diffDirBench)
                        shutil.copytree(test.diffDir, diffDirBench)
                    else:
                        shutil.copy(test.diffDir, diffDirBench)
                    suite.log.log("new diffDir: {}_{}".format(test.name, test.diffDir))

        else:   # selfTest

            if args.make_benchmarks == None:

                suite.log.log("looking for selfTest success string: {} ...".format(test.stSuccessString))

                try: of = open("%s.run.out" % (test.name), 'r')
                except IOError:
                    suite.log.warn("WARNING: no output file found")
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

                with open("%s.compare.out" % (test.name), 'w') as cf:
                    if compareSuccessful:
                        cf.write("SELF TEST SUCCESSFUL\n")
                    else:
                        cf.write("SELF TEST FAILED\n")


        #----------------------------------------------------------------------
        # do any requested visualization (2- and 3-d only) and analysis
        #----------------------------------------------------------------------
        if outputFile != "":
            if args.make_benchmarks == None:

                if test.doVis: 

                    if test.dim == 1:
                        suite.log.log("Visualization not supported for dim = {}".format(test.dim))
                    else:
                        suite.log.log("doing the visualization...")
                        tool = suite.tools["fsnapshot{}d".format(test.dim)]
                        run('{} --palette {}/Palette -cname "{}" -p "{}"'.format(
                            tool, suite.compare_tool_dir, test.visVar, outputFile))

                        # convert the .ppm files into .png files
                        ppm_file = get_recent_filename(outputDir, "", ".ppm")
                        png_file = ppm_file.replace(".ppm", ".png")
                        run("convert {} {}".format(ppm_file, png_file))

                if not test.analysisRoutine == "":

                    suite.log.log("doing the analysis...")
                    if test.useExtraBuildDir > 0:
                        tool = "{}/{}".format(suite.extra_build_dirs[test.useExtraBuildDir-1], test.analysisRoutine)
                    else:
                        tool = "{}/{}".format(suite.sourceDir, test.analysisRoutine)

                    shutil.copy(tool, os.getcwd())

                    option = eval("suite.{}".format(test.analysisMainArgs))
                    run("{} {} {}".format(os.path.basename(test.analysisRoutine),
                                          option, outputFile))

        else:
            if test.doVis or test.analysisRoutine != "":
                suite.log.warn("WARNING: no output file.  Skipping visualization")


        #----------------------------------------------------------------------
        # move the output files into the web directory
        #----------------------------------------------------------------------
        if args.make_benchmarks == None:
            shutil.copy("%s.run.out"     % (test.name), suite.full_web_dir)
            shutil.copy("%s.make.out"    % (test.name), suite.full_web_dir)
            shutil.copy("%s.compare.out" % (test.name), suite.full_web_dir)

            shutil.copy(test.inputFile, "%s/%s.%s" % (suite.full_web_dir, test.name, test.inputFile) )

            if suite.sourceTree == "C_Src":
                shutil.copy(test.probinFile, "%s/%s.%s" % (suite.full_web_dir, test.name, test.probinFile) )

            for file in test.auxFiles:

                # sometimes the auxFile was in a subdirectory under the
                # build directory.
                index = string.rfind(file, "/")
                if index > 0:
                    file = file[index+1:]

                shutil.copy(file, "%s/%s.%s" % (suite.full_web_dir, test.name, file) )

            if test.doVis:
               png_file = get_recent_filename(outputDir, "", ".png")
               if not png_file is None:
                   try: shutil.copy(png_file, suite.full_web_dir)
                   except IOError:
                       # visualization was not successful.  Reset doVis
                       test.doVis = 0

            if not test.analysisRoutine == "":
                try: shutil.copy(test.analysisOutputImage, suite.full_web_dir)
                except IOError:
                    # analysis was not successful.  Reset the output image
                    test.analysisOutputImage = ""


        else:
            shutil.copy("%s.status" % (test.name), suite.full_web_dir)


        #----------------------------------------------------------------------
        # archive (or delete) the output
        #----------------------------------------------------------------------
        suite.log.log("archiving the output...")
        for file in os.listdir(outputDir):
            if (os.path.isdir(file) and
                (file.startswith("%s_plt" % (test.name)) or
                 file.startswith("%s_chk" % (test.name)) ) ):

                if suite.purge_output == 1 and not file == outputFile:
                    # delete the plt/chk file
                    if os.path.isdir(file):
                        try: shutil.rmtree(file)
                        except:
                            suite.log.warn("WARNING: unable to remove {}".format(file))

                else:
                    # tar it up
                    try:
                        tar = tarfile.open("%s.tgz" % (file), "w:gz")
                        tar.add("%s" % (file))
                        tar.close()

                    except:
                        suite.log.warn("WARNING: unable to tar output file %s" % (file))

                    else:
                        shutil.rmtree(file)


        #----------------------------------------------------------------------
        # write the report for this test
        #----------------------------------------------------------------------
        if args.make_benchmarks == None:
            suite.log.log("creating problem test report ...")
            report_single_test(suite, test)
        

    #--------------------------------------------------------------------------
    # write the report for this instance of the test suite
    #--------------------------------------------------------------------------
    suite.log.outdent()
    suite.log.skip()
    suite.log.bold("creating new test report...")
    num_failed = report_this_test_run(suite, args.make_benchmarks, args.note,
                                      updateTime,  
                                      testList, args.input_file[0])


    # make sure that all of the files in the web directory are world readable
    for file in os.listdir(suite.full_web_dir):
       currentFile = suite.full_web_dir + file

       if os.path.isfile(currentFile):
          os.chmod(currentFile, 0o644)

    # reset the branch to what it was originally
    suite.log.skip()
    suite.log.bold("reverting git branches/hashes")
    suite.log.indent()

    for k in suite.repos:
        if suite.repos[k].update or suite.repos[k].hash_wanted:
            suite.repos[k].git_back()

    suite.log.outdent()

    # For temporary run, return now without creating suote report.
    if args.do_temp_run:
        return num_failed


    # store an output file in the web directory that can be parsed easily by
    # external program
    name = "source"
    if suite.sourceTree == "BoxLib": name = "BoxLib"
    branch = suite.repos[name].branch_wanted.strip("\"")
    
    with open("{}/suite.{}.status".format(suite.webTopDir, branch), "w") as f:
        f.write("{}; num failed: {}; source hash: {}".format(
            suite.repos[name].name, num_failed, suite.repos[name].hash_current))


    #--------------------------------------------------------------------------
    # generate the master report for all test instances
    #--------------------------------------------------------------------------
    suite.log.skip()
    suite.log.bold("creating suite report...")
    report_all_runs(suite, activeTestList)

    def emailDevelopers():
        msg = email.message_from_string(suite.emailBody)
        msg['From'] = suite.emailFrom
        msg['To'] = ",".join(suite.emailTo)
        msg['Subject'] = suite.emailSubject

        server = smtplib.SMTP('localhost')
        server.sendmail(suite.emailFrom, suite.emailTo, msg.as_string())
        server.quit()

    if num_failed > 0 and suite.sendEmailWhenFail and not args.send_no_email:
        suite.log.skip()
        suite.log.bold("sending email...")
        emailDevelopers()


    return num_failed


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
           color: black; background-color: lime; padding: 2px;}

a.passed:link {color: black; text-decoration: none;}
a.passed:visited {color: black; text-decoration: none;}
a.passed:hover {color: #ee00ee; text-decoration: underline;}

h3.failed {text-decoration: none; display: inline;
           color: yellow; background-color: red; padding: 2px;}

a.failed:link {color: yellow; text-decoration: none;}
a.failed:visited {color: yellow; text-decoration: none;}
a.failed:hover {color: #00ffff; text-decoration: underline;}

h3.benchmade {text-decoration: none; display: inline;
              color: black; background-color: orange; padding: 2px;}

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

#summary td.passed {background-color: lime;}
#summary td.failed {background-color: red;}
#summary td.benchmade {background-color: orange;}

th {background-color: grey;
    color: yellow;
    text-align: center;
    vertical-align: bottom;
    height: @TABLEHEIGHT@;
    padding-bottom: 3px;
    padding-left: 5px;
    padding-right: 5px;}

li {padding-top: 0.5em;}

ul li {color: blue;
       font-weight: bold;}
ul li ul li {color: black;
             font-weight: normal;}

ul li h3 {border: 1px solid black;}

#compare td {font-family: "Lucida Console", Monaco, monospace;
             font-size: 80%;}

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

    with open("tests.css", 'w') as cf:
        cf.write(cssC)

class HTMLList(object):
    def __init__(self, of=None):
        self.list_items = []
        self.current_indent = 0
        self.of = of

    def item(self, content):
        self.list_items.append((self.current_indent, content))

    def indent(self):
        self.current_indent += 1

    def outdent(self):
        self.current_indent -= 1

    def write_list(self):
        self.of.write("<ul>\n")
        current_indent = -1
        for i, c in self.list_items:
            if current_indent == -1:
                current_indent = i
            else:
                if i < current_indent:
                    self.of.write("</li></ul></li>\n")
                elif i > current_indent:
                    self.of.write("<ul>\n")
                else:
                    self.of.write("</li>\n")

            current_indent = i
            self.of.write("<li>{}\n".format(c))

        # finish current item
        self.of.write("</li>")

        # finish nesting
        for n in range(0,current_indent):
            self.of.write("</ul></li>\n")

        self.of.write("</ul>\n")

class HTMLTable(object):
    def __init__(self, out_file, columns=1, divs=None):
        self.hf = out_file
        self.columns = columns
        if not divs is None:
            self.divs = list(divs)
        else:
            self.divs = None

    def start_table(self):
        if not self.divs is None:
            for d in self.divs:
                self.hf.write("<div id='{}'>\n".format(d))
        self.hf.write("<p><table>\n")

    def header(self, header_list):
        n = len(header_list)
        line = "<tr>"+n*"<th>{}</th>"+"</tr>\n"
        self.hf.write(line.format(*header_list))

    def print_single_row(self, row):
        self.hf.write("<tr class='special'><td colspan={}>".format(self.columns)+row+"</td></tr>\n")

    def print_row(self, row_list, highlight=False):
        """ row_list are the individual table elements.  Note that if
        a list item is a tuple, then the first element is assumed to
        be the cell data and the second element is an html tag that
        goes in the <td >, e.g. to set the class or colspan"""

        n = len(row_list)
        if highlight:
            line = "<tr>"+n*"<td class='highlight'>{}</td>"+"</tr>\n"
        else:
            line = "<tr>"
            for d in row_list:
                if isinstance(d, tuple):
                    line += "<td {}>{}</td>".format(d[1], d[0])
                else:
                    line += "<td>{}</td>".format(d)
            line += "</tr>\n"
        self.hf.write(line.format(*row_list))

    def end_table(self):
        self.hf.write("</table>\n")
        if not self.divs is None:
            for n in range(len(self.divs)):
                self.hf.write("</div>\n")


#==============================================================================
# REPORT ROUTINES
#==============================================================================
def report_single_test(suite, test):
    """ generate a single problem's test result page """

    # get the current directory
    currentDir = os.getcwd()

    # switch to the web directory and open the report file
    os.chdir(suite.full_web_dir)

    #--------------------------------------------------------------------------
    # parse the compilation report and determine if we compiled
    #--------------------------------------------------------------------------
    compileFile = "%s.make.out" % (test.name)

    try: cf = open(compileFile, 'r')
    except IOError:
        suite.log.warn("WARNING: no compilation file found")
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
    if not test.compileTest:
        compareFile = "%s.compare.out" % (test.name)

        try: cf = open(compareFile, 'r')
        except IOError:
            suite.log.warn("WARNING: no comparison file found")
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
    with open(status_file, 'w') as sf:
        if (compileSuccessful and
            (test.compileTest or (not test.compileTest and compareSuccessful))):
            sf.write("PASSED\n")
            suite.log.success("{} PASSED".format(test.name))
        else:
            sf.write("FAILED\n")
            suite.log.testfail("{} FAILED".format(test.name))


    #--------------------------------------------------------------------------
    # generate the HTML page for this test
    #--------------------------------------------------------------------------

    # write the css file
    create_css()

    htmlFile = "%s.html" % (test.name)
    hf = open(htmlFile, 'w')

    newHead = HTMLHeader + r"""<CENTER><H1><A HREF="index.html">@TESTDIR@</A> / @TESTNAME@</H1></CENTER>"""

    newHead = newHead.replace("@TESTDIR@", os.path.normpath(suite.test_dir))
    newHead = newHead.replace("@TESTNAME@", test.name)

    hf.write(newHead)

    ll = HTMLList(of=hf)

    # build summary
    ll.item("Build/Test information:")
    ll.indent()

    ll.item("Build directory: {}".format(test.buildDir))

    if test.useExtraBuildDir > 0:
        ll.indent()
        ll.item("in {}".format(suite.extra_build_dirs[test.useExtraBuildDir-1]))
        ll.outdent()

    if not test.compileTest:

        if test.debug:
            ll.item("Debug test")

        if test.useMPI or test.useOMP:
            ll.item("Parallel run")
            ll.indent()
            if test.useMPI:
                ll.item("MPI numprocs = {}".format(test.numprocs))
            if test.useOMP:
                ll.item("OpenMP numthreads = {}".format(test.numthreads))
            ll.outdent()

        if test.restartTest:

            ll.item("Restart test")
            ll.indent()
            ll.item("Job was run as normal and then restarted from checkpoint # {}, and the two final outputs were compared".format(test.restartFileNum))
            ll.outdent()


        ll.item("Files:")
        ll.indent()

        ll.item("input file: <a href=\"{}.{}\">{}</a>".format(test.name, test.inputFile, test.inputFile))

        if suite.sourceTree == "C_Src":
            ll.item("probin file: <a href=\"{}.{}\">{}</a>".format(test.name, test.probinFile, test.probinFile))

        for i, afile in enumerate(test.auxFiles):
            # sometimes the auxFile was in a subdirectory under the
            # build directory.
            index = string.rfind(afile, "/")
            if index > 0:
                root_file = afile[index+1:]
            else:
                root_file = afile

            ll.item("auxillary file {}: <a href=\"{}.{}\">{}</a>".format(i+1, test.name, root_file, afile))

        ll.outdent()

        ll.item("Dimensionality: {}".format(test.dim))

    ll.outdent()   # end of build information

    # compilation summary
    ll.item("Compilation:")
    ll.indent()

    if compileSuccessful:
        ll.item("<h3 class=\"passed\">Successful</h3>")
    else:
        ll.item("<h3 class=\"failed\">Failed</h3>")

    ll.item("Compilation command:<br>{}".format(test.comp_string))
    ll.item("<a href=\"{}.make.out\">make output</a>".format(test.name))

    ll.outdent()


    if not test.compileTest:

        # execution summary
        ll.item("Execution:")
        ll.indent()
        ll.item("Execution time: {:.3f} s".format(test.wallTime))
        ll.item("Execution command:<br>{}".format(test.run_command))
        ll.item("<a href=\"{}.run.out\">execution output</a>".format(test.name))
        ll.outdent()

        # comparison summary
        ll.item("Comparison: ")
        ll.indent()

        if compareSuccessful:
            ll.item("<h3 class=\"passed\">Successful</h3>")
        else:
            ll.item("<h3 class=\"failed\">Failed</h3>")

    ll.write_list()


    if not test.compileTest:

        # parse the compare output and make an HTML table
        ht = HTMLTable(hf, columns=3, divs=["summary", "compare"])
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

                    hf.write(line)
                    in_diff_region = True
                    continue

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
                    if "NaN present" in line:
                        ht.print_row([fields[0], (fields[1], "colspan='2'")])
                        continue
                    else:
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
                hf.write(line)

        if in_diff_region:
            hf.write("</pre>\n")
        else:
            ht.end_table()

        if box_error:
            hf.write("<p>number of boxes do not match</p>\n")

        # show any visualizations
        if test.doVis:
            png_file = get_recent_filename(suite.full_web_dir, test.name, ".png")
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


def report_test_failure(suite, message, test):
    """ generate a simple report for an error encountered while performing
        the test """

    suite.log.testfail("aborting test")
    suite.log.testfail(message)

    # get the current directory
    currentDir = os.getcwd()

    # switch to the web directory and open the report file
    os.chdir(suite.full_web_dir)

    #--------------------------------------------------------------------------
    # write out the status file for this problem -- FAILED
    #--------------------------------------------------------------------------
    status_file = "%s.status" % (test.name)
    with open(status_file, 'w') as sf:
        sf.write("FAILED\n")

    suite.log.testfail("{} FAILED".format(test.name))


    #--------------------------------------------------------------------------
    # generate the HTML page for this test
    #--------------------------------------------------------------------------

    # check to see if the CSS file is present, if not, write it
    if (not os.path.isfile("tests.css")):
        create_css()


    htmlFile = "%s.html" % (test.name)
    hf = open(htmlFile, 'w')

    newHead = HTMLHeader + r"""<CENTER><H1><A HREF="index.html">@TESTDIR@</A> / @TESTNAME@</H1></CENTER>"""

    newHead = newHead.replace("@TESTDIR@", suite.test_dir)
    newHead = newHead.replace("@TESTNAME@", test.name)

    hf.write(newHead)

    # write out the information about the test
    hf.write("<P><b>build directory:</b> %s\n" % (test.buildDir) )

    hf.write("<P><H3 CLASS=\"failed\">Test Failed</H3></P>\n")
    hf.write("<P>%s</P>\n" % (message) )

    if (not test.comp_string == None):
        hf.write("<P>compliation command:\n %s\n" % (test.comp_string) )
        hf.write("<P><A HREF=\"%s.make.out\">make output</A>\n" % (test.name) )


    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")

    hf.close()


    # switch back to the original directory
    os.chdir(currentDir)


def report_this_test_run(suite, make_benchmarks, note, update_time,
                      testList, testFile):
    """ generate the master page for a single run of the test suite """

    # get the current directory
    currentDir = os.getcwd()

    # switch to the web directory and open the report file
    os.chdir(suite.full_web_dir)


    # keep track of the number of tests that passed and the number that failed
    num_failed = 0
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
    newHead = newHead.replace("@TESTNAME@", suite.test_dir)

    hf.write(newHead)

    if (not note == ""):
       hf.write("<p><b>Test run note:</b><br><font color=\"gray\">%s</font>\n" % (note) )

    if not make_benchmarks == None:
       hf.write("<p><b>Benchmarks updated</b><br>comment: <font color=\"gray\">{}</font>\n".format(make_benchmarks) )
       hf.write("<p>&nbsp;\n")


    hf.write("<p><b>test input parameter file:</b> <A HREF=\"%s\">%s</A>\n" %
             (testFile, testFile) )

    any_update = any([suite.repos[t].update for t in suite.repos])

    if any_update and not update_time == "":
        hf.write("<p>&nbsp;\n")
        hf.write("<p><b>Git update was done at: </b>%s\n" % (update_time) )

        hf.write("<ul>\n")
        code_str = "<li><b>{}</b><ul><li><b>branch:</b> {}; <b>hash:</b> {}</li><li><b>changelog:</b> <a href=\"{}\">{}</a></li></ul></li>"

        for k, r in suite.repos.items():
            if r.update:
                hf.write(code_str.format(r.name, r.branch_wanted, r.hash_current,
                                         "ChangeLog.{}".format(r.name),
                                         "ChangeLog.{}".format(r.name)))

        hf.write("</ul>")

    else:
        hf.write("<p>No git update done\n")

    hf.write("<p>&nbsp;\n")

    if make_benchmarks == None:
        ht = HTMLTable(hf, columns=11, divs=["summary"])
        ht.start_table()
        ht.header(["test name", "dim", "compare plotfile",
                   "# levels", "MPI (# procs)", "OMP (# threads)", "debug?",
                   "compile?", "restart?", "wall time", "result"])

    else:
        ht = HTMLTable(hf, columns=3, divs=["summary"])
        ht.start_table()
        ht.header(["test name", "result", "comment"])

    # loop over the tests and add a line for each
    for test in testList:

        if make_benchmarks == None:

            # check if it passed or failed
            status_file = "%s.status" % (test.name)

            testPassed = 0

            with open(status_file, 'r') as sf:
                for line in sf:
                    if line.find("PASSED") >= 0:
                        testPassed = 1
                        numPassed += 1
                        break

                if not testPassed:
                    num_failed += 1

            row_info = []
            row_info.append("<a href=\"{}.html\">{}</a>".format(test.name, test.name))
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
                row_info.append(("PASSED", "class='passed'"))
            else:
                row_info.append(("FAILED", "class='failed'"))

            ht.print_row(row_info)

        else:
            if test.restartTest: continue
            if test.compileTest: continue
            if test.selfTest: continue

            # the benchmark was updated -- find the name of the new benchmark file
            benchStatusFile = "%s.status" % (test.name)

            benchFile = "none"

            with open(benchStatusFile, 'r') as bf:
                for line in bf:
                    index = line.find("file:")
                    if index >= 0:
                        benchFile = line[index+5:]
                        break

            row_info = []
            row_info.append("{}".format(test.name))
            if not benchFile == "none":
                row_info.append(("BENCHMARK UPDATED", "class='benchmade'"))
                row_info.append("new benchmark file is {}".format(benchFile))
            else:
                row_info.append(("BENCHMARK NOT UPDATED", "class='failed'"))
                row_info.append("compilation or execution failed")

            ht.print_row(row_info)

    ht.end_table()

    # close
    hf.write("</BODY>\n")
    hf.write("</HTML>\n")
    hf.close()


    #--------------------------------------------------------------------------
    # write out a status file for all the tests
    #--------------------------------------------------------------------------

    status_file = os.path.normpath(suite.test_dir) + ".status"
    with open(status_file, 'w') as sf:

        if make_benchmarks == None:
            if num_failed == 0:
                sf.write("ALL PASSED\n")
            elif num_failed > 0 and numPassed > 0:
                sf.write("SOME FAILED\n")
            else:
                sf.write("ALL FAILED\n")

        else:
            sf.write("BENCHMARKS UPDATED\n")

    # switch back to the original directory
    os.chdir(currentDir)

    return num_failed


def report_all_runs(suite, activeTestList):

    tableHeight = min(max(suite.lenTestName, 4), 16)
    
    os.chdir(suite.webTopDir)

    create_css(tableHeight=tableHeight)

    validDirs = []
    allTests = []


    #--------------------------------------------------------------------------
    # start by finding the list of valid test directories
    #--------------------------------------------------------------------------
    for file in os.listdir(suite.webTopDir):

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

        for file in os.listdir(suite.webTopDir + dir):

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
            status_file = "%s/%s/%s.status" % (suite.webTopDir, dir, test)
            if os.path.isfile(status_file):
                valid = 1
                break

        if not valid: continue

        # write out the directory (date)
        hf.write("<TR><TD class='date'><SPAN CLASS='nobreak'><A class='main' HREF=\"%s/index.html\">%s&nbsp;</A></SPAN></TD>\n" %
                 (dir, dir) )

        for test in allTests:

            # look to see if the current test was part of this suite run
            status_file = "%s/%s/%s.status" % (suite.webTopDir, dir, test)
            status = 0

            if os.path.isfile(status_file):

                with open(status_file, 'r') as sf:

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
    test_suite(sys.argv[1:])
