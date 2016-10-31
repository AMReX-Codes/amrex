try: import ConfigParser as configparser
except ImportError:
    import configparser   # python 3

import getpass
import os
import socket

import repo
import suite
import test_util

def convert_type(string):
    """ return an integer, float, or string from the input string """
    if string is None:
        return None

    try: int(string)
    except: pass
    else: return int(string)

    try: float(string)
    except: pass
    else: return float(string)

    return string.strip()

def safe_get(cp, sec, opt, default=None):
    try: v = cp.get(sec, opt)
    except: v = default
    return v

def load_params(args):
    """
    reads the parameter file and creates as list of test objects as well as
    the suite object
    """

    test_list = []

    cp = configparser.ConfigParser()    # note, need strict=False for Python3
    cp.optionxform = str

    log = test_util.Log()

    log.bold("loading " + args.input_file[0])

    try: cp.read(args.input_file[0])
    except:
        log.fail("ERROR: unable to read parameter file {}".format(file))

    # "main" is a special section containing the global suite parameters.
    mysuite = suite.Suite(args)

    mysuite.log = log

    valid_options = list(mysuite.__dict__.keys())

    for opt in cp.options("main"):

        # get the value of the current option
        value = convert_type(cp.get("main", opt))

        if opt in valid_options:

            if opt == "sourceTree":
                if not value in ["C_Src", "F_Src", "BoxLib"]:
                    mysuite.log.fail("ERROR: invalid sourceTree")
                else:
                    mysuite.sourceTree = value

            elif opt == "testTopDir": mysuite.testTopDir = mysuite.check_test_dir(value)
            elif opt == "webTopDir": mysuite.webTopDir = os.path.normpath(value) + "/"

            elif opt == "emailTo": mysuite.emailTo = value.split(",")

            else:
                # generic setting of the object attribute
                setattr(mysuite, opt, value)

        else:
            mysuite.log.warn("suite parameter {} not valid".format(opt))


    # BoxLib -- this will always be defined
    rdir = mysuite.check_test_dir(safe_get(cp, "BoxLib", "dir"))

    branch = convert_type(safe_get(cp, "BoxLib", "branch"))
    rhash = convert_type(safe_get(cp, "BoxLib", "hash"))

    mysuite.repos["BoxLib"] = repo.Repo(mysuite, rdir, "BoxLib",
                                        branch_wanted=branch, hash_wanted=rhash)


    # now all the other build and source directories
    other_srcs = [s for s in cp.sections() if s.startswith("extra-")]
    if not mysuite.sourceTree == "BoxLib": other_srcs.append("source")

    for s in other_srcs:
        if s.startswith("extra-"):
            k = s.split("-")[1]
        else:
            k = "source"

        rdir = mysuite.check_test_dir(safe_get(cp, s, "dir"))
        branch = convert_type(safe_get(cp, s, "branch"))
        rhash = convert_type(safe_get(cp, s, "hash"))

        build = convert_type(safe_get(cp, s, "build", default=0))
        if s == "source": build = 1

        comp_string = safe_get(cp, s, "comp_string")

        name = os.path.basename(os.path.normpath(rdir))

        mysuite.repos[k] = repo.Repo(mysuite, rdir, name,
                                     branch_wanted=branch, hash_wanted=rhash,
                                     build=build, comp_string=comp_string)


    # BoxLib-only tests don't have a sourceDir
    mysuite.boxlib_dir = mysuite.repos["BoxLib"].dir

    if mysuite.sourceTree == "BoxLib":
        mysuite.source_dir = mysuite.repos["BoxLib"].dir
    else:
        mysuite.source_dir = mysuite.repos["source"].dir


    # now flesh out the compile strings -- they may refer to either themselves
    # or the source dir
    for r in mysuite.repos.keys():
        s = mysuite.repos[r].comp_string
        if not s is None:
            mysuite.repos[r].comp_string = \
                s.replace("@self@", mysuite.repos[r].dir).replace("@source@", mysuite.repos["source"].dir)


    # the suite needs to know any ext_src_comp_string
    for r in mysuite.repos.keys():
        if not mysuite.repos[r].build == 1:
            if not mysuite.repos[r].comp_string is None:
                mysuite.extra_src_comp_string += "{} ".format(mysuite.repos[r].comp_string)

    # checks
    if mysuite.sendEmailWhenFail and not args.send_no_email:
        if mysuite.emailTo == [] or mysuite.emailBody == "":
            mysuite.log.fail("ERROR: when sendEmailWhenFail = 1, you must specify emailTo and emailBody\n")

        if mysuite.emailFrom == "":
            mysuite.emailFrom = '@'.join((getpass.getuser(), socket.getfqdn()))

        if mysuite.emailSubject == "":
            mysuite.emailSubject = mysuite.suiteName+" Regression Test Failed"

    if mysuite.slack_post:
        if not os.path.isfile(mysuite.slack_webhookfile):
            mysuite.log.warn("slack_webhookfile invalid")
            mysuite.slack_post = 0
        else:
            print(mysuite.slack_webhookfile)
            try: f = open(mysuite.slack_webhookfile)
            except:
                mysuite.log.warn("unable to open webhook file")
                mysuite.slack_post = 0
            else:
                mysuite.slack_webhook_url = str(f.readline())
                f.close()

    if (mysuite.sourceTree == "" or mysuite.boxlib_dir == "" or
        mysuite.source_dir == "" or mysuite.testTopDir == ""):
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

        if sec in ["main", "BoxLib", "source"] or sec.startswith("extra-"): continue

        # maximum test name length -- used for HTML formatting
        mysuite.lenTestName = max(mysuite.lenTestName, len(sec))

        # create the test object for this test
        mytest = suite.Test(sec)
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
                mysuite.log.warn("unrecognized parameter {} for test {}".format(opt, sec))


        # make sure that the build directory actually exists
        if not mytest.extra_build_dir == "":
            bdir = mysuite.repos[mytest.extra_build_dir].dir + mytest.buildDir
        else:
            bdir = mysuite.source_dir + mytest.buildDir

        if not os.path.isdir(bdir):
            mysuite.log.warn("invalid build directory: {}".format(bdir))
            invalid = 1


        # make sure all the require parameters are present
        if mytest.compileTest:
            if mytest.buildDir == "":
                mysuite.log.warn("mandatory parameters for test {} not set".format(sec))
                invalid = 1

        else:
            if mytest.buildDir == "" or mytest.inputFile == "" or mytest.dim == -1:
                warn_msg = ["required params for test {} not set".format(sec),
                            "buildDir = {}".format(mytest.buildDir),
                            "inputFile = {}".format(mytest.inputFile)]
                warn_msg += ["dim = {}".format(mytest.dim)]
                mysuite.log.warn(warn_msg)

                invalid = 1

        # check the optional parameters
        if mytest.restartTest and mytest.restartFileNum == -1:
            mysuite.log.warn("restart-test {} needs a restartFileNum".format(sec))
            invalid = 1

        if mytest.selfTest and mytest.stSuccessString == "":
            mysuite.log.warn("self-test {} needs a stSuccessString".format(sec))
            invalid = 1

        if mytest.useMPI and mytest.numprocs == -1:
            mysuite.log.warn("MPI parallel test {} needs numprocs".format(sec))
            invalid = 1

        if mytest.useOMP and mytest.numthreads == -1:
            mysuite.log.warn("OpenMP parallel test {} needs numthreads".format(sec))
            invalid = 1

        if mytest.doVis and mytest.visVar == "":
            mysuite.log.warn("test {} has visualization, needs visVar".format(sec))
            invalid = 1

        if mysuite.sourceTree == "BoxLib" and mytest.testSrcTree == "":
            mysuite.log.warn("testSrcTree not set for BoxLib test {}".format(sec))
            invalid = 1


        # add the current test object to the master list
        if not invalid:
            test_list.append(mytest)
        else:
            mysuite.log.warn("test {} will be skipped".format(sec))


    # if any runs are parallel, make sure that the MPIcommand is defined
    any_MPI = any([t.useMPI for t in test_list])

    if any_MPI and mysuite.MPIcommand == "":
        mysuite.log.fail("ERROR: some tests are MPI parallel, but MPIcommand not defined")

    test_list.sort()

    return mysuite, test_list
