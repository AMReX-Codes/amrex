from __future__ import print_function

import datetime
import json
import os
import shutil
import sys
import test_util

DO_TIMINGS_PLOTS = True

try: import numpy as np
except: DO_TIMINGS_PLOTS = False

try: import matplotlib
except: DO_TIMINGS_PLOTS = False
else:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

try: import matplotlib.dates as dates
except: DO_TIMINGS_PLOTS = False


class Test(object):

    def __init__(self, name):

        self.name = name

        self.log = None

        self.buildDir = ""

        self.extra_build_dir = ""

        self.target = ""

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

        self.acc = 0
        
        self.useMPI = 0
        self.numprocs = -1

        self.useOMP = 0
        self.numthreads = -1

        self.doVis = 0
        self.visVar = ""

        self.doComparison = True

        self.analysisRoutine = ""
        self.analysisMainArgs = ""
        self.analysisOutputImage = ""

        self.png_file = None

        self.outputFile = ""
        self.compareFile = ""

        self.compare_file_used = ""

        self.diffDir = ""
        self.diffOpts = ""

        self.addToCompileString = ""

        self.runtime_params = ""

        self.reClean = 0    # set automatically, not by users

        self.wall_time = 0   # set automatically, not by users

        self.nlevels = None  # set but running fboxinfo on the output

        self.comp_string = None  # set automatically
        self.run_command = None  # set automatically

        self.job_info_field1 = ""
        self.job_info_field2 = ""
        self.job_info_field3 = ""

        self.has_jobinfo = 0  # filled automatically

        self.backtrace = []   # filled automatically

        self.has_stderr = False # filled automatically

        self.compile_successful = False  # filled automatically
        self.compare_successful = False  # filled automatically

        self.customRunCmd = None

        self.compareParticles = False
        self.particleTypes = ""

        self.keywords = []
        
    def __lt__(self, other):
        return self.value() < other.value()

    def value(self):
        return self.name

    def find_backtrace(self):
        """ find any backtrace files produced """
        return [ft for ft in os.listdir(self.output_dir)
                if os.path.isfile(ft) and ft.startswith("Backtrace.")]

    def get_last_plotfile(self, output_dir=None):
        """ Find the last plotfile written.  Note: we give an error if the
            last plotfile is 0.  If output_dir is specified, then we use
            that instead of the default
        """

        if output_dir is None:
            output_dir = self.output_dir   # not yet implemented

        plts = [d for d in os.listdir(output_dir) if \
                (os.path.isdir(d) and
                 d.startswith("{}_plt".format(self.name))) or \
                (os.path.isfile(d) and
                 d.startswith("{}_plt".format(self.name)) and d.endswith(".tgz"))]

        if len(plts) == 0:
            self.log.warn("test did not produce any output")
            return ""

        plts.sort()
        last_plot = plts.pop()

        if last_plot.endswith("00000"):
            self.log.warn("only plotfile 0 was output -- skipping comparison")
            return ""

        return last_plot


class Suite(object):

    def __init__(self, args):

        self.args = args

        # this will hold all of the Repo() objects for the AMReX, source,
        # and build directories
        self.repos = {}

        self.test_file_path = os.getcwd() + '/' + self.args.input_file[0]

        self.suiteName = "testDefault"
        self.sub_title = ""

        self.sourceTree = ""
        self.testTopDir = ""
        self.webTopDir = ""

        self.useCmake = 0

        # set automatically
        self.source_dir = ""
        self.source_build_dir ="" # Cmake build dir
        self.source_cmake_opts =""
        self.amrex_dir = ""
        self.amrex_install_dir = "" # Cmake installation dir
        self.amrex_cmake_opts = ""

        self.MPIcommand = ""
        self.MPIhost = ""

        self.FCOMP = "gfortran"
        self.COMP = "g++"

        self.add_to_f_make_command = ""
        self.add_to_c_make_command = ""

        self.summary_job_info_field1 = ""
        self.summary_job_info_field2 = ""
        self.summary_job_info_field3 = ""

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

        self.slack_post = 0
        self.slack_webhookfile = ""
        self.slack_webhook_url = None
        self.slack_channel = ""
        self.slack_username = ""

        self.globalAddToExecString = ""

        # this will be automatically filled
        self.extra_src_comp_string = ""

        # delete all plot/checkfiles but the plotfile used for comparison upon
        # completion
        self.purge_output = 0

        self.log = None

        self.do_timings_plots = DO_TIMINGS_PLOTS

        # default branch -- we use this only for display purposes --
        # if the test was run on a branch other than the default, then
        # an asterisk will appear next to the date in the main page
        self.default_branch = "master"

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

        # if we specified any keywords, only run those
        if self.args.keyword is not None:
            test_list = [t for t in test_list_old if self.args.keyword in t.keywords]
            
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
            if not self.args.make_benchmarks is None:
                os.mkdir(bench_dir)
            else:
                self.log.fail("ERROR: benchmark directory, {}, does not exist".format(bench_dir))
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
        full_web_dir = "{}/{}/".format(self.webTopDir, test_dir)

        if self.args.do_temp_run:
            if os.path.isdir(full_web_dir):
                shutil.rmtree(full_web_dir)

        os.mkdir(full_web_dir)

        # copy the test file into the web output directory
        shutil.copy(self.test_file_path, full_web_dir)

        self.test_dir = test_dir
        self.full_test_dir = full_test_dir
        self.full_web_dir = full_web_dir

    def get_run_history(self, active_test_list):
        """ return the list of output directories run over the
            history of the suite and a separate list of the tests
            run (unique names) """

        valid_dirs = []
        all_tests = []

        # start by finding the list of valid test directories
        for f in os.listdir(self.webTopDir):

            # look for a directory of the form 20* (this will work up until 2099
            if f.startswith("20") and os.path.isdir(f):

                # look for the status file
                status_file = f + '/' + f + '.status'
                if os.path.isfile(status_file):
                    valid_dirs.append(f)

        valid_dirs.sort()
        valid_dirs.reverse()

        # now find all of the unique problems in the test directories
        for d in valid_dirs:

            for f in os.listdir(self.webTopDir + d):
                if f.endswith(".status") and not f.startswith("20"):
                    index = f.rfind(".status")
                    test_name = f[0:index]

                    if all_tests.count(test_name) == 0:
                        if (not self.reportActiveTestsOnly) or (test_name in active_test_list):
                            all_tests.append(test_name)

        all_tests.sort()

        return valid_dirs, all_tests

    def make_timing_plots(self, active_test_list):
        """ plot the wallclock time history for all the valid tests """

        valid_dirs, all_tests = self.get_run_history(active_test_list)

        # store the timings in NumPy arrays in a dictionary
        timings = {}
        N = len(valid_dirs)
        for t in all_tests:
            timings[t] = np.zeros(N, dtype=np.float64)

        # now get the timings from the web output
        for n, d in enumerate(valid_dirs):
            for t in all_tests:
                ofile = "{}/{}/{}.html".format(self.webTopDir, d, t)
                try: f = open(ofile)
                except:
                    timings[t][n] = 0.0
                    continue

                found = False
                for line in f:
                    if "Execution time" in line:
                        found = True
                        # this is of the form: <li>Execution time: 412.930 s
                        timings[t][n] = float(line.split(":")[1].strip().split(" ")[0])
                        break

                    elif "(seconds)" in line:
                        found = True
                        # this is the older form -- split on "="
                        # form: <p><b>Execution Time</b> (seconds) = 399.414828
                        timings[t][n] = float(line.split("=")[1])
                        break

                f.close()
                if not found:
                    timings[t][n] = 0.0

        # make the plots
        for t in all_tests:
            _d = valid_dirs[:]
            _t = list(timings[t])

            days = []
            times = []
            for n, ttime in enumerate(_t):
                if not ttime == 0.0:
                    # sometimes the date is of the form YYYY-MM-DD-NNN, where NNN
                    # is the run -- remove that
                    date = _d[n]
                    if len(date) > 10:
                        date = date[:date.rfind("-")]

                    days.append(dates.datestr2num(date))
                    times.append(ttime)


            if len(times) == 0: continue

            plt.clf()
            plt.plot_date(days, times, "o", xdate=True)

            years = dates.YearLocator()   # every year
            months = dates.MonthLocator()
            years_fmt = dates.DateFormatter('%Y')

            ax = plt.gca()
            ax.xaxis.set_major_locator(years)
            ax.xaxis.set_major_formatter(years_fmt)
            ax.xaxis.set_minor_locator(months)

            plt.ylabel("time (seconds)")
            plt.title(t)

            if max(times)/min(times) > 10.0:
                ax.set_yscale("log")

            fig = plt.gcf()
            fig.autofmt_xdate()

            plt.savefig("{}/{}-timings.png".format(self.webTopDir, t))


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
                    if line.find("FAILED") >= 0 or line.find("CRASHED") >= 0:
                        failed.append(test)

        os.chdir(cwd)
        return failed

    def make_realclean(self, repo="source"):
        build_comp_string = ""
        if self.repos[repo].build == 1:
            if not self.repos[repo].comp_string is None:
                build_comp_string = self.repos[repo].comp_string

        extra_src_comp_string = ""
        if not self.extra_src_comp_string is None:
            extra_src_comp_string = self.extra_src_comp_string

        cmd = "{} AMREX_HOME={} {} {} realclean".format(
            self.MAKE, self.amrex_dir,
            extra_src_comp_string, build_comp_string)

        test_util.run(cmd)

    def build_f(self, test=None, opts="", target="", outfile=None):
        """ build an executable with the Fortran AMReX build system """

        build_opts = ""
        if test is not None:
            build_opts += "NDEBUG={} ".format(f_flag(test.debug, test_not=True))
            build_opts += "ACC={} ".format(f_flag(test.acc))
            build_opts += "MPI={} ".format(f_flag(test.useMPI))
            build_opts += "OMP={} ".format(f_flag(test.useOMP))

            if not test.extra_build_dir == "":
                build_opts += self.repos[test.extra_build_dir].comp_string + " "

            if not test.addToCompileString == "":
                build_opts += test.addToCompileString + " "

        all_opts = "{} {} {}".format(self.extra_src_comp_string, build_opts, opts)

        comp_string = "{} -j{} AMREX_HOME={} COMP={} {} {} {}".format(
            self.MAKE, self.numMakeJobs, self.amrex_dir,
            self.FCOMP, self.add_to_f_make_command, all_opts, target)

        self.log.log(comp_string)
        stdout, stderr, rc = test_util.run(comp_string, outfile=outfile)

        # make returns 0 if everything was good
        if not rc == 0:
            self.log.warn("build failed")

        return comp_string, rc


    def build_c(self, test=None, opts="", outfile=None):

        build_opts = ""

        if test is not None:
            build_opts += "DEBUG={} ".format(c_flag(test.debug))
            build_opts += "USE_ACC={} ".format(c_flag(test.acc))
            build_opts += "USE_MPI={} ".format(c_flag(test.useMPI))
            build_opts += "USE_OMP={} ".format(c_flag(test.useOMP))
            build_opts += "DIM={} ".format(test.dim)

            if not test.extra_build_dir == "":
                build_opts += self.repos[test.extra_build_dir].comp_string + " "

            if "source" in self.repos:
                if not self.repos["source"].comp_string is None:
                    build_opts += self.repos["source"].comp_string + " "
                
            if not test.addToCompileString == "":
                build_opts += test.addToCompileString + " "

        all_opts = "{} {} {}".format(self.extra_src_comp_string, build_opts, opts)

        comp_string = "{} -j{} AMREX_HOME={} {} COMP={} {}".format(
            self.MAKE, self.numMakeJobs, self.amrex_dir,
            all_opts, self.COMP, self.add_to_c_make_command)

        self.log.log(comp_string)
        stdout, stderr, rc = test_util.run(comp_string, outfile=outfile)

        # make returns 0 if everything was good
        if not rc == 0:
            self.log.warn("build failed")

        return comp_string, rc

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
        sout, serr, ierr = test_util.run(test_run_command, stdin=True,
                                         outfile="{}.run.out".format(test.name),
                                         errfile="{}.err.out".format(test.name),
                                         env=test_env)
        test.run_command = test_run_command

    def copy_backtrace(self, test):
        """
        if any backtrace files were output (because the run crashed), find them
        and copy them to the web directory
        """
        backtrace = test.find_backtrace()

        for btf in backtrace:
            ofile = "{}/{}.{}".format(self.full_web_dir, test.name, btf)
            shutil.copy(btf, ofile)
            test.backtrace.append("{}.{}".format(test.name, btf))


    def build_tools(self, test_list):

        self.log.skip()
        self.log.bold("building tools...")
        self.log.indent()

        self.tools = {}

        self.f_compare_tool_dir = "{}/Tools/Postprocessing/F_Src/".format(
            os.path.normpath(self.amrex_dir))

        os.chdir(self.f_compare_tool_dir)

        self.make_realclean(repo="AMReX")

        ftools = ["fcompare", "fboxinfo"]
        if any([t for t in test_list if t.dim == 2]): ftools.append("fsnapshot2d")
        if any([t for t in test_list if t.dim == 3]): ftools.append("fsnapshot3d")

        for t in ftools:
            self.log.log("building {}...".format(t))
            comp_string, rc = self.build_f(target="programs={}".format(t), opts="NDEBUG=t MPI= ")
            if not rc == 0:
                self.log.fail("unable to continue, tools not able to be built")

            exe = test_util.get_recent_filename(self.f_compare_tool_dir, t, ".exe")
            self.tools[t] = "{}/{}".format(self.f_compare_tool_dir, exe)

        self.c_compare_tool_dir = "{}/Tools/Postprocessing/C_Src/".format(
            os.path.normpath(self.amrex_dir))

        os.chdir(self.c_compare_tool_dir)

        ctools = ["particle_compare"]

        for t in ctools:
            self.log.log("building {}...".format(t))
            comp_string, rc = self.build_c(opts="NDEBUG=t MPI= ")
            if not rc == 0:
                self.log.fail("unable to continue, tools not able to be built")

            exe = test_util.get_recent_filename(self.c_compare_tool_dir, t, ".exe")

            self.tools[t] = "{}/{}".format(self.c_compare_tool_dir, exe)

        self.log.outdent()

    def slack_post_it(self, message):

        payload = {}

        # make sure there are no quotes in the strings
        payload["channel"] = self.slack_channel.replace('"', '')
        payload["username"] = self.slack_username.replace('"', '')
        payload["text"] = message.replace("'", "")  # apostrophes

        s = json.dumps(payload)
        cmd = "curl -X POST --data-urlencode 'payload={}' {}".format(s, self.slack_webhook_url)
        test_util.run(cmd)



    #######################################################
    #        CMake utilities                              #
    #######################################################
    def cmake_config( self, name, path, configOpts="",  install = 0, env = None):
        "Generate Cmake configuration"

        self.log.outdent()
        self.log.skip()
        self.log.bold("configuring " + name +  " build...")
        self.log.indent()

        # Setup dir names
        builddir   = path + 'builddir'
        if install:
            installdir = path + 'installdir'
        else:
            installdir = None

        # Define enviroment
        ENV = {}
        ENV =  dict(os.environ) # Copy of current enviroment
        ENV['FC']  = self.FCOMP
        ENV['CXX'] = self.COMP

        if env is not None: ENV.update(env)

        # remove build and installation directories if present and re-make them
        if os.path.isdir(builddir):
            shutil.rmtree(builddir)
        self.log.log("mkdir " + builddir)
        os.mkdir(builddir)

        if install:
            if os.path.isdir(installdir):
                shutil.rmtree(installdir)
            self.log.log("mkdir " + installdir)
            os.mkdir(installdir)

        # Logfile
        coutfile = '{}{}.cmake.log'.format( self.full_test_dir, name )

        # Run cmake
        cmd = 'cmake {} -H{} -B{} '.format(configOpts, path, builddir)
        if install:
            cmd += '-DCMAKE_INSTALL_PREFIX:PATH='+installdir

        self.log.log(cmd)
        stdout, stderr, rc = test_util.run(cmd, outfile=coutfile, env=ENV)

        # Check exit condition
        if not rc == 0:
            errstr  = "\n \nERROR! Cmake configuration failed for " + name + " \n"
            errstr += "Check " + coutfile + " for more information."
            sys.exit(errstr)

        return builddir, installdir


    def cmake_clean( self, name, path ):
        "Clean Cmake build and install directories"

        self.log.outdent()
        self.log.skip()
        self.log.bold("cleaning " + name +  " Cmake directories...")
        self.log.indent()

        # Setup dir names
        builddir   = path + 'builddir'
        installdir = path + 'installdir'

        # remove build and installation directories if present
        if os.path.isdir(builddir):
            shutil.rmtree(builddir)

        if os.path.isdir(installdir):
                shutil.rmtree(installdir)

        return




    def cmake_build( self, name, target, path, opts = '', env = None, outfile = None ):
        "Build target for a repo configured via cmake"

        self.log.outdent()
        self.log.skip()
        self.log.bold("building " + name +  "...")
        self.log.indent()

        # Set enviroment
        ENV =  dict(os.environ) # Copy of current enviroment
        if env is not None: ENV.update(env)

        if outfile is not None:
            coutfile = outfile
        else:
            coutfile = '{}{}.{}.make.log'.format( self.full_test_dir, name, target )

        cmd = '{} -j{} {} {}'.format( self.MAKE, self.numMakeJobs, opts, target )
        self.log.log(cmd)
        stdout, stderr, rc = test_util.run(cmd, outfile=coutfile, cwd=path, env=ENV )

        # make returns 0 if everything was good
        if not rc == 0:
            self.log.warn("build failed")

        comp_string = cmd

        return rc, comp_string



    def build_test_cmake(self, test, opts="",  outfile=None):
        """ build an executable with CMake build system """

        env = {"AMREX_HOME":self.amrex_install_dir}

        rc, comp_string = self.cmake_build( name    = test.name,
                                            target  = test.target,
                                            path    = self.source_build_dir,
                                            opts    = opts,
                                            env     = env,
                                            outfile = outfile)

        # make returns 0 if everything was good
        if rc == 0:
            # Find location of executable
            for root, dirnames, filenames in os.walk(self.source_build_dir):
                if test.target in filenames:
                    path_to_exe = os.path.join(root, test.target)
                    break

            # Copy and rename executable to test dir
            shutil.move("{}".format(path_to_exe),
                        "{}/{}/{}.ex".format(self.source_dir,test.buildDir,test.name))
        else:
          self.log.warn("build failed")


        return comp_string, rc



def f_flag(opt, test_not=False):
    """ convert a test parameter into t if true for the Fortran build system """
    if test_not:
        if opt: return " "
        else: return "t"
    else:
        if opt: return "t"
        else: return " "

def c_flag(opt, test_not=False):
    """ convert a test parameter into t if true for the Fortran build system """
    if test_not:
        if opt: return "FALSE"
        else: return "TRUE"
    else:
        if opt: return "TRUE"
        else: return "FALSE"
