from __future__ import print_function

import argparse
import os
import shlex
import subprocess
import sys

usage = """
The test suite and its tests are defined through an input file in an INI
configuration file format.

The "main" block specifies the global test suite parameters:

  [main]

  testTopDir     = < full path to test output directory >
  webTopDir      = < full path to test web output directory >

  sourceTree = < C_Src, F_Src, or BoxLib -- what type is it? >

  suiteName = < descriptive name (i.e. Castro) >

  reportActiveTestsOnly = <0: web shows every test ever run;
                           1: just current tests >

  goUpLink = <1: add "Go UP" link at top of the web page >

  FCOMP = < name of Fortran compiler >
  COMP  = < name of C/C++ compiler >

  add_to_f_make_command = < any additional defines to add to the make invocation for F_Src BoxLib >
  add_to_c_make_command = < any additional defines to add to the make invocation for C_Src BoxLib >

  purge_output = <0: leave all plotfiles in place;
                  1: delete plotfiles after compare >

  MAKE = < name of make >
  numMakeJobs = < number of make jobs >

  MPIcommand = < MPI run command, with holders for host, # of proc, command >

     This should look something like:

          mpiexec -host @host@ -n @nprocs@ @command@ >

  MPIhost = < host for MPI job -- depends on MPI implementation >

  sendEmailWhenFail = < 1: send email when any tests fail >

  emailTo = < list of email addresses separated by commas, such as,
              foo@example.com, bar@example.com >

  emailBody = < email body >


The source git repositories are defined in separate blocks.  There
will always be a "BoxLib" block, and usually a "source" block which is
the default directory used for compiling the tests.  Any extra repos
(including those where additional tests are to be build) are defined
in their own block starting with "extra-"

The general form is:

  [name]

  dir = < full path to git repo >

  branch = < desired branch in the git repo >

  build = < 1: this is a directory that tests will be compiled in >

  comp_string = < a string that is added to the make line >

      comp_string can refer to both the main source directory (as @source@)
      and its own directory (as @self@), for example:

      comp_string = CASTRO_DIR=@source@ WDMERGER_HOME=@self@


Each test is given its own block, with the general form:

  [Sod-x]

  buildDir = < relative path (from sourceDir) for this problem >

  inputFile = < input file name >
  probinFile = < probin file name >

  dim = < dimensionality: 1, 2, or 3 >

  aux?File = < name of additional file needed by the test >
  link?File = < name of additional file needed by the test >

      Here "?" is 1, 2, or 3, allowing for several files per test

  restartTest = < is this a restart test? 0 for no, 1 for yes >
  restartFileNum = < # of file to restart from (if restart test) >

  useMPI = <is this a parallel (MPI) job? 0 for no, 1 for yes) >
  numprocs = < # of processors to run on (if parallel job) >

  useOMP = <is this an OpenMP job? 0 for no, 1 for yes) >
  numthreads = < # of threads to us with OpenMP (if OpenMP job) >

  acc = < 0 for normal run, 1 if we want OpenACC >

  debug = < 0 for normal run, 1 if we want debugging options on >

  compileTest = < 0 for normal run, 1 if we just test compilation >

  selfTest = < 0 for normal run, 1 if test self-diagnoses if it succeeded >
  stSuccessString = < string to find in self-test output to determine success >

  doVis = < 0 for no visualization, 1 if we do visualization >
  visVar = < string of the variable to visualize >

  analysisRoutine = < name of the script to run on the output >

      The script is run as:

        analysisRoutine [options] plotfile

  analysisMainArgs = < commandline arguments to pass to the analysisRoutine --
                       these should refer to options from the [main] block >

  analysisOutputImage = < name on analysis result image to show on web page >

  compareFile = < explicit output file to do the comparison with -- this is
                  assumed to be prefixed with the test name when output by
                  the code at runtime, e.g. test_plt00100 >

  outputFile = < explicit output file to compare with -- exactly as it will
                 be written.  Not prefix of the test name will be done >

  diffDir = < directory/file to do a plain text diff on (recursive, if dir) >

  diffOpts = < options to use with the diff command for the diffDir comparison >


Getting started:

To set up a test suite, it is probably easiest to write the
testfile.ini as described above and then run the test routine with the
--make_benchmarks option to create the benchmark directory.
Subsequent runs can be done as usual, and will compare to the newly
created benchmarks.  If differences arise in the comparisons due to
(desired) code changes, the benchmarks can be updated using
--make_benchmarks to reflect the new ``correct'' solution.

"""


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

    def warn(self, warn_msg):
        """
        output a warning.  It is always prefix with 'WARNING:'
        For multi-line warnings, send in a list of strings
        """
        prefix = self.indent_str + "WARNING: "
        filler = self.indent_str + "         "
        if isinstance(warn_msg, list):
            msg = [prefix + warn_msg[0]] + [filler + x for x in warn_msg[1:]]
            omsg = "\n".join(msg).strip()
        else:
            omsg = prefix + warn_msg
        nstr = self.warn_color + omsg + self.end_color
        print(nstr)

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


def get_args(arg_string=None):
    """ parse the commandline arguments.  If arg_string is present, we
        parse from there, otherwise we use the default (sys.argv) """

    parser = argparse.ArgumentParser(description=usage,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

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
    parser.add_argument("--with_valgrind", action="store_true",
                        help="run with valgrind")
    parser.add_argument("--valgrind_options", type=str, default="--leak-check=yes --log-file=vallog.%p",
                        help="valgrind options", metavar="'valgrind options'")
    parser.add_argument("--boxlib_git_hash", type=str, default=None, metavar="hash",
                        help="git hash of a version of BoxLib.  If provided, this version will be used to run tests.")
    parser.add_argument("--source_git_hash", type=str, default=None, metavar="hash",
                        help="git hash of a version of the main source code.  For BoxLib tests, this will be ignored.")
    parser.add_argument("--extra_src_git_hash", type=str, default=None, metavar="hash",
                        help="git hash of a version of the extra source code.  For BoxLib tests, this will be ignored.")
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


def run(string, stdin=False, outfile=None, store_command=False, env=None,
        outfile_mode="a", errfile=None, log=None):

    # shlex.split will preserve inner quotes
    prog = shlex.split(string)
    sin = None
    if stdin: sin = subprocess.PIPE
    p0 = subprocess.Popen(prog, stdin=sin, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, env=env)

    stdout0, stderr0 = p0.communicate()
    if stdin: p0.stdin.close()
    rc = p0.returncode
    p0.stdout.close()
    p0.stderr.close()

    if outfile is not None:
        try: cf = open(outfile, outfile_mode)
        except IOError:
            log.fail("  ERROR: unable to open file for writing")
        else:
            if store_command:
                cf.write(string)
            for line in stdout0:
                cf.write(line)

            if errfile is None:
                for line in stderr0:
                    cf.write(line)

            cf.close()

    if errfile is not None and stderr0 is not None:
        write_err = True
        if isinstance(stderr0, str):
            if stderr0.strip() == "":
                write_err = False
        if write_err:
            try: cf = open(errfile, outfile_mode)
            except IOError:
                log.fail("  ERROR: unable to open file for writing")
            else:
                for line in stderr0:
                    cf.write(line)
                cf.close()

    return stdout0, stderr0, rc


def get_recent_filename(fdir, base, extension):
    """ find the most recent file matching the base and extension """

    files = [f for f in os.listdir(fdir) if (f.startswith(base) and
                                             f.endswith(extension))]

    files.sort(key=lambda x: os.path.getmtime(x))

    try: return files.pop()
    except: return None
