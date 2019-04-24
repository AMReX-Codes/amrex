.. role:: cpp(code)
   :language: c++

Continuous Compilation Testing
==============================

As a first line of testing, on every commit to the repository, we verify that we can compile
AMReX as a library for a common set of configuration options. This operation is performed
through Travis-CI. This layer of testing is deliberately limited, so that it can be run
quickly on every commit. For more extensive testing, we rely on the nightly regression results.

              
Nightly Regression Testing
==========================

Each night, we automically run a suite of tests, both on AMReX itself, and on a most of the major
application codes that use it as a framework. We use an in-house test runner script to manage this
operation, originally developed by Michael Zingale for the Castro code, and later expanded to other
application codes as well. The results for each night are collected and stored on a web page; see
https://ccse.lbl.gov/pub/RegressionTesting/ for the latest set of results.

Running the test suite locally
==============================

The test suite is mostly used internally by AMReX developers. However,
f you are making a pull request to AMReX, it can be useful to run the test suite
on your local machine, to reduce the liklihood that your changes break some existing functionality.
To run the test suite on locally, you must first obtain a copy of the test runner source, available
on Github here: https://github.com/AMReX-Codes/regression_testing. The test runner requires Python
version 2.7 or greater.

Next, you need a configuration file that defines which tests to run, which amrex repository to test,
which branch to use, etc. A sample configuration file for AMReX is distributed with the amrex source
code at :cpp:`amrex/Tools/RegressionTesting/AMReX-tests.ini`. You will need to modify a few of the entries
to, for example, point the test runner to the clone of amrex on your local machine. Entries you will
likely want to change include:

.. highlight:: bash

::

   testTopDir = /path/to/test/output # the tests results and benchmarks will stored here
   webTopDir  = /path/to/web/output  # a web page with the test results will be written here

to control where the generated output will be written, and

::
   
   [AMReX]
   dir = /path/to/amrex  # the path to the amrex repository you want to test
   branch = "development"

to control which repository and branch to test.

The test runner is a Python script and can be invoked like so:

::

   python regtest.py <options> AMReX-Tests.ini
   
Before you can use it, you must first generate a set of "benchmarks" - i.e. known "good" answers to the
tests that will be run. If you are testing a pull request, you can generate these by running the script
with the a recent version of the :cpp:`development` branch of AMReX. You can generate the benchmarks like so:

::
   
   python regtest.py --make_benchmarks 'generating initial benchmarks' AMReX-Tests.ini

Once that is finished, you can switch over to the branch you want to test in :cpp:`AMReX-Tests.ini`, and then
re-run the script without the :cpp:`--make_benchmarks` option:

::
   
   python regtest.py --make_benchmarks 'generating initial benchmarks' AMReX-Tests.ini

The script will generate a set of html pages in the directory specified in your :cpp:`AMReX-Tests.ini`
file that you can examine using the browser of your choice.

For a complete set of script options, run

::

   python regtest.py --help

A particularly useful option lets you run just a subset of the complete test suite. To run only one test, you can do:
   
::
   
   python regtest.py --single_test <TestName> AMReX-Tests.ini

To run an enumerated list of tests, do:
   
::
   
   python regtest.py --tests '<TestName1> <TestName2> <TestName3>' AMReX-Tests.ini


Adding a new test
=================

New tests can be added to the suite by modifying the :cpp:`AMReX-Tests.ini` file. The easiest thing to
do is start from an existing test and modify it. For example, this entry:

::

   [MLMG_FI_PoisCom] 
   buildDir = Tutorials/LinearSolvers/ABecLaplacian_F
   inputFile = inputs-rt-poisson-com
   dim = 3
   restartTest = 0
   useMPI = 1
   numprocs = 2
   useOMP = 1
   numthreads = 3
   compileTest = 0
   doVis = 0
   outputFile = plot
   testSrcTree = C_Src

defines a test called :cpp:`MLMG_FI_PoisCom` by specifying the apppropriate build directory, inputs file,
and a set of configuration options. The above options are the most commonly changed; for a full list
of options, see the example configuration file at https://github.com/AMReX-Codes/regression_testing/blob/master/example-tests.ini.


