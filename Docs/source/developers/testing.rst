.. _developers-testing:

Testing the code
================

When adding a new feature, you want to make sure that (i) you did not break the existing code and (ii) your contribution gives correct results. While existing capabilities are tested regularly remotely (when commits are pushed to an open PR on TravisCI, and every night on local clusters), it can also be useful to run tests on your custom input file. This section details how to use both automated and custom tests.

Continuous Integration in WarpX
-------------------------------

WarpX-tests.ini files
^^^^^^^^^^^^^^^^^^^^^

.. note::
   section empty!

TravisCI tests
^^^^^^^^^^^^^^

.. note::
   section empty!

Local tests every night
^^^^^^^^^^^^^^^^^^^^^^^

.. note::
   section empty!

Run the test suite locally
--------------------------

Once your new feature is ready, you can check that you did not break anything. WarpX has automated tests running for each Pull Request. For easier debugging, it can be convenient to run the tests on your local machine with

.. code-block:: sh

    ./run_test.sh

from WarpX root folder. It will compile all necessary executables and run all tests. The tests can be influenced by environment variables:

* ``export WARPX_TEST_DIM=3``, ``export WARPX_TEST_DIM=2`` or ``export WARPX_TEST_DIM=RZ`` in order to select only the tests that correspond to this dimensionality
* ``export WARPX_TEST_ARCH=CPU`` or ``export WARPX_TEST_ARCH=GPU`` in order to run the tests on CPU or GPU respectively.
* ``export WARPX_TEST_COMMIT=...`` in order to test a specific commit.

The command above (without command line arguments) runs all the tests defined in [Regression/WarpX-tests.ini](./Regression/WarpX-tests.ini). In order to run single tests, pass the test names as command line arguments:

.. code-block:: sh

    ./run_test.sh test1 test2
    # For instance
    ./run_test.sh PlasmaAccelerationBoost3d Larmor

Add a test to the suite
-----------------------

There are three steps to follow to add a new automated test (illustrated here for PML boundary conditions):

* An input file for your test, in folder `Example/Tests/...`. For the PML test, the input file is at ``Examples/Tests/PML/inputs_2d``. You can also re-use an existing input file (even better!) and pass specific parameters at runtime (see below).
* A Python script that reads simulation output and tests correctness versus theory or calibrated results. For the PML test, see ``Examples/Tests/PML/analysis_pml_yee.py``. It typically ends with Python statement `assert( error<0.01 )`.
* Add an entry to [Regression/WarpX-tests.ini](./Regression/WarpX-tests.ini), so that a WarpX simulation runs your test in the continuous integration process on [Travis CI](https://docs.travis-ci.com/user/tutorial/), and the Python script is executed to assess the correctness. For the PML test, the entry is

.. code-block::

   [pml_x_yee]
   buildDir = .
   inputFile = Examples/Tests/PML/inputs2d
   runtime_params = warpx.do_dynamic_scheduling=0 algo.maxwell_fdtd_solver=yee
   dim = 2
   addToCompileString =
   restartTest = 0
   useMPI = 1
   numprocs = 2
   useOMP = 1
   numthreads = 2
   compileTest = 0
   doVis = 0
   analysisRoutine = Examples/Tests/PML/analysis_pml_yee.py

If you re-use an existing input file, you can add arguments to ``runtime_params``, like ``runtime_params = amr.max_level=1 amr.n_cell=32 512 max_step=100 plasma_e.zmin=-200.e-6``.

Useful tool for plotfile comparison: ``fcompare``
--------------------------------------------------

AMReX provides ``fcompare``, an executable that takes two ``plotfiles`` as input and returns the absolute and relative difference for each field between these two plotfiles. For some changes in the code, it is very convenient to run the same input file with an old and your current version, and ``fcompare`` the plotfiles at the same iteration. To use it:

.. code-block:: sh

   # Compile the executable
   cd <path to AMReX>/Tools/Plotfile/ # This may change
   make -j 8
   # Run the executable to compare old and new versions
   <path to AMReX>/Tools/Plotfile/fcompare.gnu.ex old/plt00200 new/plt00200

which should return something like

.. code-block:: sh

             variable name             absolute error            relative error
                                          (||A - B||)         (||A - B||/||A||)
   ----------------------------------------------------------------------------
   level = 0
   jx                                 1.044455105e+11               1.021651316
   jy                                  4.08631977e+16               7.734299273
   jz                                 1.877301764e+14               1.073458933
   Ex                                 4.196315448e+10               1.253551615
   Ey                                 3.330698083e+12               6.436470137
   Ez                                 2.598167798e+10              0.6804387128
   Bx                                     273.8687473               2.340209782
   By                                     152.3911863                1.10952567
   Bz                                     37.43212767                 2.1977289
   part_per_cell                                   15                    0.9375
   Ex_fp                              4.196315448e+10               1.253551615
   Ey_fp                              3.330698083e+12               6.436470137
   Ez_fp                              2.598167798e+10              0.6804387128
   Bx_fp                                  273.8687473               2.340209782
   By_fp                                  152.3911863                1.10952567
   Bz_fp                                  37.43212767                 2.1977289
