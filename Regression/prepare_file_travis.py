# Copyright 2018-2019 Andrew Myers, Luca Fedeli, Maxence Thevenet
# Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script modifies `WarpX-test.ini` (which is used for nightly builds)
# and creates the file `travis-test.ini` (which is used for continous
# integration on TravisCI (https://travis-ci.org/)
# The subtests that are selected are controlled by WARPX_TEST_DIM
# The architecture (CPU/GPU) is selected by WARPX_TEST_ARCH
import re
import os
# Get relevant environment variables
arch = os.environ.get('WARPX_TEST_ARCH', 'CPU')

ci_regular_cartesian = os.environ.get('WARPX_CI_REGULAR_CARTESIAN') == 'TRUE'
ci_psatd = os.environ.get('WARPX_CI_PSATD') == 'TRUE'
ci_python_main = os.environ.get('WARPX_CI_PYTHON_MAIN') == 'TRUE'
ci_single_precision = os.environ.get('WARPX_CI_SINGLE_PRECISION') == 'TRUE'
ci_rz_or_nompi = os.environ.get('WARPX_CI_RZ_OR_NOMPI') == 'TRUE'
ci_qed = os.environ.get('WARPX_CI_QED') == 'TRUE'
ci_openpmd = os.environ.get('WARPX_CI_OPENPMD') == 'TRUE'

# Find the directory in which the tests should be run
current_dir = os.getcwd()
test_dir = re.sub('warpx/Regression', '', current_dir )

with open('WarpX-tests.ini') as f:
    text = f.read()

# Replace default folder name
text = re.sub('/home/regtester/AMReX_RegTesting', test_dir, text)
# Remove the web directory
text = re.sub('[\w\-\/]*/web', '', text)

# Add doComparison = 0 for each test
text = re.sub( '\[(?P<name>.*)\]\nbuildDir = ',
               '[\g<name>]\ndoComparison = 0\nbuildDir = ', text )

# Change compile options when running on GPU
if arch == 'GPU':
    text = re.sub( 'addToCompileString =',
                    'addToCompileString = USE_GPU=TRUE USE_OMP=FALSE USE_ACC=TRUE ', text)
    text = re.sub( 'COMP\s*=.*', 'COMP = pgi', text )
print('Compiling for %s' %arch)

# Extra dependencies
if ci_openpmd:
    text = re.sub('addToCompileString =',
                  'addToCompileString = USE_OPENPMD=TRUE ', text)

# Add runtime option: crash for unused variables
text = re.sub('runtime_params =',
              'runtime_params = amrex.abort_on_unused_inputs=1 ',
              text)

# Use only 2 cores for compiling
text = re.sub( 'numMakeJobs = \d+', 'numMakeJobs = 2', text )
# Use only 1 OMP thread for running
text = re.sub( 'numthreads = \d+', 'numthreads = 1', text)
# Prevent emails from being sent
text = re.sub( 'sendEmailWhenFail = 1', 'sendEmailWhenFail = 0', text )

# Remove Python test (does not compile)
text = re.sub( '\[Python_Langmuir\]\n(.+\n)*', '', text)

# Remove Langmuir_x/y/z test (too long; not that useful)
text = re.sub( '\[Langmuir_[xyz]\]\n(.+\n)*', '', text)

# Select the tests to be run
# --------------------------

# - Extract test blocks (they are identified by the fact that they contain "inputFile")
select_test_regex = r'(\[(.+\n)*inputFile(.+\n)*)'
test_blocks =  [ match[0] for match in re.findall(select_test_regex, text) ]
# - Remove the test blocks from `text` (only the selected ones will be added back)
text = re.sub( select_test_regex, '', text )

def select_tests(blocks, match_string_list, do_test):
    """Remove or keep tests from list in WarpX-tests.ini according to do_test variable"""
    if do_test not in [True, False]:
        raise ValueError("do_test must be True or False")
    if (do_test == False):
        for match_string in match_string_list:
            print('Selecting tests without ' + match_string)
            blocks = [ block for block in blocks if not match_string in block ]
    else:
        for match_string in match_string_list:
            print('Selecting tests with ' + match_string)
            blocks = [ block for block in blocks if match_string in block ]
    return blocks

if ci_regular_cartesian:
    test_blocks = select_tests(test_blocks, ['USE_RZ=TRUE'], False)
    test_blocks = select_tests(test_blocks, ['USE_PSATD=TRUE'], False)
    test_blocks = select_tests(test_blocks, ['PYTHON_MAIN=TRUE'], False)
    test_blocks = select_tests(test_blocks, ['PRECISION=FLOAT', 'USE_SINGLE_PRECISION_PARTICLES=TRUE'], False)
    test_blocks = select_tests(test_blocks, ['useMPI = 0'], False)
    test_blocks = select_tests(test_blocks, ['QED=TRUE'], False)

if ci_psatd:
    test_blocks = select_tests(test_blocks, ['USE_PSATD=TRUE'], True)
    # Remove PSATD single-precision, which is done in ci_single_precision
    test_blocks = select_tests(test_blocks, ['PRECISION=FLOAT'], False)

if ci_python_main:
    test_blocks = select_tests(test_blocks, ['PYTHON_MAIN=TRUE'], True)

if ci_single_precision:
    test_blocks = select_tests(test_blocks, ['PRECISION=FLOAT', 'USE_SINGLE_PRECISION_PARTICLES=TRUE'], True)

if ci_rz_or_nompi:
    test_blocks = select_tests(test_blocks, ['PYTHON_MAIN=TRUE'], False)
    block1 = select_tests(test_blocks, ['USE_RZ=TRUE'], True)
    block2 = select_tests(test_blocks, ['useMPI = 0'], True)
    test_blocks = block1 + block2

if ci_qed:
    test_blocks = select_tests(test_blocks, ['QED=TRUE'], True)

# - Add the selected test blocks to the text
text = text + '\n' + '\n'.join(test_blocks)

with open('travis-tests.ini', 'w') as f:
    f.write(text)
