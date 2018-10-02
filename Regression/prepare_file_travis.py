# This script modifies `WarpX-test.ini` (which is used for nightly builds)
# and creates the file `travis-test.ini` (which is used for continous
# integration on TravisCI (https://travis-ci.org/)
# The subtests that are selected are controlled by the environement WARPX_DIM
import re
import os
# Get relevant environment variables
dim = os.environ.get('WARPX_DIM', None)

with open('WarpX-tests.ini') as f:
    text = f.read()

# Replace default folder name
text = re.sub('/home/regtester/AMReX_RegTesting/', '/home/travis/', text)

# Add doComparison = 0 for each test
text = re.sub( '\[(?P<name>.*)\]\nbuildDir = ',
               '[\g<name>]\ndoComparison = 0\nbuildDir = ', text )

# Use only 2 cores for compiling
text = re.sub( 'numMakeJobs = \d+', 'numMakeJobs = 2', text )

# Remove Python test (does not compile)
text = re.sub( '\[Python_Langmuir\]\n(.+\n)*', '', text)

# Remove Langmuir_x/y/z test (too long; not that useful)
text = re.sub( '\[Langmuir_[xyz]\]\n(.+\n)*', '', text)
# Skip unit tests (too long; not that useful)
text = re.sub( '\[UnitTest_[a-zA-Z]+\]\n(.+\n)*', '', text)

# Remove tests that do not have the right dimension
if dim is not None:
    print('Selecting tests with dim = %s' %dim)
    text = re.sub('\[.+\n(.+\n)*dim = [^%s]\n(.+\n)*' %dim, '', text)

# Prevent emails from being sent
text = re.sub( 'sendEmailWhenFail = 1', 'sendEmailWhenFail = 0', text )

with open('travis-tests.ini', 'w') as f:
    f.write(text)
