# This script modifies `WarpX-test.ini` (which is used for nightly builds)
# and creates the file `shippable-test.ini` (which is used for continous
# integration on Shippable (https://app.shippable.com/))
import re

with open('WarpX-tests.ini') as f:
    text = f.read()

# Add doComparison = 0 for each test
text = re.sub( '\[(?P<name>.*)\]\nbuildDir = ',
               '[\g<name>]\ndoComparison = 0\nbuildDir = ', text )

# Use ccache in order to save compilation time
text = re.sub( 'addToCompileString = (?P<strings>.*)',
               'addToCompileString = USE_CCACHE=TRUE \g<strings>', text )

# Use only 2 cores for compiling
text = re.sub( 'numMakeJobs = \d+', 'numMakeJobs = 2', text )

# Remove Python test (does not compile in the Docker container)
text = re.sub( '\[Python_Langmuir\]\n(.+\n)*', '', text)

# Remove Langmuir_x/y/z test (too long; not that useful)
text = re.sub( '\[Langmuir_[xyz]\]\n(.+\n)*', '', text)
# Skip unit tests (too long; not that useful)
text = re.sub( '\[UnitTest_[a-zA-Z]+\]\n(.+\n)*', '', text)

# Prevent emails from being sent
text = re.sub( 'sendEmailWhenFail = 1', 'sendEmailWhenFail = 0', text )

with open('shippable-tests.ini', 'w') as f:
    f.write(text)
