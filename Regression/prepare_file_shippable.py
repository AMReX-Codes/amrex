# TODO: This file prepares the .ini file for tests on shippable
import re

with open('WarpX-tests.ini') as f:
    text = f.read()

# Add doCompare=0 for each test
text = re.sub( '\[(?P<name>.*)\]\nbuildDir = ',
               '[\g<name>]\ndoCompare=0\nbuildDir = ', text )

# Remove Python test (does not compile in the Docker container)
text = re.sub( '\[Python_Langmuir\]\n(.+\n)*', '', text)

# Prevent emails from being sent
text = re.sub( 'sendEmailWhenFail = 1', 'sendEmailWhenFail = 0', text )

# Modify paths (which are customized for the test server at CRD)
text = re.sub( 'regtester/AMReX_RegTesting/', '', text )

with open('shippable-tests.ini', 'w') as f:
    f.write(text)
