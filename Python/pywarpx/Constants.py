# Copyright 2018-2019 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket

class Constants(Bucket):
    """
    The purpose of this class is to be hold user defined constants
    """
    def __init__(self):
        Bucket.__init__(self, 'my_constants')

    def __setattr__(self, name, value):
        # Make sure that any constants redefined have a consistent value
        if name in self.argvattrs:
            assert self.argvattrs[name] == value, Exception('An consistent values given for user defined constants')
        Bucket.__setattr__(self, name, value)

my_constants = Constants()
