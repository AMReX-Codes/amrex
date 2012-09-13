"""Various Fabric (fabfile.org) utilities for launching and managing
BoxLib runs.

These routines assume

"""

import glob
import os

from fabric.api import *
from fabric.colors import *
from pyboxlib.utils import *


def submit(exe, dname, probin):

    local('cd %s && %s %s' % (dname, exe, probin))


def mkdir(dname):

    try:
        os.mkdir(dname)
    except OSError:
        pass


def fstr(x):
    """Convert x to a string, appending 'd0' to floating point
    numbers."""

    if not isinstance(x, float):
        return str(x)

    s = str(x)
    s = s.replace('e', 'd')
    if s.find('d') == -1:
        s = s + 'd0'

    return s


def convergence(name, param, values, **params):
    """Perform a convergence run..."""

    mkdir(name)

    exe = glob.glob(env.exe)[0]
    exe = os.path.abspath(exe)
    
    # read probin template
    with open(env.probin, 'r') as f:
        probin = f.read()

    # run convergence tests
    for value in values:
        dname = name + '/' + param + str(value)
        mkdir(dname)

        # write probin
        with open(dname + '/inputs', 'w') as f:
            params[param] = fstr(value)
            f.write(probin.format(**params))

        # submit
        submit(exe, dname, 'inputs')

            

                


    
    

    
