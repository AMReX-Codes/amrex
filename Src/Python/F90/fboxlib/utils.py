"""Various PyBoxLib utilities."""

from __future__ import with_statement

import os
import re

import numpy as np

from numpy import isnan
from pyboxlib.plotfile import compare


###############################################################################
# helpers

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


def _key(k):
    if len(k) == 1:
        kk = k[0]
    else:
        kk = k
    return kk


###############################################################################
# probin

class Probin(object):

    def __init__(self, fname):
        self.params = {}

        f = open(fname, 'r')
        self.probin = f.read()
        f.close()

    def update(self, **params):
        self.params.update(**params)

    def write(self, fname):
        params = {}

        for param in self.params:
            params[param] = fstr(self.params[param])

        f = open(fname, 'w')
        f.write(self.probin.format(**params))
        f.close()

    def parse(self):

        for line in self.probin.splitlines():
            try:
                param, value = map(lambda x: x.strip(), line.split('='))
            except:
                continue
            
            try:
                self.params[param] = int(value)
                continue
            except:
                pass

            try:
                v = value.replace('d', 'e')
                self.params[param] = float(v)
                continue
            except:
                pass

            self.params[param] = value
        

###############################################################################
# auto compare

def auto_compare(name, params, reference):

    # find last plotfile of reference
    for root, dirs, files in os.walk(reference):
        if 'probin.nml' in files:
            last_plotfile = sorted(dirs)[-1]
            reference = os.path.join(root, last_plotfile)

    # find runs and compute errors etc
    comps = []

    for root, dirs, files in os.walk(name):
        if 'probin.nml' in files:
            probin = Probin(os.path.join(root, 'probin.nml'))
            probin.parse()

            with open(os.path.join(root, 'out'), 'r') as f:
                out = f.read()

            m = re.search('run time =\s*(\S*)', out, flags=re.I)
            if m:
                runtime = float(m.group(1))
            else:
                runtime = 0.0

            last_plotfile = sorted(dirs)[-1]
            errs, _ = compare(reference, os.path.join(root, last_plotfile))

            comp = { 'abs': max([ errs[x][0] for x in errs ]), 
                     'rel': max([ errs[x][1] for x in errs ]),
                     'run': runtime }

            for p in params:
                comp[p[0]] = probin.params[p[0]]
                comp[p[1]] = probin.params[p[0]]

            comps.append(comp)

    return comps

def to_xy(series):

    series = [ s for s in series if not isnan(s[0]) and not isnan(s[1]) ]

    x = np.asarray([ s[0] for s in series ])
    y = np.asarray([ s[1] for s in series ])

    idx = np.argsort(x)
    
    return x[idx], y[idx]
