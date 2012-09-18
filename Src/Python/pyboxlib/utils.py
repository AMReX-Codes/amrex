"""Various PyBoxLib utilities."""

import os

from pyboxlib.plotfile import compare

# def has_plotfiles(dirs):
#     for d in dirs:
#         if d.find('plt') > -1:
#             return True
#     return False

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
        


def auto_errors(name, params, reference):

    # find last plotfile of reference
    for root, dirs, files in os.walk(os.path.join(name, reference)):
        if 'probin.nml' in files:
            last_plotfile = sorted(dirs)[-1]
            reference = os.path.join(root, last_plotfile)

    # find runs and compute errors
    errors = {}

    for root, dirs, files in os.walk(name):
        if 'probin.nml' in files:
            probin = Probin(os.path.join(root, 'probin.nml'))
            probin.parse()

            last_plotfile = sorted(dirs)[-1]
            errs, _ = compare(reference, os.path.join(root, last_plotfile))

            values = tuple([ probin.params[p] for p in params ])
            errors[values] = (max([ errs[x][0] for x in errs ]), max([ errs[x][1] for x in errs ]))

    return errors


def plot_convergence(errors, params, order=None):

    import collections
    import matplotlib.pylab as plt

    x = collections.defaultdict(list)
    y = collections.defaultdict(list)

    for key in sorted(errors):
        k = key[:-1]
        x[k].append(key[-1])
        y[k].append(errors[key][1]) # relative error

    for k in sorted(y):
        label = [ params[i] + ' ' + str(k[i]) for i in range(len(k)) ]
        label = ', '.join(label)

        plt.loglog(x[k], y[k], label=label, marker='o')

        if order:
            e0 = y[k][0]
            xx = [ x[k][0], x[k][-1] ]
            if len(k) == 1:
                kk = k[0]
            else:
                kk = k
            yy = [ e0, e0*(x[k][-1]/x[k][0])**order[kk] ]
            plt.loglog(xx, yy, ':k', label=None)

    plt.xlabel(params[-1])
    plt.ylabel('Maximum relative error')
    plt.legend(loc='best')
