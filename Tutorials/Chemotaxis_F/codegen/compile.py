

import subprocess
import ctypes

def compile(name, src=None):

    if src:
        with open(name + '.f90', 'w') as f:
            f.write(src)

    subprocess.check_call('gfortran -pg -g -fbounds-check -fbacktrace -shared -fPIC {name}.f90 -o {name}.so'.format(
            name=name), shell=True)

    so = ctypes.CDLL('./{name}.so'.format(name=name))

    return so
    
    

