
from pyboxlib import *

pybl.open()

mfab = multifab()
mfab.read('pltTEST', 'test')

if mfab.nboxes > 1:
    fab = mfab.fab(1)
    print fab[0,3]

pybl.close()
