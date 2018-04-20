import os, sys
import yt
import shutil
import numpy as np
import scipy.constants as scc

yt.funcs.mylog.setLevel(0)
os.system('./main2d.gnu.MPI.OMP.ex inputs2d')

filename = 'plt02000'
ds = yt.load( filename )
ad = ds.all_data()
ex = np.reshape(ad['boxlib', 'Ex'].v,(128,128))
ez = np.reshape(ad['boxlib', 'Ez'].v,(128,128))
by = np.reshape(ad['boxlib', 'By'].v,(128,128))
ener = np.sum(ex**2 + ez**2 + scc.c**2*by**2)

print('#################')
print('#### 2d test ####')
print('#################')
print('Energy = ' + str(ener*1.e-15))
print('Should be: ')
print('FILTER OFF: 2.6570911530826467e+17')
print('FILTER ON : 1.3963228479075573')

os.system('rm -r ./plt*')
