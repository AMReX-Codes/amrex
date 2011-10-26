
from pyboxlib import *

fboxlib.open()

la = layout()
la.create(boxes=[ [(1,1), (3,3)],
                  [(7,7), (100,100)],
                  ])

lmfab = lmultifab()
lmfab.create(la)

b = lmfab.fab(2)                        # get the second box of lmfab
b[2:8,2:8] = True                       # direct access (local indexing) to the fortran data array!

print b.shape
print 'bx: ', b.bx                      # lo/hi bounds of cell indexes (boxes)
print 'ibx:', b.ibx                     # lo/hi bounds of indexes
print 'pbx:', b.pbx                     # lo/hi of physical indexes

#lmfab.echo()

la = layout()
la.from_regrid(lmfab)                   # regrid based on tagged cells in lmfab

la.echo()

mfab = multifab()
mfab.create(la)
#mfab.echo()

fboxlib.close()

