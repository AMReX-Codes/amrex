
from pyboxlib import *

pybl.open()

N = 8
n = N / pybl.mpi_size()

boxes = []
for k in range(pybl.mpi_size()):
    boxes.append( ((1,k*n+1), (N,(k+1)*n)) )

la = layout()
la.create(boxes=boxes)

mfab = multifab()
mfab.create(la, components=1, ghost_cells=2)

if pybl.mpi_rank() == 1:
    for n in la.local_boxes:
        fab = mfab.fab(n)
        fab[2,5] = 22.0

mfab.fill_boundary()
mfab.echo()
mfab.write('pltTEST', 'test')

pybl.close()
