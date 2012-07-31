
import numpy as np
import matplotlib.pylab as plt


from pyboxlib import pybl, layout, multifab


nc = 2
ng = 2
nx = 32
dx = 2.0/nx



pybl.open()

la = layout()
la.create(boxes=[[ (0,0), (nx-1, nx-1) ]])

q = multifab()
q.create(la, components=nc, ghost_cells=ng)


for n in la.local_boxes:
    fab = q.fab(n)

    for j in fab.bxrange(2):
        y = -1.0 + j * dx

        for i in fab.bxrange(1):
            x = -1.0 + i * dx

            fab[i, j, 0] = np.sin(x)

q.fill_boundary()

fab = q.fab(1)
plt.imshow(fab.array[:, :, 0])
plt.show()

pybl.close()


