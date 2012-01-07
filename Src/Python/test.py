
from pyboxlib import *

fboxlib.open()

la = layout()
la.create(boxes=[ [(1,1), (3,3)],
                  [(4,1), (6,3)],
                  ])

mfab1 = multifab()
mfab1.create(la)

fab = mfab1.fab(1)
fab[0,0] = 22.0

fab = mfab1.fab(2)
fab[-1,-1] = 44.0

mfab2 = multifab()
mfab2.create_from_bbox(mfab1)

mfab1.copy(mfab2)

mfab1.echo()
mfab2.echo()

fboxlib.close()
