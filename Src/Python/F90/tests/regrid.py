
import fboxlib
fboxlib.open()

ba = fboxlib.boxarray([[(0,0), (100,100)]])
la = fboxlib.layout(ba)
mf = fboxlib.multifab(la, nc=3, ng=0)

print "#"*80
print "# before regridding"
la.echo()

fab = mf.fab(1)
fab.array[...] = 1.0
fab.array[10:20,10:20] = 2.0
fab.array[50:60,50:60] = 2.0

def tag_boxes(mf, tb, dx, lev):
    if lev > 1:
        return
    mf = fboxlib.multifab(cptr=mf)
    tb = fboxlib.lmultifab(cptr=tb)
    mfab = mf.fab(1)
    tfab = tb.fab(1)
    tfab.array[mfab.array[:,:,0] > 1.0] = 1

mfs = fboxlib.regrid([la], [mf], [0.5], tag_boxes)

print "#"*80
print "# after regridding"
for mf in mfs:
    mf.layout.echo()
