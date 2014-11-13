
import fboxlib.fcboxlib as fcboxlib
import numpy as np

from fboxlib.multifab import multifab

def regrid(layouts, multifabs, dxs, tag_boxes, max_levs=8, amr_buf_width=2, max_grid_size=64):

    # XXX: passing layouts here is a bit redundant...

    nlevs = len(layouts)
    lacptrs = np.zeros((max_levs,), np.long)
    mfcptrs = np.zeros((max_levs,), np.long)
    for i in range(nlevs):
        lacptrs[i] = layouts[i].cptr
        mfcptrs[i] = multifabs[i].cptr

    nlevs = fcboxlib.regrid(lacptrs, mfcptrs, nlevs, max_levs, tag_boxes)

    return [ multifab(cptr=x) for x in mfcptrs[:nlevs] ]
