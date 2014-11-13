
import fboxlib

from fboxlib.pybl import bl

from ctypes import byref, c_void_p, c_int, c_double, CFUNCTYPE, POINTER

# see tag_boxes.f90
TAGBOXES = CFUNCTYPE(None, c_void_p, c_void_p, c_double, c_int)

def regrid(layouts, multifabs, dxs, tag_boxes, max_levs=8, amr_buf_width=2, max_grid_size=64):

    # type(ml_layout), intent(inout) :: mla
    # type(multifab) , intent(inout) :: phi(:)
    # integer        , intent(inout) :: nlevs, max_levs
    # real(dp_t)     , intent(in   ) :: dx(:)
    # type(bc_tower) , intent(inout) :: the_bc_tower
    # integer        , intent(in   ) :: amr_buf_width, max_grid_size

    nlevs = len(layouts)

    lacptrs = (max_levs*c_void_p)()
    for i, l in enumerate(layouts):
      lacptrs[i] = l.cptr

    mfcptrs = (max_levs*c_void_p)()
    for i, l in enumerate(multifabs):
      mfcptrs[i] = l.cptr

    nlevs = c_int(nlevs)
    print "regridding..."
    bl.pybl_regrid(lacptrs, mfcptrs, byref(nlevs), c_int(max_levs), TAGBOXES(tag_boxes))
    print "regridding... done."

    mfs = []
    for i in range(nlevs.value):
        mf = fboxlib.multifab(mfcptrs[i])
        mf.get_info()
        mfs.append(mf)

    return mfs
