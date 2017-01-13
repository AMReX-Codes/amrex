#include <WarpXUtil.H>

void shiftMF(MultiFab& mf, const Geometry& geom, int num_shift, 
	     int dir, const IntVect& nodalflag) {

  // create tmp copy with num_shift ghost cells
  BoxArray ba = mf.boxArray();
  MultiFab tmpmf(ba, mf.nComp(), num_shift, Fab_allocate, nodalflag);
  tmpmf.setVal(0.0);
  MultiFab::Copy(tmpmf, mf, 0, 0, 1, 0);
  tmpmf.FillBoundary(geom.periodicity());
  
  // copy from tmpmf to mf using the shifted boxes
  for (MFIter mfi(mf); mfi.isValid(); ++mfi ) {
    const Box& dstBox = mfi.validbox();
    Box srcBox(dstBox.smallEnd(), dstBox.bigEnd());
    srcBox.shift(dir, num_shift);
    mf[mfi].copy(tmpmf[mfi], srcBox, 0, dstBox, 0, 1);
  }

  mf.FillBoundary(geom.periodicity());
}
