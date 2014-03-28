
import boxlib

def test_multifab():

  ndim = 3
  size = 64
  dof  = 1

  bx = boxlib.Box(lo=[1]*ndim, hi=[size]*ndim)
  ba = boxlib.BoxArray(boxes=[bx])
  ba.maxSize(32)

  mf = boxlib.MultiFab(ba, ncomp=dof, nghost=2)

  mf.FillBoundary(0, mf.nComp())


if __name__ == '__main__':
  test_multifab()
