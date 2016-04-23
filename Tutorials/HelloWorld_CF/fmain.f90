
subroutine fmain () bind(c)

  use baselib_module

  implicit none

  integer, parameter :: ncell=16

  integer :: lo(3), hi(3), ncomp, nghost, i
  type(Box)      :: domain, bx
  type(BoxArray) :: ba
  type(MultiFab) :: mf
  type(MFIter)   :: mfi
  double precision, pointer, dimension(:,:,:,:) :: p

  lo = 0
  hi = ncell-1
  domain = Box(lo,hi)

  ba = BoxArray(domain)

  ncomp = 1
  nghost = 0
  mf = MultiFab(ba, ncomp, nghost)

  mfi = MFIter(mf, tiling=.true.)

  do while (mfi%next())
     bx = mfi%tilebox()
     p=> mf%dataPtr(mfi)
     print *, "tile: " , bx%lo, bx%hi
     print *, "lbound, ubound of p, " ,lbound(p), ubound(p)
  end do

end subroutine fmain
