
subroutine fmain () bind(c)

  use baselib_module

  implicit none

  integer, parameter :: ncell=16

  integer :: lo(3), hi(3), ncomp, nghost
  type(Box)      :: domain, bx
  type(BoxArray) :: ba, ba2
  type(MultiFab) :: mf
  type(MFIter)   :: mfi
  double precision, pointer, dimension(:,:,:,:) :: p

  lo = 0
  hi = ncell-1
  domain = Box(lo,hi)

  call boxarray_build(ba, domain)

  ba2 = ba
  print *, ba2%owner, ba%owner

  ncomp = 1
  nghost = 0
  call multifab_build(mf, ba, ncomp, nghost)

  call mfiter_build(mfi, mf, tiling=.true.)
  do while (mfi%next())
     bx = mfi%tilebox()
     p=> mf%dataPtr(mfi)
     print *, "tile: " , bx%lo, bx%hi
     print *, "lbound, ubound of p, " ,lbound(p), ubound(p)
  end do

end subroutine fmain
