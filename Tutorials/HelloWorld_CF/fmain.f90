
subroutine fmain () bind(c)

  use boxlib_module

  implicit none

  integer, parameter :: ncell=64, sz=32

  integer :: lo(3), hi(3), i, j, k
  double precision :: rmin, rmax, rnm0, rnm1, rnm2
  type(Box)      :: domain, bx
  type(BoxArray) :: ba
  type(MultiFab) :: mf
  type(MFIter)   :: mfi
  double precision, pointer, dimension(:,:,:,:) :: p

  ! define the lower and upper corner of a 3D domain
  lo = 0
  hi = ncell-1

  ! build a box for the domain
  domain = Box(lo,hi)

  ! build a box array from the domain box
  call boxarray_build(ba, domain)
  ! break the box array into smaller boxes
  call ba%maxSize(sz)

  ! build a multifab on the box array with 1 component, 0 ghost cells
  call multifab_build(mf, ba, 1, 0)

  !$OMP PARALLEL PRIVATE(mfi,bx,p,i,j,k)
  !
  ! build a multifab iterator with tiling
  call mfiter_build(mfi, mf, tiling=.true.)
  !
  ! loop over tiles
  do while (mfi%next())
     bx = mfi%tilebox()
     p => mf%dataPtr(mfi)
     do                   k = bx%lo(3) , bx%hi(3)
        do                j = bx%lo(2) , bx%hi(2)
           do concurrent (i = bx%lo(1) : bx%hi(1))
              p(i,j,k,1) = exp(-dble(i*i+j*j+k*k)/512.d0)
           end do
        end do
     end do
  end do
  !$OMP END PARALLEL

  rmin = mf%min(1)    ! min of component 1
  rmax = mf%max(1)    ! max of component 1
  rnm0 = mf%norm0()
  rnm1 = mf%norm1()
  rnm2 = mf%norm2()
  
  if (parallel_ioprocessor()) then
     print *, "min      = ", rmin
     print *, "max      = ", rmax
     print *, "max norm = ", rnm0
     print *, "L1  norm = ", rnm1
     print *, "L2  norm = ", rnm2
  end if

end subroutine fmain
