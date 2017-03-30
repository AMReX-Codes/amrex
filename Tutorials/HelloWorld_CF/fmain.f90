
subroutine amrex_fmain () bind(c)

  use amrex_base_module

  implicit none

  integer, parameter :: ncell=64, sz=32

  integer :: lo(3), hi(3), i, j, k
  real(amrex_real)      :: rmin, rmax, rnm0, rnm1, rnm2
  type(amrex_box)       :: domain, bx
  type(amrex_boxarray)  :: ba
  type(amrex_distromap) :: dm
  type(amrex_multifab)  :: mf
  type(amrex_mfiter)    :: mfi
  real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: p

  ! define the lower and upper corner of a 3D domain
  lo = 0
  hi = ncell-1

  ! build a box for the domain
  domain = amrex_box(lo,hi)

  ! build a box array from the domain box
  call amrex_boxarray_build(ba, domain)
  ! break the box array into smaller boxes
  call ba%maxSize(sz)

  call amrex_distromap_build(dm, ba)

  ! build a multifab on the box array with 1 component, 0 ghost cells
  call amrex_multifab_build(mf, ba, dm, 1, 0)

  call amrex_distromap_destroy(dm)
  call amrex_boxarray_destroy(ba)

  !$OMP PARALLEL PRIVATE(mfi,bx,p,i,j,k)
  !
  ! build a multifab iterator with tiling
  call amrex_mfiter_build(mfi, mf, tiling=.true.)
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

  call amrex_mfiter_destroy(mfi)
  !$OMP END PARALLEL

  rmin = mf%min(1)    ! min of component 1
  rmax = mf%max(1)    ! max of component 1
  rnm0 = mf%norm0()
  rnm1 = mf%norm1()
  rnm2 = mf%norm2()
  
  call amrex_multifab_destroy(mf)

  if (amrex_parallel_ioprocessor()) then
     print *, "min      = ", rmin
     print *, "max      = ", rmax
     print *, "max norm = ", rnm0
     print *, "L1  norm = ", rnm1
     print *, "L2  norm = ", rnm2
  end if

end subroutine amrex_fmain
