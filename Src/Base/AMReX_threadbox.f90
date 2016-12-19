module threadbox_module

  implicit none

  integer, save :: nthreads = 0
  integer, save :: nfacs
  integer, allocatable, save :: nthreads_facs(:)
  double precision, allocatable, save :: ntfinv(:)

  private

  public :: build_threadbox_2d, build_threadbox_3d, get_lo_hi

contains

  ! this subroutine should not be called inside OMP region
  subroutine build_threadbox_3d(nthreads_in, boxsize, blocksize_min, nb)
    integer, intent(in) :: nthreads_in, boxsize(3), blocksize_min
    integer, intent(out) :: nb(3)
    
    integer :: ifac
    double precision :: dncells(3), dntry(3)
    
    if (nthreads .ne. nthreads_in) then
       nthreads = nthreads_in
       call fa()
    end if

    dncells = boxsize
    do ifac=1,nfacs
       dntry = dncells*ntfinv(ifac)
       if (dntry(2) .lt. blocksize_min .and. dntry(3) .lt. blocksize_min) then
          ! split in direction 1 if we can
          if (dntry(1) .ge. blocksize_min) then
             dncells(1) = dntry(1)
          else
             exit
          end if
       else if (dncells(2) .gt. dncells(3)) then
          ! split in direction 2
          dncells(2) = dntry(2)
       else
          ! split in direction 3
          dncells(3) = dntry(3)
       end if
    end do

    nb = nint(boxsize/dncells)

    if (nb(1)*nb(2)*nb(3) .ne. nthreads) nb = 0

  end subroutine build_threadbox_3d

  ! this subroutine should not be called inside OMP region
  subroutine build_threadbox_2d(nthreads_in, boxsize, blocksize_min, nb)
    integer, intent(in) :: nthreads_in, boxsize(2), blocksize_min
    integer, intent(out) :: nb(2)
    
    integer :: ifac
    double precision :: dncells(2), dntry(2)
    
    if (nthreads .ne. nthreads_in) then
       nthreads = nthreads_in
       call fa()
    end if

    dncells = boxsize
    do ifac=1,nfacs
       dntry = dncells*ntfinv(ifac)
       if (dntry(2) .lt. blocksize_min) then
          ! split in direction 1 if we can
          if (dntry(1) .ge. blocksize_min) then
             dncells(1) = dntry(1)
          else
             exit
          end if
       else 
          ! split in direction 2
          dncells(2) = dntry(2)
       end if
    end do

    nb = nint(boxsize/dncells)

    if (nb(1)*nb(2) .ne. nthreads) nb = 0

  end subroutine build_threadbox_2d

  subroutine fa()
    integer :: facs(100), f, rem

    nfacs = 0
    f = 2
    rem = nthreads
    do while (rem .ne. 1)
       if (mod(rem, f) .eq. 0) then
          rem = rem/f
          nfacs = nfacs+1
          facs(nfacs) = f
       else
          f = f + 1
       end if
    end do
    
    allocate(nthreads_facs(nfacs))
    allocate(ntfinv(nfacs))
    nthreads_facs(1:nfacs) = facs(nfacs:1:-1)
    ntfinv = 1.d0/nthreads_facs
  end subroutine fa


  subroutine get_lo_hi(ncells, nblks, lo, hi)
    integer, intent(in) :: ncells, nblks
    integer, intent(out) :: lo(0:nblks-1), hi(0:nblks-1)
    integer ::  blksize, nleft, i
    blksize = ncells/nblks
    nleft = ncells - nblks * blksize  ! Note nleft < nblks
    do i=0,nblks-1
       lo(i) = i*blksize + min(nleft, i)
    end do
    if (nblks .ge. 2) then
       do i=0,nblks-2
          hi(i) = lo(i+1)-1
       end do
    end if
    hi(nblks-1) = ncells-1
  end subroutine get_lo_hi

end module threadbox_module
