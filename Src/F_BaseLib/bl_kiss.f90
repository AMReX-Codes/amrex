module  bl_kiss_module

  ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
  ! (2) A 3-shift shift-register generator, period 2^32-1,
  ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497 > 2^59
  !  Overall period > 2^123;  Default seeds x,y,z,w.
  !  Set your own seeds with statement i=kisset(ix,iy,iz,iw).

  implicit none

  private
  public :: kiss, kiss_init, kiss_rn

  type :: kiss_rn
     private
     integer :: x = 123456789
     integer :: y = 362436069
     integer :: z = 521288629
     integer :: w = 916191069
  end type kiss_rn
  type(kiss_rn), save :: the_kiss

  interface kiss
     module procedure kiss_d
     module procedure kiss_k
  end interface

  interface kiss_init
     module procedure kiss_init_d
     module procedure kiss_init_k
  end interface

  interface kiss_fill
     module procedure kiss_fill_0d
     module procedure kiss_fill_1d
     module procedure kiss_fill_2d
     module procedure kiss_fill_3d
     module procedure kiss_fill_4d
     module procedure kiss_fill_5d
     module procedure kiss_fill_6d
     module procedure kiss_fill_7d
  end interface

contains

  function kiss_k (k, range) result(r)
    type(kiss_rn), intent(inout) :: k
    integer, intent(in), optional :: range
    integer :: r
    k%x = 69069 * k%x + 1327217885
    k%y = m (m (m (k%y, 13), - 17), 5)
    k%z = 18000 * iand (k%z, 65535) + ishft (k%z, - 16)
    k%w = 30903 * iand (k%w, 65535) + ishft (k%w, - 16)
    r = k%x + k%y + ishft (k%z, 16) + k%w
    if ( present(range) ) r = mod(r,range)
  contains
    function m(k, n)
      integer :: m
      integer, intent(in) :: k, n
      m = ieor (k, ishft (k, n) )
    end function m
  END FUNCTION kiss_k

  function kiss_d (range) result(r)
    integer :: r
    integer, intent(in), optional :: range
    r = kiss_k(the_kiss, range)
  end function kiss_d

  subroutine kiss_init_k (k, ix, iy, iz, iw)
    type(kiss_rn), intent(out) :: k
    integer, intent(in) :: ix, iy, iz, iw
    k%x = ix
    k%y = iy
    k%z = iz
    k%w = iw
  end subroutine kiss_init_k

  subroutine kiss_init_d (ix, iy, iz, iw)
    integer, intent(in) :: ix, iy, iz, iw
    call kiss_init_k(the_kiss, ix, iy, iz, iw)
  end subroutine kiss_init_d

  subroutine kiss_fill_0d(a, range)
    integer, intent(out) :: a
    integer, intent(in), optional :: range
    a = kiss(range)
  end subroutine kiss_fill_0d

  subroutine kiss_fill_1d(a, range)
    integer, intent(out)  :: a(:)
    integer, intent(in), optional :: range
    integer :: i
    do i = 1, size(a,1)
       a(i) = kiss(range)
    end do

  end subroutine kiss_fill_1d

  subroutine kiss_fill_2d(a, range)
    integer, intent(out)  :: a(:,:)
    integer, intent(in), optional :: range
    integer :: i, j
    do j = 1, size(a,2)
       do i = 1, size(a,1)
          a(i,j) = kiss(range)
       end do
    end do

  end subroutine kiss_fill_2d

  subroutine kiss_fill_3d(a, range)
    integer, intent(out)  :: a(:,:,:)
    integer, intent(in), optional :: range
    integer :: i, j, k
    do k = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             a(i,j,k) = kiss(range)
          end do
       end do
    end do

  end subroutine kiss_fill_3d

  subroutine kiss_fill_4d(a, range)
    integer, intent(out)  :: a(:,:,:,:)
    integer, intent(in), optional :: range
    integer :: i, j, k, l

    do l = 1, size(a,4)
       do k = 1, size(a,3)
          do j = 1, size(a,2)
             do i = 1, size(a,1)
                a(i,j,k,l) = kiss(range)
             end do
          end do
       end do
    end do

  end subroutine kiss_fill_4d

  subroutine kiss_fill_5d(a, range)
    integer, intent(out)  :: a(:,:,:,:,:)
    integer, intent(in), optional :: range
    integer :: i, j, k, l, m

    do m = 1, size(a,5)
       do l = 1, size(a,4)
          do k = 1, size(a,3)
             do j = 1, size(a,2)
                do i = 1, size(a,1)
                   a(i,j,k,l,m) = kiss(range)
                end do
             end do
          end do
       end do
    end do

  end subroutine kiss_fill_5d

  subroutine kiss_fill_6d(a, range)
    integer, intent(out)  :: a(:,:,:,:,:,:)
    integer, intent(in), optional :: range
    integer :: i, j, k, l, m, n

    do n = 1, size(a,6)
       do m = 1, size(a,5)
          do l = 1, size(a,4)
             do k = 1, size(a,3)
                do j = 1, size(a,2)
                   do i = 1, size(a,1)
                      a(i,j,k,l,m,n) = kiss(range)
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine kiss_fill_6d

  subroutine kiss_fill_7d(a, range)
    integer, intent(out)  :: a(:,:,:,:,:,:,:)
    integer, intent(in), optional :: range
    integer :: i, j, k, l, m, n, o

    do o = 1, size(a,7)
       do n = 1, size(a,6)
          do m = 1, size(a,5)
             do l = 1, size(a,4)
                do k = 1, size(a,3)
                   do j = 1, size(a,2)
                      do i = 1, size(a,1)
                         a(i,j,k,l,m,n,o) = kiss(range)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine kiss_fill_7d

end module bl_kiss_module
