subroutine SetIC(mf, lo, hi, dptr) bind(C, name="FSetIC")

  use amrex_error_module
  use, intrinsic :: iso_c_binding

  implicit none

  integer, intent(in) :: lo(3),hi(3)
  real(c_double), intent(inout) :: mf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
  real(c_double), intent(inout) :: dptr((hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(2)-lo(2)+1))
!  real(c_double), intent(inout) :: dptr(0:((hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(2)-lo(2)+1)-1))

  ! local variables
  integer i,j,k,n 
  integer(c_long) :: offset, tile_size(3)
  tile_size = hi - lo + 1
  n=1

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        offset = tile_size(1) * (j - lo(2)) + tile_size(1) * tile_size(2) * (k - lo(3)) + &
                      tile_size(1) * tile_size(2) * tile_size(3) * (n - 1) + 1 - lo(1)
        do i=lo(1),hi(1)
!           mf(i,j,k) = dptr(offset+1)
           ! Set initial conditions for the ODE for this cell.
           dptr(offset+i) = real(i+j+k, c_double)
        end do
     end do
  end do

end subroutine SetIC

subroutine SetIC_mfab(mf, lo, hi) bind(C, name="FSetIC_mfab")

  use amrex_error_module
  use, intrinsic :: iso_c_binding

  implicit none

  integer, intent(in) :: lo(3),hi(3)
  real(c_double), intent(inout) :: mf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))

  ! local variables
  integer i,j,k,n 
  integer(c_long) :: offset, tile_size(3)

  !gpu
  tile_size = hi - lo + 1
  n=1

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           ! Set initial conditions for the ODE for this cell.
           mf(i,j,k) = real(i+j+k, c_double)

        end do
     end do
  end do

end subroutine SetIC_mfab

subroutine SetSol(mf, lo, hi, dptr) bind(C, name="FSetSol")

  use amrex_error_module
  use, intrinsic :: iso_c_binding

  implicit none

  integer, intent(in) :: lo(3),hi(3)
  real(c_double), intent(inout) :: mf(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))
!  real(c_double), intent(inout) :: dptr(0:((hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(2)-lo(2)+1)-1))
  real(c_double), intent(inout) :: dptr((hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(2)-lo(2)+1))

  ! local variables
  integer i,j,k,n 
  integer(c_long) :: offset, tile_size(3)
  !gpu
  tile_size = hi - lo + 1
  n=1

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        offset = tile_size(1) * (j - lo(2)) + tile_size(1) * tile_size(2) * (k - lo(3)) + &
                      tile_size(1) * tile_size(2) * tile_size(3) * (n - 1) + 1 - lo(1)
        do i=lo(1),hi(1)

           ! Set solution for the ODE for this cell.
           mf(i,j,k) = dptr(offset+i)

        end do
     end do
  end do

end subroutine SetSol


subroutine fort_fab_copytoreal (lo, hi, bx_lo, bx_hi, dst, src, slo, shi, ncomp, src_comp) &
       bind(c, name='fort_fab_copytoreal')
    use iso_c_binding, only : c_long, c_double
!    use amrex_fort_module, only: amrex_real
    integer, intent(in) :: lo(3), hi(3), bx_lo(3), bx_hi(3), slo(3), shi(3)
    integer, intent(in), value :: ncomp
    integer, intent(in) :: src_comp(1:ncomp)
    real(c_double)             :: dst(1:((shi(1)-slo(1)+1)*(shi(2)-slo(2)+1)*(shi(2)-slo(2)+1)*ncomp))
!    real(c_double)             :: dst(*)
    real(c_double), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)

    integer :: i, j, k, n
    integer(c_long) :: offset, tile_size(3)

!!!    !gpu

    tile_size = bx_hi - bx_lo + 1
!    do       k = slo(3), shi(3)
!       do    j = slo(2), shi(2)
!          do i = slo(1), shi(1)
!             offset = tile_size(1) * (i - bx_lo(1)) + tile_size(1) * tile_size(2) * (j - bx_lo(2)) + &
!                      tile_size(1) * tile_size(2) * tile_size(3) * (k - bx_lo(3)) + i - 1
!             do n = 1, ncomp
!                print*, "real location:", offset+n, "ijkn:",i,j,k,n
!                dst(i,j,k,n)  = src(offset+i)
!!                if(src_comp(n) .gt. 0) then
!!                   dst(offset+n) = src(i,j,k,src_comp(n))
!!                endif
!             end do
!          end do
!       end do
!    end do
    do n = 1, ncomp
       do       k = slo(3), shi(3)
          do    j = slo(2), shi(2)
             offset = tile_size(1) * (j - bx_lo(2)) + tile_size(1) * tile_size(2) * (k - bx_lo(3)) + &
                      tile_size(1) * tile_size(2) * tile_size(3) * (n - 1) + 1 - bx_lo(1)
          do i = slo(1), shi(1)
!                 print*, "real location:", offset+i, "ijkn:",i,j,k,n
                dst(offset+i)  = src(i,j,k,n)
!                if(src_comp(n) .gt. 0) then
!                   dst(offset+n) = src(i,j,k,src_comp(n))
!                endif
             end do
          end do
       end do
    end do

end subroutine fort_fab_copytoreal


subroutine fort_fab_copyfromreal (lo, hi, bx_lo, bx_hi, dst, dlo, dhi, ncomp, src) &
       bind(c, name='fort_fab_copyfromreal')
    use iso_c_binding, only : c_long, c_double
!    use amrex_fort_module, only: amrex_real
    integer, intent(in) :: lo(3), hi(3), bx_lo(3), bx_hi(3), dlo(3), dhi(3)
    integer, intent(in), value :: ncomp
    real(c_double), intent(in   ) :: src(1:((dhi(1)-dlo(1)+1)*(dhi(2)-dlo(2)+1)*(dhi(2)-dlo(2)+1)*ncomp))
!    real(c_double), intent(in   ) :: src(*)
    real(c_double), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i, j, k, n
    integer(c_long) :: offset, tile_size(3)

!!!    !gpu

    tile_size = bx_hi - bx_lo + 1

    do n = 1, ncomp
       do       k = dlo(3), dhi(3)
          do    j = dlo(2), dhi(2)
             offset = tile_size(1) * (j - bx_lo(2)) + tile_size(1) * tile_size(2) * (k - bx_lo(3)) + &
                      tile_size(1) * tile_size(2) * tile_size(3) * (n - 1) + 1 - bx_lo(1)
             do i = dlo(1), dhi(1)
                dst(i,j,k,n)  = src(offset+i)
             end do
          end do
       end do
    end do

end subroutine fort_fab_copyfromreal
