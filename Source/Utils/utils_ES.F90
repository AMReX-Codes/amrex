module warpx_ES_utils

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  implicit none

contains

  subroutine warpx_sum_fine_to_crse_nodal_3d (lo, hi, lrat, crse, clo, chi, fine, flo, fhi) &
       bind(c, name="warpx_sum_fine_to_crse_nodal_3d")

    integer, intent(in)             ::   lo(3),  hi(3)
    integer, intent(in)             ::  clo(3), chi(3)
    integer, intent(in)             ::  flo(3), fhi(3)
    integer, intent(in)             ::  lrat(3)
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(amrex_real), intent(in)    :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

    integer :: i, j, k, ii, jj, kk

    do k        = lo(3), hi(3)
       kk       = k * lrat(3)
       do j     = lo(2), hi(2)
          jj    = j * lrat(2)
          do i  = lo(1), hi(1)
             ii = i * lrat(1)
             crse(i,j,k)  =  fine(ii,jj,kk)                              + &
! These six fine nodes are shared by two coarse nodes...
                  0.5d0   * (fine(ii-1,jj,kk)     + fine(ii+1,jj,kk)     + &
                             fine(ii,jj-1,kk)     + fine(ii,jj+1,kk)     + &
                             fine(ii,jj,kk-1)     + fine(ii,jj,kk+1))    + &
! ... these twelve are shared by four...
                  0.25d0  * (fine(ii,jj-1,kk-1)   + fine(ii,jj+1,kk-1)   + &
                             fine(ii,jj-1,kk+1)   + fine(ii,jj+1,kk+1)   + &
                             fine(ii-1,jj,kk-1)   + fine(ii+1,jj,kk-1)   + &
                             fine(ii-1,jj,kk+1)   + fine(ii+1,jj,kk+1)   + &
                             fine(ii-1,jj-1,kk)   + fine(ii+1,jj-1,kk)   + &
                             fine(ii-1,jj+1,kk)   + fine(ii+1,jj+1,kk))  + &
! ... and these eight are shared by eight
                  0.125d0 * (fine(ii-1,jj-1,kk-1) + fine(ii-1,jj-1,kk+1) + &
                             fine(ii-1,jj+1,kk-1) + fine(ii-1,jj+1,kk+1) + &
                             fine(ii+1,jj-1,kk-1) + fine(ii+1,jj-1,kk+1) + &
                             fine(ii+1,jj+1,kk-1) + fine(ii+1,jj+1,kk+1))
! ... note that we have 27 nodes in total...
             crse(i,j,k) = crse(i,j,k) / 8.d0
          end do
       end do
    end do

  end subroutine warpx_sum_fine_to_crse_nodal_3d

  subroutine warpx_sum_fine_to_crse_nodal_2d (lo, hi, lrat, crse, clo, chi, fine, flo, fhi) &
       bind(c, name="warpx_sum_fine_to_crse_nodal_2d")

    integer, intent(in)             ::   lo(2),  hi(2)
    integer, intent(in)             ::  clo(2), chi(2)
    integer, intent(in)             ::  flo(2), fhi(2)
    integer, intent(in)             ::  lrat(2)
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2))
    real(amrex_real), intent(in)    :: fine(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i, j, ii, jj

    do j     = lo(2), hi(2)
       jj    = j * lrat(2)
       do i  = lo(1), hi(1)
          ii = i * lrat(1)
          crse(i,j)  =  fine(ii,jj)                             + &
! These four fine nodes are shared by two coarse nodes...
               0.5d0   * (fine(ii-1,jj)     + fine(ii+1,jj)     + &
               fine(ii,jj-1)     + fine(ii,jj+1))               + &
! ... and these four are shared by four...
               0.25d0  * (fine(ii-1,jj-1)   + fine(ii-1,jj+1)   + &
               fine(ii-1,jj+1)   + fine(ii+1,jj+1))
! ... note that we have 9 nodes in total...
             crse(i,j) = crse(i,j) / 4.d0
       end do
    end do

  end subroutine warpx_sum_fine_to_crse_nodal_2d

  subroutine warpx_zero_out_bndry_3d (lo, hi, input_data, bndry_data, mask) &
       bind(c,name='warpx_zero_out_bndry_3d')

    integer(c_int),   intent(in   ) :: lo(3), hi(3)
    double precision, intent(inout) :: input_data(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision, intent(inout) :: bndry_data(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
    integer(c_int),   intent(in   ) :: mask (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (mask(i,j,k) .eq. 1) then
                bndry_data(i,j,k) = input_data(i,j,k)
                input_data(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine warpx_zero_out_bndry_3d

  subroutine warpx_zero_out_bndry_2d (lo, hi, input_data, bndry_data, mask) &
       bind(c,name='warpx_zero_out_bndry_2d')

    integer(c_int),   intent(in   ) :: lo(2), hi(2)
    double precision, intent(inout) :: input_data(lo(1):hi(1),lo(2):hi(2))
    double precision, intent(inout) :: bndry_data(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)
    integer(c_int),   intent(in   ) :: mask (lo(1):hi(1),lo(2):hi(2))

    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (mask(i,j) .eq. 1) then
             bndry_data(i,j) = input_data(i,j)
             input_data(i,j) = 0.d0
          end if
       end do
    end do

  end subroutine warpx_zero_out_bndry_2d

  subroutine warpx_build_mask_3d (lo, hi, tmp_mask, mask, ncells) &
       bind(c,name='warpx_build_mask_3d')
    integer(c_int),   intent(in   ) :: lo(3), hi(3)
    integer(c_int),   intent(in   ) :: ncells
    integer(c_int),   intent(in   ) :: tmp_mask(lo(1)-ncells:hi(1)+ncells,lo(2)-ncells:hi(2)+ncells,lo(3)-ncells:hi(3)+ncells)
    integer(c_int),   intent(inout) :: mask (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    integer :: i, j, k, ii, jj, kk, total

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             total = 0
             do ii = i-ncells, i+ncells
                do jj = j-ncells, j+ncells
                   do kk = k-ncells, k+ncells
                      total = total + tmp_mask(ii, jj, kk)
                   end do
                end do
             end do

             if (total .gt. 0) then
                mask(i,j,k) = 1
             else
                mask(i,j,k) = 0
             end if

          end do
       end do
    end do

  end subroutine warpx_build_mask_3d

  subroutine warpx_build_mask_2d (lo, hi, tmp_mask, mask, ncells) &
       bind(c,name='warpx_build_mask_2d')
    integer(c_int),   intent(in   ) :: lo(2), hi(2)
    integer(c_int),   intent(in   ) :: ncells
    integer(c_int),   intent(in   ) :: tmp_mask(lo(1)-ncells:hi(1)+ncells,lo(2)-ncells:hi(2)+ncells)
    integer(c_int),   intent(inout) :: mask (lo(1):hi(1),lo(2):hi(2))

    integer :: i, j, ii, jj, total

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          total = 0
          do ii = i-ncells, i+ncells
             do jj = j-ncells, j+ncells
                total = total + tmp_mask(ii, jj)
             end do
          end do

          if (total .gt. 0) then
             mask(i,j) = 1
          else
             mask(i,j) = 0
          end if

       end do
    end do

  end subroutine warpx_build_mask_2d

end module warpx_ES_utils
