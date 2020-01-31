! Copyright 2019 Maxence Thevenet, Remi Lehe
!
! This file is part of WarpX.
!
! License: BSD-3-Clause-LBNL

module warpx_ES_utils

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  implicit none

contains

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
