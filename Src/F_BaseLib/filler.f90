module filler_module

  use plotfile_module

  implicit none

  interface blow_out_to_fab
     module procedure blow_out_to_fab_2d
     module procedure blow_out_to_fab_3d
  end interface

  interface blow_out_to_sub_fab
     module procedure blow_out_to_sub_fab_2d
     module procedure blow_out_to_sub_fab_3d
  end interface

contains

  ! wrongo
  ! Performs a piecewise constant interpolation of coarse data onto
  ! the fine cells of a fab.

  subroutine blow_out_to_fab_2d(f_fab, tlo, pf, comps, max_level)
    type(plotfile) :: pf
    integer, intent(in) :: tlo(:)
    integer, intent(in) :: comps(:)
    real(dp_t), intent(out) :: f_fab(tlo(1):,tlo(2):,:)
    integer, intent(in), optional :: max_level
    real(kind=dp_t), pointer :: cb(:,:,:,:)
    integer :: level
    integer :: i, rr(2), j, n
    integer :: ii, jj
    integer :: iii, jjj
    integer :: ir, jr
    integer, dimension(pf%dim) :: lo, hi
    integer :: nc

    nc = ubound(f_fab, dim = 3)

    level = plotfile_nlevels(pf)
    if ( present(max_level) ) level = min(max_level, level)
    do i = 1, level
       rr = 1
       do j = i, level-1
          rr = rr*plotfile_refrat_n(pf, j)
       end do
       do j = 1, plotfile_nboxes_n(pf, i)
          call fab_bind_comp_vec(pf, i, j, comps)
          cb => plotfile_dataptr(pf, i, j)
          do n = 1, 2
             lo(n) = lbound(cb, dim=n)
             hi(n) = ubound(cb, dim=n)
          end do
          do n = 1, nc
             do jj = lo(2), hi(2)
                do ii = lo(1), hi(1)
                   do jjj = 0, rr(2) - 1
                      jr = jj*rr(2) + jjj
                      do iii = 0, rr(1) - 1
                         ir = ii*rr(1) + iii
                         f_fab(ir, jr, n) = cb(ii, jj, 1, n)
                      end do
                   end do
                end do
             end do
          end do
          call fab_unbind(pf, i, j)
       end do
    end do
  end subroutine blow_out_to_fab_2d

  subroutine blow_out_to_fab_3d(f_fab, tlo, pf, comps, max_level)
    type(plotfile), intent(inout) :: pf
    integer, intent(in) :: tlo(:)
    integer, intent(in) :: comps(:)
    real(dp_t), intent(out) :: f_fab(tlo(1):,tlo(2):,tlo(3):,:)
    integer, intent(in), optional :: max_level
    integer :: level
    real(kind=dp_t), pointer :: cb(:,:,:,:)
    integer :: i, rr(3), j, n
    integer :: ii, jj, kk
    integer :: iii, jjj, kkk
    integer :: ir, jr, kr
    integer, dimension(pf%dim) :: lo, hi
    integer :: nc

    nc = ubound(f_fab, dim = 4)

    level = plotfile_nlevels(pf)
    if ( present(max_level) ) level = min(max_level, level)
    do i = 1, level
       rr = 1
       do j = i, level-1
          rr = rr*plotfile_refrat_n(pf, j)
       end do
       do j = 1, plotfile_nboxes_n(pf, i)
          call fab_bind_comp_vec(pf, i, j, comps)
          cb => plotfile_dataptr(pf, i, j)
          do n = 1, 3
             lo(n) = lbound(cb, dim=n)
             hi(n) = ubound(cb, dim=n)
          end do
          do n = 1, nc
             do kk = lo(3), hi(3)
                do jj = lo(2), hi(2)
                   do ii = lo(1), hi(1)
                      do kkk = 0, rr(3) - 1
                         kr = kk*rr(3) + kkk
                         do jjj = 0, rr(2) - 1
                            jr = jj*rr(2) + jjj
                            do iii = 0, rr(1) - 1
                               ir = ii*rr(1) + iii
                               f_fab(ir, jr, kr, n) = cb(ii, jj, kk, n)
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
          call fab_unbind(pf, i, j)
       end do
    end do
  end subroutine blow_out_to_fab_3d

  !! blows out to a subdomain of the fine fab, bounded by
  !! tlo and thi

  subroutine blow_out_to_sub_fab_2d(f_fab, tlo, thi, pf, comps, max_level)
    type(plotfile), intent(inout) :: pf
    ! tlo and thi are the indices specifying the bounds of the subdomain
    ! we are blowing out to.
    integer, intent(in) :: tlo(:), thi(:)
    integer, intent(in) :: comps(:)
    real(dp_t), intent(out) :: f_fab(tlo(1):,tlo(2):,:)
    integer, intent(in), optional :: max_level
    integer :: level
    real(kind=dp_t), pointer :: cb(:,:,:,:)
    integer :: i, rr(2), j, ii, jj, iii, jjj, n, ir, jr
    integer :: lo(2), hi(2)

    ! probably should put in some error checking here to make sure that
    ! we are not specifying a subdomain outside of the full domain

    level = plotfile_nlevels(pf)
    if ( present(max_level) ) level = min(max_level, level)

    do i = 1, level
       rr = 1
       do j = i, level-1
          rr = rr*plotfile_refrat_n(pf, j)
       end do

       do j = 1, plotfile_nboxes_n(pf, i)
          call fab_bind_comp_vec(pf, i, j, comps)
          cb => plotfile_dataptr(pf, i, j)

          ! get the limits of the current patch
          do n = 1, 2
             lo(n) = lbound(cb, dim=n)
             hi(n) = ubound(cb, dim=n)
          end do

          ! loop over the limits of the current patch -- in the future, we 
          ! should actually check here if the patch falls completely out of 
          ! our subdomain to speed things up.
          do jj = lo(2), hi(2)
             do ii = lo(1), hi(1)

                ! scale this patch up by rr to match the 
                ! resolution of f_fab
                do jjj = 0, rr(2) - 1
                   jr = jj*rr(2) + jjj
                   if (jr < tlo(2) .or. jr > thi(2)) cycle
                   
                   do iii = 0, rr(1) - 1
                      ir = ii*rr(1) + iii
                      if (ir < tlo(1) .or. ir > thi(1)) cycle
                      
                      f_fab(ir, jr, :) = cb(ii, jj, 1, :)
                   end do
                end do

             end do
          end do

          call fab_unbind(pf, i, j)
       end do

    end do
  end subroutine blow_out_to_sub_fab_2d

  subroutine blow_out_to_sub_fab_3d(f_fab, tlo, thi, pf, comps, max_level)
    type(plotfile), intent(inout) :: pf
    ! tlo and thi are the indices specifying the bounds of the subdomain
    ! we are blowing out to.
    integer, intent(in) :: tlo(:), thi(:)
    integer, intent(in) :: comps(:)
    real(dp_t), intent(out) :: f_fab(tlo(1):,tlo(2):,tlo(3):,:)
    integer, intent(in), optional :: max_level
    integer :: level
    real(kind=dp_t), pointer :: cb(:,:,:,:)
    integer :: i, rr(3), j, ii, jj, kk, iii, jjj, kkk, n, ir, jr, kr
    integer :: lo(3), hi(3)

    ! probably should put in some error checking here to make sure that
    ! we are not specifying a subdomain outside of the full domain

    level = plotfile_nlevels(pf)
    if ( present(max_level) ) level = min(max_level, level)
    do i = 1, level
       rr = 1
       do j = i, level-1
          rr = rr*plotfile_refrat_n(pf, j)
       end do
       do j = 1, plotfile_nboxes_n(pf, i)
          call fab_bind_comp_vec(pf, i, j, comps)
          cb => plotfile_dataptr(pf, i, j)

          ! get the limits of the current patch
          do n = 1, 3
             lo(n) = lbound(cb, dim=n)
             hi(n) = ubound(cb, dim=n)
          end do

          ! loop over the limits of the current patch -- in the future, we 
          ! should actually check here if the patch falls completely out of 
          ! our subdomain to speed things up.
          do kk = lo(3), hi(3)
             do jj = lo(2), hi(2)
                do ii = lo(1), hi(1)

                   ! scale this patch up by rr to match the 
                   ! resolution of f_fab
                   do kkk = 0, rr(3)-1
                      kr = kk*rr(3) + kkk
                      if ( kr < tlo(3) .or. kr > thi(3) ) cycle

                      do jjj = 0, rr(2) - 1
                         jr = jj*rr(2) + jjj
                         if (jr < tlo(2) .or. jr > thi(2)) cycle

                         do iii = 0, rr(1) - 1
                            ir = ii*rr(1) + iii
                            if (ir < tlo(1) .or. ir > thi(1)) cycle

                            f_fab(ir, jr, kr, :) = cb(ii, jj, kk, :)
                         end do
                      end do
                   end do

                end do
             end do
          end do

          call fab_unbind(pf, i, j)
       end do
    end do
  end subroutine blow_out_to_sub_fab_3d

end module filler_module
