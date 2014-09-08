module ml_restrict_fill_module

  use define_bc_module
  use ml_cc_restriction_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  
  implicit none

contains

  ! This subroutine restricts cell-centered multi-level multifab data, 
  ! and fill all boundaries (ghost cells covered the same level, 
  ! coarse-fine, periodic and physical domain boundaries).
  ! 
  ! This subroutine works on single level too.  It will skip restriction and perform
  ! boundary filling.

  subroutine ml_restrict_and_fill(nlevs, mf, rr, bc, icomp, bcomp, nc, ng, &
       time_in, dx_in, prob_lo_in, prob_hi_in, stencil_width_input,fourth_order_input, &
       same_boundary)

    integer, intent(in) :: nlevs
    integer, intent(in) :: rr(:,:)
    type(bc_level), intent(in) :: bc(:)
    type(multifab), intent(inout) :: mf(:)
    integer, intent(in), optional :: icomp, bcomp, nc, ng, stencil_width_input
    logical, intent(in), optional :: fourth_order_input, same_boundary
    real(kind=dp_t), intent(in), optional :: time_in, dx_in(:), prob_lo_in(:), prob_hi_in(:)
    
    integer :: licomp, lbcomp, lnc, lng
    logical :: lsameb
    integer :: n, i
    
    licomp = 1;              if (present(icomp)) licomp = icomp
    lbcomp = 1;              if (present(bcomp)) lbcomp = bcomp
    lnc    = ncomp (mf(1));  if (present(nc   )) lnc    = nc
    lng    = nghost(mf(1));  if (present(ng   )) lng    = ng
    lsameb = .false.;        if (present(same_boundary)) lsameb = same_boundary

    do n=nlevs,2,-1
       call ml_cc_restriction_c(mf(n-1), licomp, mf(n), licomp, rr(n-1,:), lnc)
    end do

    n = 1
    call multifab_fill_boundary_c(mf(n), licomp, lnc, lng)
    
    if (lsameb) then
       do i=1, lnc
          call multifab_physbc(mf(n), licomp+i-1, lbcomp, 1, bc(n), &
               time_in,dx_in,prob_lo_in,prob_hi_in)
       end do
    else
       call multifab_physbc(mf(n), licomp, lbcomp, lnc, bc(n), &
            time_in,dx_in,prob_lo_in,prob_hi_in)
    end if
    
    do n=2,nlevs
       if (lsameb) then
          do i=1, lnc
             call multifab_fill_ghost_cells(mf(n), mf(n-1), lng, rr(n-1,:), bc(n-1), bc(n), &
                                            licomp+i-1, lbcomp, 1, &
                                            stencil_width_input      = stencil_width_input, &
                                            fourth_order_input       = fourth_order_input, &
                                            fill_crse_boundary_input = .false., &
                                            fill_crse_physbc_input   = .false.)
          end do
       else
          call multifab_fill_ghost_cells(mf(n), mf(n-1), lng, rr(n-1,:), bc(n-1), bc(n), &
                                         licomp, lbcomp, lnc, &
                                         stencil_width_input      = stencil_width_input, &
                                         fourth_order_input       = fourth_order_input, &
                                         fill_crse_boundary_input = .false., &
                                         fill_crse_physbc_input   = .false.)
       end if
    end do

  end subroutine ml_restrict_and_fill

end module ml_restrict_fill_module
