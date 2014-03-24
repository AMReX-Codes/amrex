module cc_applyop_module

  use bl_constants_module
  use mg_module
  use ml_layout_module
  use bndry_reg_module

  implicit none

contains

  subroutine ml_cc_applyop(mla, mgt, res, full_soln, ref_ratio)

    use bl_prof_module
    use ml_restriction_module , only : ml_restriction
    use ml_prolongation_module, only : ml_interp_bcs
    use cc_ml_resid_module    , only : crse_fine_residual_cc

    type(ml_layout), intent(in)    :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: res(:)
    type( multifab), intent(inout) :: full_soln(:)
    integer        , intent(in   ) :: ref_ratio(:,:)

    integer :: nlevs
    type(multifab), allocatable  ::      soln(:)
    type(multifab), allocatable  ::        uu(:)
    type(multifab), allocatable  ::   uu_hold(:)
    type(multifab), allocatable  ::        rh(:) ! this will be set to zero
    type(multifab), allocatable  ::  temp_res(:)

    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box) :: pd, pdc
    type(layout) :: la, lac
    integer :: i, n, ng_fill, dm
    integer :: mglev, mglev_crse

    type(bl_prof_timer), save :: bpt
    integer                   :: lo(get_dim(res(1))),hi(get_dim(res(1))),ng
    real(kind=dp_t),  pointer :: resp(:,:,:,:)

    call build(bpt, "ml_cc_applyop")

    nlevs = mla%nlevel

    allocate(soln(nlevs), uu(nlevs), rh(nlevs), temp_res(nlevs))
    allocate(uu_hold(2:nlevs-1))

    allocate(brs_flx(2:nlevs))
    allocate(brs_bcs(2:nlevs))

    do n = 2,nlevs-1
       la = mla%la(n)
       call build(uu_hold(n),la,1,1)
       call setval( uu_hold(n), ZERO,all=.true.)
    end do

    do n = nlevs, 1, -1

       la = mla%la(n)
       call build(    soln(n), la, 1, 1)
       call build(      uu(n), la, 1, 1)
       call build(      rh(n), la, 1, 0)
       call build(temp_res(n), la, 1, 0)
       call setval(    soln(n), ZERO,all=.true.)
       call setval(      uu(n), ZERO,all=.true.)
       call setval(      rh(n), ZERO,all=.true.)
       call setval(temp_res(n), ZERO,all=.true.)

       ! zero residual just to be safe
       call setval(     res(n), ZERO,all=.true.)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2, other = .false.)
       call bndry_reg_rr_build(brs_flx(n), la, lac, ref_ratio(n-1,:), pdc, width = 0)

    end do

    dm = get_dim(rh(1))

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in 
    !  its ghost cells before we evaluate the initial residual  
    do n = 2,nlevs
       ng_fill = nghost(full_soln(n))
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, &
                             ref_ratio(n-1,:), ng_fill, -i)
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, &
                             ref_ratio(n-1,:), ng_fill, +i)
       end do
       call multifab_fill_boundary(full_soln(n))
    end do

    !   Make sure all periodic and internal boundaries are filled as well
    do n = 1,nlevs   
       call multifab_fill_boundary(full_soln(n))
    end do


    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n), mgt(n)%mm(mglev), &
                           mgt(n)%stencil_type, mgt(n)%lcross)
    end do

    ! Make sure to correct the coarse cells immediately next to fine grids
    !   using the averaged fine grid fluxes
    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc, &
                                  ref_ratio(n-1,:))

       call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))
    enddo

    ! still need to multiply residual by -1 to get (alpha - del dot beta grad)
    do n=1,nlevs
       call multifab_mult_mult_s_c(res(n),1,-1.d0,1,0)
    end do

    do n = 2,nlevs-1
       call destroy(uu_hold(n))
    end do

    do n = nlevs, 1, -1
       call destroy(soln(n))
       call destroy(uu(n))
       call destroy(rh(n))
       call destroy(temp_res(n))
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_bcs(n))
       call bndry_reg_destroy(brs_flx(n))
    end do

    call destroy(bpt)

 end subroutine ml_cc_applyop
!
! ******************************************************************************************
!
 subroutine ml_cc_n_applyop(mla, mgt, res, full_soln, ref_ratio)

    use bl_prof_module
    use ml_restriction_module , only : ml_restriction
    use ml_prolongation_module, only : ml_interp_bcs
    use cc_ml_resid_module    , only : crse_fine_residual_n_cc

    type(ml_layout), intent(in)    :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: res(:)
    type( multifab), intent(inout) :: full_soln(:)
    integer        , intent(in   ) :: ref_ratio(:,:)

    integer :: nlevs
    type(multifab), allocatable  ::      soln(:)
    type(multifab), allocatable  ::        uu(:)
    type(multifab), allocatable  ::   uu_hold(:)
    type(multifab), allocatable  ::        rh(:) ! this will be set to zero
    type(multifab), allocatable  ::  temp_res(:)

    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box) :: pd, pdc
    type(layout) :: la, lac
    integer :: i, n, ng_fill, dm, nComp
    integer :: mglev, mglev_crse

    type(bl_prof_timer), save :: bpt
    integer                   :: lo(get_dim(res(1))),hi(get_dim(res(1))),ng
    real(kind=dp_t),  pointer :: resp(:,:,:,:)

    call build(bpt, "ml_cc_applyop")

    nlevs = mla%nlevel

    allocate(soln(nlevs), uu(nlevs), rh(nlevs), temp_res(nlevs))
    allocate(uu_hold(2:nlevs-1))

    allocate(brs_flx(2:nlevs))
    allocate(brs_bcs(2:nlevs))

    do n = 2,nlevs-1
       la = mla%la(n)
       call build(uu_hold(n),la,1,1)
       call setval( uu_hold(n), ZERO,all=.true.)
    end do

    dm    = 2
    nComp = 2

    do n = nlevs, 1, -1

       la = mla%la(n)
       call build(    soln(n), la, 1, 1)
       call build(      uu(n), la, 1, 1)
       call build(      rh(n), la, 1, 0)
       call build(temp_res(n), la, 1, 0)
       call setval(    soln(n), ZERO,all=.true.)
       call setval(      uu(n), ZERO,all=.true.)
       call setval(      rh(n), ZERO,all=.true.)
       call setval(temp_res(n), ZERO,all=.true.)

       ! zero residual just to be safe
       call setval(     res(n), ZERO,all=.true.)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2, other = .false.)
       call bndry_reg_rr_build(brs_flx(n), la, lac, ref_ratio(n-1,:), pdc, nc = nComp, width = 0)

    end do

    dm = get_dim(rh(1))

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in 
    !  its ghost cells before we evaluate the initial residual  
    do n = 2,nlevs
       ng_fill = nghost(full_soln(n))
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, &
                             ref_ratio(n-1,:), ng_fill, -i)
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, &
                             ref_ratio(n-1,:), ng_fill, +i)
       end do
       call multifab_fill_boundary(full_soln(n))
    end do

    !   Make sure all periodic and internal boundaries are filled as well
    do n = 1,nlevs   
       call multifab_fill_boundary(full_soln(n))
    end do


    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n), mgt(n)%mm(mglev), &
                           mgt(n)%stencil_type, mgt(n)%lcross) 
    end do

    ! Make sure to correct the coarse cells immediately next to fine grids
    !   using the averaged fine grid fluxes
    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_n_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc, &
                                    ref_ratio(n-1,:))
       call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))
    enddo

    ! still need to multiply residual by -1 to get (alpha - del dot beta grad)
    do n=1,nlevs
       call multifab_mult_mult_s_c(res(n),1,-1.d0,1,0)
    end do

    do n = 2,nlevs-1
       call destroy(uu_hold(n))
    end do

    do n = nlevs, 1, -1
       call destroy(soln(n))
       call destroy(uu(n))
       call destroy(rh(n))
       call destroy(temp_res(n))
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_bcs(n))
       call bndry_reg_destroy(brs_flx(n))
    end do

    call destroy(bpt)

  end subroutine ml_cc_n_applyop

end module cc_applyop_module
