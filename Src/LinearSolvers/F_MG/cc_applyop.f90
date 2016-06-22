module cc_applyop_module

  use bl_constants_module
  use mg_module
  use ml_layout_module
  use bndry_reg_module
  use define_bc_module

  implicit none

  private

  public :: cc_applyop, ml_cc_applyop

contains

  subroutine cc_applyop(mla,res,phi,alpha,beta,dx,the_bc_tower,bc_comp,base_level_in,crse_ratio_in,stencil_order_in)

    use mg_module             , only: mg_tower, mg_tower_build, mg_tower_destroy
    use cc_stencil_fill_module, only: stencil_fill_cc
    use stencil_types_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: res(:), phi(:)
    type(multifab) , intent(in   ) :: alpha(:), beta(:,:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: bc_comp
    integer, intent(in), optional  :: base_level_in
    integer, intent(in), optional  :: crse_ratio_in
    integer, intent(in), optional  :: stencil_order_in

    type(layout  ) :: la
    type(box     ) :: pd

    type(multifab) :: cell_coeffs
    type(multifab) :: edge_coeffs(mla%dim)

    type(mg_tower)  :: mgt(mla%nlevel)
    integer         :: d, dm, nlevs
    integer         :: stencil_order
    integer         :: base_level
    integer         :: crse_ratio

    ! MG solver defaults
    integer    :: n
    real(dp_t) ::  xa(mla%dim),  xb(mla%dim)
    real(dp_t) :: pxa(mla%dim), pxb(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "cc_applyop")

    dm = mla%dim
    nlevs = mla%nlevel

    stencil_order = 2
    if (present(stencil_order_in)) stencil_order = stencil_order_in

    base_level = 1
    if (present(base_level_in)) base_level = base_level_in

    crse_ratio = 2
    if (present(crse_ratio_in)) crse_ratio = crse_ratio_in

    do n = nlevs, 1, -1

       pd = layout_get_pd(mla%la(n))

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp),&
                           stencil_type = CC_CROSS_STENCIL, &
                           dh = dx(n,:), &
                           nodal = nodal_flags(res(nlevs)))

    end do

    !! Fill coefficient array

    do n = nlevs,1,-1

       la = mla%la(n)

       call multifab_build(cell_coeffs, la,          nc=1,ng=nghost(alpha(n)))
       call multifab_copy_c(cell_coeffs,1,alpha(n),1,nc=1,ng=nghost(alpha(n)))

       do d = 1,dm
          call multifab_build_edge(edge_coeffs(d), la,      nc=1,ng=nghost(beta(n,d)),dir=d)
          call multifab_copy_c(edge_coeffs(d),1,beta(n,d),1,nc=1,ng=nghost(beta(n,d)))
       end do

       ! xa and xb tell the stencil how far away the ghost cell values are in physical
       ! dimensions from the edge of the grid
       if (n > 1) then
          xa = HALF*mla%mba%rr(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
          xb = HALF*mla%mba%rr(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
       else
          if (base_level .gt. 1) then
             xa = HALF*crse_ratio*mgt(n)%dh(:,mgt(n)%nlevels)
             xb = HALF*crse_ratio*mgt(n)%dh(:,mgt(n)%nlevels)
          else
             xa = ZERO
             xb = ZERO
          end if
       end if

       pxa = ZERO
       pxb = ZERO

       call stencil_fill_cc(mgt(n)%ss(mgt(n)%nlevels), cell_coeffs, edge_coeffs, mgt(n)%dh(:,mgt(n)%nlevels), &
                            mgt(n)%mm(mgt(n)%nlevels), xa, xb, pxa, pxb, stencil_order, &
                            the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp))

       call multifab_destroy(cell_coeffs)

       do d = 1,dm
          call multifab_destroy(edge_coeffs(d))
       end do

    end do

    call ml_cc_applyop(mla, mgt, res, phi, mla%mba%rr)

    do n = 1, nlevs
       call mg_tower_destroy(mgt(n))
    end do

    call destroy(bpt)

  end subroutine cc_applyop

  subroutine ml_cc_applyop(mla, mgt, res, full_soln, ref_ratio)

    use bl_prof_module
    use ml_cc_restriction_module , only : ml_cc_restriction
    use ml_prolongation_module   , only : ml_interp_bcs
    use cc_ml_resid_module       , only : crse_fine_residual_cc

    type(ml_layout), intent(in)    :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: res(:)
    type( multifab), intent(inout) :: full_soln(:)
    integer        , intent(in   ) :: ref_ratio(:,:)

    integer :: nlevs
    type(multifab) , allocatable :: rh(:) ! this will be set to zero
    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box) :: pd, pdc
    type(layout) :: la, lac
    integer :: n, ng_fill
    integer :: mglev, mglev_crse

    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_cc_applyop")

    nlevs = mla%nlevel

    allocate(rh(nlevs))

    allocate(brs_flx(2:nlevs))
    allocate(brs_bcs(2:nlevs))

    do n = nlevs, 1, -1

       la = mla%la(n)
       call multifab_build(rh(n), la, 1, 0)
       call setval(rh(n), ZERO,all=.true.)

       ! zero residual just to be safe
       call setval(     res(n), ZERO,all=.true.)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2)
       call  flux_reg_build   (brs_flx(n), la, lac, ref_ratio(n-1,:), pdc)

    end do

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in 
    !  its ghost cells before we evaluate the initial residual  
    do n = 2,nlevs
       ng_fill = nghost(full_soln(n))
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(1,0), pd, &
            ref_ratio(n-1,:), ng_fill, brs_bcs(n)%facemap, brs_bcs(n)%indxmap, &
            brs_bcs(n)%uncovered)
    end do

    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n),mgt(n)%mm(mglev), &
                           mgt(n)%stencil_type, mgt(n)%lcross)
    end do

    ! Make sure to correct the coarse cells immediately next to fine grids
    !   using the averaged fine grid fluxes
    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc, &
                                  ref_ratio(n-1,:), filled=.true.)

       call ml_cc_restriction(res(n-1), res(n), ref_ratio(n-1,:))
    enddo

    ! still need to multiply residual by -1 to get (alpha - del dot beta grad)
    do n=1,nlevs
       call multifab_mult_mult_s_c(res(n),1,-one,1,0)
    end do

    do n = nlevs, 1, -1
       call multifab_destroy(rh(n))
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_bcs(n))
       call bndry_reg_destroy(brs_flx(n))
    end do

    call destroy(bpt)

 end subroutine ml_cc_applyop
!
! ******************************************************************************************
!

 ! subroutine ml_cc_n_applyop(mla, mgt, res, full_soln, ref_ratio)

 !    use bl_prof_module
 !    use ml_cc_restriction_module , only : ml_cc_restriction
 !    use ml_prolongation_module   , only : ml_interp_bcs
 !    use cc_ml_resid_module       , only : crse_fine_residual_n_cc

 !    type(ml_layout), intent(in)    :: mla
 !    type(mg_tower) , intent(inout) :: mgt(:)
 !    type( multifab), intent(inout) :: res(:)
 !    type( multifab), intent(inout) :: full_soln(:)
 !    integer        , intent(in   ) :: ref_ratio(:,:)

 !    integer :: nlevs
 !    type(multifab) , allocatable :: rh(:) ! this will be set to zero
 !    type(bndry_reg), allocatable :: brs_flx(:)
 !    type(bndry_reg), allocatable :: brs_bcs(:)

 !    type(box) :: pd, pdc
 !    type(layout) :: la, lac
 !    integer :: n, ng_fill, nComp
 !    integer :: mglev, mglev_crse

 !    type(bl_prof_timer), save :: bpt

 !    call build(bpt, "ml_cc_n_applyop")

 !    nlevs = mla%nlevel

 !    allocate(rh(nlevs))

 !    allocate(brs_flx(2:nlevs))
 !    allocate(brs_bcs(2:nlevs))

 !    nComp = 2

 !    do n = nlevs, 1, -1

 !       la = mla%la(n)
 !       call multifab_build(rh(n), la, 1, 0)
 !       call multifab_setval(rh(n), ZERO,all=.true.)

 !       ! zero residual just to be safe
 !       call setval(res(n), ZERO,all=.true.)

 !       if ( n == 1 ) exit

 !       ! Build the (coarse resolution) flux registers to be used in computing
 !       !  the residual at a non-finest AMR level.

 !       pdc = layout_get_pd(mla%la(n-1))
 !       lac = mla%la(n-1)
 !       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2)
 !       call  flux_reg_build   (brs_flx(n), la, lac, ref_ratio(n-1,:), pdc, nc = nComp)

 !    end do

 !    !  Make sure full_soln at fine grid has the correct coarse grid bc's in 
 !    !  its ghost cells before we evaluate the initial residual  
 !    do n = 2,nlevs
 !       ng_fill = nghost(full_soln(n))
 !       pd = layout_get_pd(mla%la(n))
 !       call multifab_fill_boundary(full_soln(n-1))
 !       call bndry_reg_copy(brs_bcs(n), full_soln(n-1), filled=.true.)
 !       call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(1,0), pd, &
 !            ref_ratio(n-1,:), ng_fill, brs_bcs(n)%facemap, brs_bcs(n)%indxmap, &
 !            brs_bcs(n)%uncovered)
 !    end do

 !    call multifab_fill_boundary(full_soln(nlevs))

 !    do n = 1,nlevs,1
 !       mglev = mgt(n)%nlevels
 !       call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n),mgt(n)%mm(mglev), &
 !                           mgt(n)%stencil_type, mgt(n)%lcross, filled=.true.) 
 !    end do

 !    ! Make sure to correct the coarse cells immediately next to fine grids
 !    !   using the averaged fine grid fluxes
 !    do n = nlevs,2,-1
 !       mglev      = mgt(n  )%nlevels
 !       mglev_crse = mgt(n-1)%nlevels

 !       pdc = layout_get_pd(mla%la(n-1))
 !       call crse_fine_residual_n_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc, &
 !                                    ref_ratio(n-1,:), filled=.true.)
 !       call ml_cc_restriction(res(n-1), res(n), ref_ratio(n-1,:))
 !    enddo

 !    ! still need to multiply residual by -1 to get (alpha - del dot beta grad)
 !    do n=1,nlevs
 !       call multifab_mult_mult_s_c(res(n),1,-one,1,0)
 !    end do

 !    do n = nlevs, 1, -1
 !       call multifab_destroy(rh(n))
 !       if ( n == 1 ) exit
 !       call bndry_reg_destroy(brs_bcs(n))
 !       call bndry_reg_destroy(brs_flx(n))
 !    end do

 !    call destroy(bpt)

 !  end subroutine ml_cc_n_applyop

end module cc_applyop_module
