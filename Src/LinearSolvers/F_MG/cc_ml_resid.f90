module cc_ml_resid_module

  use bl_constants_module
  use mg_module
  use ml_layout_module
  use bndry_reg_module

  implicit none

  private :: ml_fill_fluxes, ml_fill_n_fluxes

contains

  subroutine crse_fine_residual_cc(n, mgt, uu, crse_res, brs_flx, pdc, ref_ratio)

      use cc_interface_stencil_module, only : ml_interface

      integer        , intent(in   ) :: n
      type(mg_tower) , intent(inout) :: mgt(:)
      type(bndry_reg), intent(inout) :: brs_flx
      type(multifab) , intent(inout) :: uu(:)
      type(multifab) , intent(inout) :: crse_res
      type(box)      , intent(in   ) :: pdc
      integer        , intent(in   ) :: ref_ratio(:)

      integer :: i, dm, mglev

      dm = brs_flx%dim
      mglev = mgt(n)%nlevels

      call multifab_fill_boundary(uu(n))

      do i = 1, dm
         call ml_fill_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,0), &
              uu(n), mgt(n)%mm(mglev), ref_ratio(i), -1, i)
         call ml_fill_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,1), &
              uu(n), mgt(n)%mm(mglev), ref_ratio(i), 1, i)
      end do
      call bndry_reg_copy_to_other(brs_flx)
      do i = 1, dm
         call ml_interface(crse_res, brs_flx%obmf(i,0), uu(n-1), &
              mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, -1, i, ONE)
         call ml_interface(crse_res, brs_flx%obmf(i,1), uu(n-1), &
              mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, +1, i, ONE)
      end do

  end subroutine crse_fine_residual_cc

  subroutine crse_fine_residual_n_cc(n, mgt, uu, crse_res, brs_flx, pdc, ref_ratio)

      use cc_interface_stencil_module, only : ml_interface

      integer        , intent(in   ) :: n
      type(mg_tower) , intent(inout) :: mgt(:)
      type(bndry_reg), intent(inout) :: brs_flx
      type(multifab) , intent(inout) :: uu(:)
      type(multifab) , intent(inout) :: crse_res
      type(box)      , intent(in   ) :: pdc
      integer        , intent(in   ) :: ref_ratio(:)

      integer :: i, dm, mglev

      dm = brs_flx%dim
      mglev = mgt(n)%nlevels

      call multifab_fill_boundary(uu(n))

      do i = 1, dm
         call ml_fill_n_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,0), &
              uu(n), mgt(n)%mm(mglev), ref_ratio(i), -1, i)
         call ml_fill_n_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,1), &
              uu(n), mgt(n)%mm(mglev), ref_ratio(i), 1, i)
      end do
      call bndry_reg_copy_to_other(brs_flx)
      do i = 1, dm
         call ml_interface(crse_res, brs_flx%obmf(i,0), uu(n-1), &
              mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, -1, i, ONE)
         call ml_interface(crse_res, brs_flx%obmf(i,1), uu(n-1), &
              mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, +1, i, ONE)
      end do

  end subroutine crse_fine_residual_n_cc

!
! ******************************************************************************************
!

  subroutine ml_resid(mla, mgt, rh, res, full_soln, ref_ratio)

    use ml_restriction_module  , only : ml_restriction
    use ml_prolongation_module , only : ml_interp_bcs
    use stencil_defect_module  , only : compute_defect

    type(ml_layout), intent(in)    :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: rh(:)
    type( multifab), intent(inout) :: res(:)
    type( multifab), intent(inout) :: full_soln(:)
    integer        , intent(in   ) :: ref_ratio(:,:)

    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box)    :: pd, pdc
    type(layout) :: la, lac
    integer      :: i, n, ng_fill, dm, nlevs, mglev, mglev_crse

    dm = get_dim(rh(1))

    nlevs = mla%nlevel

    allocate(brs_bcs(2:nlevs))
    allocate(brs_flx(2:nlevs))

    do n = nlevs, 1, -1

       la = mla%la(n)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_flx(n), la, lac, ref_ratio(n-1,:), pdc, width = 0)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2, other = .false.)

    end do

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in its ghost cells 
    !   before we evaluate the initial residual  
    do n = 2,nlevs
       ng_fill = nghost(full_soln(n))
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, ref_ratio(n-1,:), ng_fill, -i)
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, ref_ratio(n-1,:), ng_fill, +i)
       end do
       call multifab_fill_boundary(full_soln(n))
    end do

    !   Make sure all periodic and internal boundaries are filled as well
    do n = 1,nlevs   
       call multifab_fill_boundary(full_soln(n))
    end do

    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call compute_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n),mgt(n)%mm(mglev), &
                           mgt(n)%stencil_type,mgt(n)%lcross)
    end do

    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc,ref_ratio(n-1,:))

       call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))
    enddo

    do n = nlevs, 1, -1
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_flx(n))
       call bndry_reg_destroy(brs_bcs(n))
    end do

  end subroutine ml_resid

!
! ******************************************************************************************
!
  subroutine ml_fill_fluxes(ss, flux, uu, mm, ratio, face, dim)

    use bl_prof_module
    use cc_stencil_apply_module

    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    integer :: ratio
    integer :: face, dim
    integer :: i, n, dm
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: ng
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_fill_fluxes")

    ng = nghost(uu)

    if ( ncomp(uu) /= ncomp(flux) ) then
       call bl_error("ML_FILL_FLUXES: uu%nc /= flux%nc")
    end if

    dm = get_dim(ss)

    !$OMP PARALLEL DO PRIVATE(i,fp,up,sp,mp,n)
    do i = 1, nfabs(flux)
       fp => dataptr(flux, i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       do n = 1, ncomp(uu)
          select case(dm)
          case (1)
             call stencil_flux_1d(sp(:,:,1,1), fp(:,1,1,n), up(:,1,1,n), &
                  mp(:,1,1,1), ng, ratio, face, dim)
          case (2)
             call stencil_flux_2d(sp(:,:,:,1), fp(:,:,1,n), up(:,:,1,n), &
                  mp(:,:,1,1), ng, ratio, face, dim)
          case (3)
             call stencil_flux_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
                  mp(:,:,:,1), ng, ratio, face, dim)
          end select
       end do
    end do
    !$OMP END PARALLEL DO
    call destroy(bpt)
  end subroutine ml_fill_fluxes
!
! ******************************************************************************************
!
  subroutine ml_fill_fluxes_c(ss, flux, cf, uu, cu, mm, ratio, face, dim, lcross)

    use bl_prof_module
    use cc_stencil_apply_module

    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: cf, cu
    integer, intent(in) :: ratio
    integer, intent(in) :: face, dim
    logical, intent(in) :: lcross

    integer                   :: i,ng
    type(bl_prof_timer), save :: bpt

    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)

    call build(bpt, "ml_fill_fluxes_c")

    ng = nghost(uu)

    call multifab_fill_boundary(uu, cross = lcross)

    !$OMP PARALLEL DO PRIVATE(i,fp,up,sp,mp)
    do i = 1, nfabs(flux)
       fp => dataptr(flux, i, cf)
       up => dataptr(uu, i, cu)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       select case(get_dim(ss))
       case (1)
          call stencil_flux_1d(sp(:,:,1,1), fp(:,1,1,1), up(:,1,1,1), &
               mp(:,1,1,1), ng, ratio, face, dim)
       case (2)
          call stencil_flux_2d(sp(:,:,:,1), fp(:,:,1,1), up(:,:,1,1), &
               mp(:,:,1,1), ng, ratio, face, dim)
       case (3)
          call stencil_flux_3d(sp(:,:,:,:), fp(:,:,:,1), up(:,:,:,1), &
               mp(:,:,:,1), ng, ratio, face, dim)
       end select
    end do
    !$OMP END PARALLEL DO
    call destroy(bpt)
  end subroutine ml_fill_fluxes_c
!
! ******************************************************************************************
!
  subroutine ml_fill_n_fluxes(ss, flux, uu, mm, ratio, face, dim)

    use bl_prof_module
    use cc_stencil_apply_module

    type(multifab), intent(inout) :: flux
    type(multifab), intent(in)    :: ss
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in)   :: mm
    integer :: ratio
    integer :: face, dim
    integer :: i, n
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: ng
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_fill_fluxes")

    ng = nghost(uu)

    !$OMP PARALLEL DO PRIVATE(i,fp,up,sp,mp,n)
    do i = 1, nfabs(flux)
       fp => dataptr(flux, i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       do n = 1, ncomp(uu)
          select case(get_dim(ss))
          case (1)
             call stencil_flux_1d(sp(:,:,1,1), fp(:,1,1,n), up(:,1,1,n), &
                                  mp(:,1,1,1), ng, ratio, face, dim)
          case (2)
             if ( ncomp(flux) > 1 ) then
                call stencil_flux_n_2d(sp(:,:,:,1), fp(:,:,1,:), up(:,:,1,n), &
                                     mp(:,:,1,1), ng, ratio, face, dim)
             else
                call stencil_flux_2d(sp(:,:,1,:), fp(:,:,1,n), up(:,:,1,n), &
                                     mp(:,:,1,1), ng, ratio, face, dim)
             end if
          case (3)
             call stencil_flux_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
                                  mp(:,:,:,1), ng, ratio, face, dim)
          end select
       end do
    end do
    !$OMP END PARALLEL DO
    call destroy(bpt)
  end subroutine ml_fill_n_fluxes

end module cc_ml_resid_module
