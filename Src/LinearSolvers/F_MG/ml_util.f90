module ml_util_module

  use stencil_module
  use stencil_nodal_module

  implicit none

contains

  subroutine ml_fill_fluxes(ss, flux, uu, mm, ratio, face, dim)
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    integer :: ratio
    integer :: face, dim
    integer :: i, n
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: ng, nc


    nc = uu%nc
    ng = uu%ng

    if ( uu%nc /= flux%nc ) then
       call bl_error("ML_FILL_FLUXES: uu%nc /= flux%nc")
    end if

    call multifab_fill_boundary(uu)
    do i = 1, flux%nboxes
       if ( multifab_remote(flux, i) ) cycle
       fp => dataptr(flux, i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       do n = 1, uu%nc
          select case(ss%dim)
          case (1)
             call stencil_flux_1d(sp(:,1,1,:), fp(:,1,1,n), up(:,1,1,n), &
                  mp(:,1,1,1), ng, ratio, face, dim)
          case (2)
             call stencil_flux_2d(sp(:,:,1,:), fp(:,:,1,n), up(:,:,1,n), &
                  mp(:,:,1,1), ng, ratio, face, dim)
          case (3)
             call stencil_flux_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
                  mp(:,:,:,1), ng, ratio, face, dim)
          end select
       end do
    end do
  end subroutine ml_fill_fluxes

  subroutine ml_fill_fluxes_c(ss, flux, cf, uu, cu, mm, ratio, face, dim)
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: cf, cu
    integer :: ratio
    integer :: face, dim
    integer :: i, n
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: ng


    ng = uu%ng

    if ( uu%nc /= flux%nc ) then
       call bl_error("ML_FILL_FLUXES: uu%nc /= flux%nc")
    end if

    call multifab_fill_boundary(uu)
    do i = 1, flux%nboxes
       if ( multifab_remote(flux, i) ) cycle
       fp => dataptr(flux, i, cf)
       up => dataptr(uu, i, cu)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       select case(ss%dim)
       case (1)
          call stencil_flux_1d(sp(:,1,1,:), fp(:,1,1,1), up(:,1,1,1), &
               mp(:,1,1,1), ng, ratio, face, dim)
       case (2)
          call stencil_flux_2d(sp(:,:,1,:), fp(:,:,1,1), up(:,:,1,1), &
               mp(:,:,1,1), ng, ratio, face, dim)
       case (3)
          call stencil_flux_3d(sp(:,:,:,:), fp(:,:,:,1), up(:,:,:,1), &
               mp(:,:,:,1), ng, ratio, face, dim)
       end select
    end do
  end subroutine ml_fill_fluxes_c

  subroutine ml_fine_contrib(ss, flux, uu, res, mm, ratio, crse_domain, side)
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(multifab), intent(inout) :: res
    type(imultifab), intent(in) :: mm
    type(box) :: crse_domain
    type(box) :: fbox
    integer :: side
    integer :: ratio(:)
    integer :: lof(flux%dim)
    integer :: lo_dom(flux%dim), hi_dom(flux%dim)
    integer :: i, n, dir
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: nc, ng

    nc = uu%nc
    ng = uu%ng

    if ( uu%nc /= flux%nc ) then
       call bl_error("ML_FILL_FLUXES: uu%nc /= flux%nc")
    end if

    lo_dom = lwb(crse_domain)
    hi_dom = upb(crse_domain)
    if ( nodal_q(uu) ) hi_dom = hi_dom + 1
    dir = iabs(side)

    call multifab_fill_boundary(uu)
    do i = 1, flux%nboxes
       if ( multifab_remote(flux, i) ) cycle
       fbox   = get_ibox(flux,i)
       lof = lwb(fbox)
       fp => dataptr(flux, i)
       up => dataptr(uu, i)
       rp => dataptr(res, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       do n = 1, nc
          if (lof(dir) /= lo_dom(dir) .and. lof(dir) /= hi_dom(dir)) then
             select case(flux%dim)
             case (1)
                call fine_edge_resid_1d(sp(:,1,1,:), fp(:,1,1,n), up(:,1,1,n), &
                     rp(:,1,1,1), mp(:,1,1,1), ng, ratio, side)
             case (2)
                call fine_edge_resid_2d(sp(:,:,1,:), fp(:,:,1,n), up(:,:,1,n), &
                     rp(:,:,1,1), mp(:,:,1,1), ng, ratio, side, lof)
             case (3)
                call fine_edge_resid_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
                     rp(:,:,:,1), mp(:,:,:,1), ng, ratio, side, lof)
             end select
          end if
       end do
    end do
  end subroutine ml_fine_contrib

end module ml_util_module
