module ml_util_module

  use stencil_module
  use mg_module

  implicit none

contains

  subroutine ml_fill_fluxes(mgt, ss, flux, uu, mm, ratio, face, dim)
    type(mg_tower), intent(inout) :: mgt
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

    call multifab_fill_boundary(uu)
    do i = 1, flux%nboxes
       if ( multifab_remote(flux, i) ) cycle
       fp => dataptr(flux, i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       n = 1
       select case(mgt%dim)
       case (1)
          call stencil_flux_1d(sp(:,1,1,:), fp(:,1,1,n), up(:,1,1,n), &
               mp(:,1,1,1), mgt%ng, ratio, face, dim)
       case (2)
          call stencil_flux_2d(sp(:,:,1,:), fp(:,:,1,n), up(:,:,1,n), &
               mp(:,:,1,1), mgt%ng, ratio, face, dim)
       case (3)
          call stencil_flux_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
               mp(:,:,:,1), mgt%ng, ratio, face, dim)
       end select
    end do
  end subroutine ml_fill_fluxes

  subroutine ml_fine_contrib(mgt, ss, flux, uu, res, mm, ratio, crse_domain, side)
    type(mg_tower), intent(inout) :: mgt
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(multifab), intent(inout) :: res
    type(imultifab), intent(in) :: mm
    type(box) :: crse_domain
    type(box) :: fbox
    integer :: side
    integer :: ratio(mgt%dim)
    integer :: lof(mgt%dim)
    integer :: lo_dom(mgt%dim), hi_dom(mgt%dim)
    integer :: i, n, dir
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)

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
       n = 1
       if (lof(dir) /= lo_dom(dir) .and. lof(dir) /= hi_dom(dir)) then
       select case(mgt%dim)
       case (1)
          call fine_edge_resid_1d(sp(:,1,1,:), fp(:,1,1,n), up(:,1,1,n), &
               rp(:,1,1,1), mp(:,1,1,1), mgt%ng, ratio, side)
       case (2)
          call fine_edge_resid_2d(sp(:,:,1,:), fp(:,:,1,n), up(:,:,1,n), &
               rp(:,:,1,1), mp(:,:,1,1), mgt%ng, ratio, side, lof)
       case (3)
          call fine_edge_resid_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
               rp(:,:,:,1), mp(:,:,:,1), mgt%ng, ratio, side, lof)
       end select
       end if
    end do
  end subroutine ml_fine_contrib

end module ml_util_module
