
module compute_flux_module

  use multifab_module
  use ml_layout_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use ml_cc_restriction_module
  use bndry_reg_module
  use define_bc_module
  use bc_module

  implicit none

  private

  public :: compute_flux, compute_flux_single_level

contains

  subroutine compute_flux(mla,phi,velocity,flux,dx,dt)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi(:)
    type(multifab) , intent(in   ) :: velocity(:,:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:), dt(:)

    ! local variables
    integer :: dm, ng_p, ng_f, i, n, nlevs

    dm    = mla%dim
    nlevs = mla%nlevel

    ng_p = phi(1)%ng
    ng_f = flux(1,1)%ng

    do n=1,nlevs
       call compute_flux_single_level(mla,phi(n),velocity(n,:),flux(n,:),dx(n),dt(n))
    end do

    ! set level n-1 fluxes to be the average of the level n fluxes covering it
    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       do i=1,dm
          call ml_edge_restriction_c(flux(n-1,i),1,flux(n,i),1,mla%mba%rr(n-1,:),i,1)
       end do
    end do

  end subroutine compute_flux

  subroutine compute_flux_single_level(mla,phi,velocity,flux,dx,dt)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi
    type(multifab) , intent(in   ) :: velocity(:)
    type(multifab) , intent(inout) :: flux(:)
    real(kind=dp_t), intent(in   ) :: dx, dt

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, ng_u, i

    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)
    real(kind=dp_t), pointer ::  up(:,:,:,:)
    real(kind=dp_t), pointer ::  vp(:,:,:,:)
    real(kind=dp_t), pointer ::  wp(:,:,:,:)

    dm    = mla%dim

    ng_p = phi%ng
    ng_f = flux(1)%ng
    ng_u = velocity(1)%ng

    do i=1,nfabs(phi)
       pp  => dataptr(phi,i)
       fxp => dataptr(flux(1),i)
       fyp => dataptr(flux(2),i)
       up  => dataptr(velocity(1),i)
       vp  => dataptr(velocity(2),i)
       lo = lwb(get_box(phi,i))
       hi = upb(get_box(phi,i))
       select case(dm)
       case (2)
          call compute_flux_2d(pp(:,:,1,1), ng_p, &
                               up(:,:,1,1), vp(:,:,1,1), ng_u, &
                               fxp(:,:,1,1),  fyp(:,:,1,1), ng_f, &
                               lo, hi, dx, dt)

       case (3)
          fzp => dataptr(flux(3),i)
          wp  => dataptr(velocity(3),i)
          call compute_flux_3d(pp(:,:,:,1), ng_p, &
                               fxp(:,:,:,1),  fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                               lo, hi, dx, dt)
       end select
    end do

  end subroutine compute_flux_single_level

  subroutine compute_flux_2d(phi, ng_p, umac, vmac, ng_u, fluxx, fluxy, ng_f, &
                             lo, hi, dx, dt)
    
    use prob_module , only : mu
    use slope_module, only : slope_2d

    integer          :: lo(2), hi(2), ng_p, ng_f, ng_u
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision ::  umac(lo(1)-ng_u:,lo(2)-ng_u:)
    double precision ::  vmac(lo(1)-ng_u:,lo(2)-ng_u:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx, dt

    ! local variables
    integer          :: i,j
    double precision :: hdtdx
    double precision, allocatable :: slope(:,:)
    double precision, allocatable :: phix(:,:)
    double precision, allocatable :: phiy(:,:)
    double precision, allocatable :: phix_temp(:,:)
    double precision, allocatable :: phiy_temp(:,:)

    allocate(     slope(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(      phix(lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(      phiy(lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1))
    allocate(phix_temp(lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(phiy_temp(lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1))


    hdtdx = 0.5d0*dt/dx

    ! ****************************************************************
    ! Advective Flxues
    ! ****************************************************************

    ! Compute the slopes in the x-direction
    call slope_2d(phi  ,lo(1)-ng_p,lo(2)-ng_p,hi(1)+ng_p,hi(2)+ng_p, &
                  slope,lo(1)-1   ,lo(2)-1   ,hi(1)+1   ,hi(2)+1, &
                  lo(1),lo(2),hi(1),hi(2),1,1) 

    ! compute phi on x edges using umac to upwind
    do j=lo(2)-1,hi(2)+1
       do i=lo(1),hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix(i,j) = phi(i  ,j) - (0.5d0 + hdtdx*umac(i,j))*slope(i  ,j)
          else
             phix(i,j) = phi(i-1,j) + (0.5d0 - hdtdx*umac(i,j))*slope(i-1,j)
          end if

          phix_temp(i,j) = phix(i,j)

       end do
    end do

    ! compute the slopes in the y-direction
    call slope_2d(phi  ,lo(1)-ng_p,lo(2)-ng_p,hi(1)+ng_p,hi(2)+ng_p, &
                  slope,lo(1)-1   ,lo(2)-1   ,hi(1)+1   ,hi(2)+1, &
                  lo(1),lo(2),hi(1),hi(2),1,2) 

    ! compute phi on y edges using umac to upwind
    do j=lo(2),hi(2)+1
       do i=lo(1)-1,hi(1)+1

          if (vmac(i,j) .lt. 0.d0) then
             phiy(i,j) = phi(i,j  ) - (0.5d0 + hdtdx*vmac(i,j))*slope(i,j  )
          else
             phiy(i,j) = phi(i,j-1) + (0.5d0 - hdtdx*vmac(i,j))*slope(i,j-1)
          end if

          phiy_temp(i,j) = phiy(i,j)

       end do
    end do

    ! Use the fluxes on y-edges to add transverse contributions on x-edges
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix(i,j) = phix(i,j) - hdtdx*( 0.5d0*(vmac(i  ,j+1)+vmac(i  ,j)) * (phiy_temp(i  ,j+1)-phiy_temp(i  ,j)) )
          else
             phix(i,j) = phix(i,j) - hdtdx*( 0.5d0*(vmac(i-1,j+1)+vmac(i-1,j)) * (phiy_temp(i-1,j+1)-phiy_temp(i-1,j)) )
          end if

          ! compute final x-fluxes
          fluxx(i,j) = -phix(i,j)*umac(i,j)

       end do
    end do

    ! Use the fluxes on x-edges to add transverse contributions on y-edges
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)

          if (vmac(i,j) .lt. 0.d0) then
             phiy(i,j) = phiy(i,j) - hdtdx*( 0.5d0*(umac(i+1,j  )+umac(i,j  )) * (phix_temp(i+1,j  )-phix_temp(i,j  )) )
          else
             phiy(i,j) = phiy(i,j) - hdtdx*( 0.5d0*(umac(i+1,j-1)+umac(i,j-1)) * (phix_temp(i+1,j-1)-phix_temp(i,j-1)) )
          end if

          ! compute final y-fluxes
          fluxy(i,j) = -phiy(i,j)*vmac(i,j)

       end do
    end do

    deallocate(slope,phix,phiy,phix_temp,phiy_temp)

  end subroutine compute_flux_2d

  subroutine compute_flux_3d(phi, ng_p, fluxx, fluxy, fluxz, ng_f, &
                             lo, hi, dx, dt)

    use prob_module , only : mu, uadv, vadv, wadv
    use slope_module, only : slope_3d

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx, dt

    ! local variables
    integer          :: i,j,k

  end subroutine compute_flux_3d

end module compute_flux_module

