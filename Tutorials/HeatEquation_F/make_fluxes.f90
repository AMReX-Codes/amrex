module make_flux_module

  use multifab_module
  use layout_module

  implicit none

  private

  public :: make_fluxes

contains
  
  subroutine make_fluxes(nlevs,data,flux,dx,ref_ratio)

    use ml_restriction_module , only : ml_edge_restriction

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: data(:)
    type(multifab) , intent(inout) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: ref_ratio(:,:)

    ! local variables
    integer :: lo(data(1)%dim), hi(data(1)%dim)
    integer :: dm, ng, ng_f, i, n

    real(kind=dp_t), pointer ::    dp(:,:,:,:)
    real(kind=dp_t), pointer :: fluxx(:,:,:,:)
    real(kind=dp_t), pointer :: fluxy(:,:,:,:)
    real(kind=dp_t), pointer :: fluxz(:,:,:,:)

    ! Set these here so we don't have to pass them into the subroutine
    dm   = data(1)%dim
    ng   = data(1)%ng
    ng_f = flux(1,1)%ng

    do n = 1, nlevs
       do i=1,nboxes(data(n))
          if ( multifab_remote(data(n),i) ) cycle
          dp => dataptr(data(n),i)
          fluxx => dataptr(flux(n,1),i)
          fluxy => dataptr(flux(n,2),i)
          lo = lwb(get_box(data(n),i))
          hi = upb(get_box(data(n),i))
          select case(dm)
          case (2)
             call make_flux_2d(dp(:,:,1,1), ng, fluxx, fluxy, ng_f, lo, hi, dx(n))
          case (3)
             fluxz => dataptr(flux(n,3),i)
             call make_flux_3d(dp(:,:,:,1), ng, fluxx, fluxy, fluxz, ng_f, lo, hi, dx(n))
          end select
       end do
    end do

    do n = nlevs,2,-1
       do i = 1,dm
          call ml_edge_restriction(flux(n-1,i),flux(n,i),ref_ratio(n-1,:),i)
       end do
    end do
    
  end subroutine make_fluxes

  subroutine make_flux_2d(U, ng, fluxx, fluxy, ng_f, lo, hi, dx)

    integer          :: lo(2), hi(2), ng, ng_f
    double precision ::     U(lo(1)-ng  :hi(1)+ng    ,lo(2)-ng  :hi(2)+ng)
    double precision :: fluxx(lo(1)-ng_f:hi(1)+ng_f+1,lo(2)-ng_f:hi(2)+ng_f  )
    double precision :: fluxy(lo(1)-ng_f:hi(1)+ng_f  ,lo(2)-ng_f:hi(2)+ng_f+1)
    double precision :: dx
    
    integer          :: i,j
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          fluxx(i,j) = (U(i,j) - U(i-1,j))/dx
       end do
    end do

    do j = lo(2),hi(2)+1
       do i = lo(1),hi(1)
          fluxy(i,j) = (U(i,j) - U(i,j-1))/dx
       end do
    end do

  end subroutine make_flux_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_flux_3d(U, ng, fluxx, fluxy, fluxz, ng_f, lo, hi, dx)

    integer          :: lo(3), hi(3), ng, ng_f
    double precision ::     U(lo(1)-ng  :hi(1)+ng    ,lo(2)-ng  :hi(2)+ng    ,lo(3)-ng  :hi(3)+ng)
    double precision :: fluxx(lo(1)-ng_f:hi(1)+ng_f+1,lo(2)-ng_f:hi(2)+ng_f  ,lo(3)-ng_f:hi(3)+ng_f  )
    double precision :: fluxy(lo(1)-ng_f:hi(1)+ng_f  ,lo(2)-ng_f:hi(2)+ng_f+1,lo(3)-ng_f:hi(3)+ng_f  )
    double precision :: fluxz(lo(1)-ng_f:hi(1)+ng_f  ,lo(2)-ng_f:hi(2)+ng_f  ,lo(3)-ng_f:hi(3)+ng_f+1)
    double precision :: dx
    
    integer          :: i,j,k

    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          fluxx(i,j,k) = (U(i,j,k) - U(i-1,j,k))/dx
       end do
    end do
    end do

    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          fluxy(i,j,k) = (U(i,j,k) - U(i,j-1,k))/dx
       end do
    end do
    end do

    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)+1
          fluxz(i,j,k) = (U(i,j,k) - U(i,j,k-1))/dx
       end do
    end do
    end do

  end subroutine make_flux_3d

end module make_flux_module
