module advance_module

  use multifab_module
  use layout_module

  implicit none

  private

  public :: advance

contains
  
  subroutine advance(nlevs,data,flux,dx,dt,coef)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: data(:)
    type(multifab) , intent(in   ) :: flux(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt
    real(kind=dp_t), intent(in   ) :: coef

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
          call advance_2d(dp(:,:,1,1), ng, fluxx(:,:,1,1), fluxy(:,:,1,1), ng_f, lo, hi, dx(n), dt, coef)
       case (3)
          fluxz => dataptr(flux(n,3),i)
          call advance_3d(dp(:,:,:,1), ng, fluxx(:,:,:,1), fluxy(:,:,:,1), fluxz(:,:,:,1), ng_f, lo, hi, dx(n), dt, coef)
       end select
     end do
    end do
    
    ! fill ghost cells
    ! this only fills periodic ghost cells and ghost cells for neighboring
    ! grids at the same level.  Physical boundary ghost cells are filled
    ! using multifab_physbc.  But this problem is periodic, so this
    ! call is sufficient.
    do n = 1, nlevs
       call multifab_fill_boundary(data(n))
    end do

  end subroutine advance

  subroutine advance_2d(U, ng, fluxx, fluxy, ng_f, lo, hi, dx, dt, coef)

    integer          :: lo(2), hi(2), ng, ng_f
    double precision ::     U(lo(1)-ng  :,lo(2)-ng  :)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx, dt, coef
    
    integer          :: i,j
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          U(i,j) = U(i,j) + coef * dt * ( &
              (fluxx(i+1,j)-fluxx(i,j)) / dx + &
              (fluxy(i,j+1)-fluxy(i,j)) / dx )  

       end do
    end do

  end subroutine advance_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine advance_3d(U, ng, fluxx, fluxy, fluxz, ng_f, lo, hi, dx, dt, coef)

    integer          :: lo(3), hi(3), ng, ng_f
    double precision ::     U(lo(1)-ng  :,lo(2)-ng  :,lo(3)-ng:  )
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx, dt, coef
    
    integer          :: i,j,k
    
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          U(i,j,k) = U(i,j,k) + coef * dt * ( &
              (fluxx(i+1,j,k)-fluxx(i,j,k)) / dx + &
              (fluxy(i,j+1,k)-fluxy(i,j,k)) / dx + &
              (fluxz(i,j,k+1)-fluxz(i,j,k)) / dx )  

       end do
    end do
    end do

  end subroutine advance_3d


end module advance_module

