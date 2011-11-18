module advance_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: advance

contains
  
  subroutine advance(mla,dx,dt,data)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: data(mla%nlevel)
    real(kind=dp_t), intent(in   ) :: dx(mla%nlevel,mla%dim)
    real(kind=dp_t), intent(in   ) :: dt

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, n, i

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    ! set these here so we don't have to pass them into the subroutine
    nlevs = mla%nlevel
    dm    = mla%dim
    ng    = nghost(data(1))

    do n=1,nlevs

       do i=1,nboxes(data(n))
          if ( multifab_remote(data(n),i) ) cycle
          dp => dataptr(data(n),i)
          lo = lwb(get_box(data(n),i))
          hi = upb(get_box(data(n),i))
          select case(dm)
          case (2)
             call advance_2d(dp(:,:,1,:), ng, lo, hi, dx(n,1), dt)
          case (3)
             call advance_3d(dp(:,:,:,:), ng, lo, hi, dx(n,1), dt)
          end select
       end do

       ! fill ghost cells
       ! this only fills periodic ghost cells and ghost cells for neighboring
       ! grids at the same level.  Physical boundary ghost cells are filled
       ! using multifab_physbc.  But this problem is periodic, so this
       ! call is sufficient.
       call multifab_fill_boundary(data(n))

    end do

  end subroutine advance

  subroutine advance_2d(U, ng, lo, hi, dx, dt)

    integer          :: lo(2), hi(2), ng
    double precision :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,2)
    double precision :: dx, dt

    double precision, allocatable :: dU(:,:,:), Unew(:,:,:)

    double precision, parameter :: third     = 1.d0/3.d0
    double precision, parameter :: twothirds = 2.d0/3.d0
    double precision, parameter :: fac1      = -1.d0/12.d0
    double precision, parameter :: fac2      =  4.d0/3.d0
    double precision, parameter :: fac3      = -5.0d0
    
    integer          :: i,j,ndo

    allocate(  dU(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,2))
    allocate(Unew(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,2))

    if (ng < 6) then
       print *,'NOT ENOUGH GHOST CELLS IN ADVANCE ',ng
       stop
    end if
    
    ! First call to get RHS = Lap(U)
    ndo = 4

    do j = lo(2)-ndo,hi(2)+ndo
       do i = lo(1)-ndo,hi(1)+ndo

          dU(i,j,1) = U(i,j,2)

          dU(i,j,2) = (  fac1 * ( U(i-2,j,1) + U(i+2,j,1) + &
                                  U(i,j-2,1) + U(i,j+2,1) ) &
                       + fac2 * ( U(i-1,j,1) + U(i+1,j,1) + &
                                  U(i,j-1,1) + U(i,j+1,1) ) &
                       + fac3 *   U(i,j,1)  ) / (dx*dx)
       end do
    end do

    ! First update
    ! This Unew lives at t^{n+1}
    Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2) = & 
                 U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2) &
         + dt * dU(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2)

    ! Second call to get RHS = Lap(Unew)
    ndo = 2

    do j = lo(2)-ndo,hi(2)+ndo
       do i = lo(1)-ndo,hi(1)+ndo

          dU(i,j,1) = Unew(i,j,2)

          dU(i,j,2) = (  fac1 * ( Unew(i-2,j,1) + Unew(i+2,j,1) + &
                                  Unew(i,j-2,1) + Unew(i,j+2,1) ) &
                       + fac2 * ( Unew(i-1,j,1) + Unew(i+1,j,1) + &
                                  Unew(i,j-1,1) + Unew(i,j+1,1) ) &
                       + fac3 *   Unew(i,j,1)  ) / (dx*dx)
       end do
    end do

    ! Second update
    ! This Unew lives at t^{n+1/2}
    Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2) =        &
                .75d0 *    U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2)   &
              + .25d0 * Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2)   &
         + dt * .25d0 *   dU(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2)

    ! Third call to get RHS = Lap(Unew)
    ndo = 0

    do j = lo(2)-ndo,hi(2)+ndo
       do i = lo(1)-ndo,hi(1)+ndo

          dU(i,j,1) = Unew(i,j,2)

          dU(i,j,2) = (  fac1 * ( Unew(i-2,j,1) + Unew(i+2,j,1) + &
                                  Unew(i,j-2,1) + Unew(i,j+2,1) ) &
                       + fac2 * ( Unew(i-1,j,1) + Unew(i+1,j,1) + &
                                  Unew(i,j-1,1) + Unew(i,j+1,1) ) &
                       + fac3 *   Unew(i,j,1)  ) / (dx*dx)
       end do
    end do

    ! Third update
    ! This Unew lives at t^{n+1}
    U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2) =        &
                third     *    U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2)   &
              + twothirds * Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2)   &
         + dt * twothirds *   dU(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,1:2)

    deallocate(dU,Unew)

  end subroutine advance_2d

  subroutine advance_3d(U, ng, lo, hi, dx, dt)

    integer          :: lo(3), hi(3), ng
    double precision :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,2)
    double precision :: dx, dt

    double precision, allocatable :: dU(:,:,:,:), Unew(:,:,:,:)

    double precision, parameter :: third     = 1.d0/3.d0
    double precision, parameter :: twothirds = 2.d0/3.d0
    double precision, parameter :: fac1      = -1.d0/12.d0
    double precision, parameter :: fac2      =  4.d0/3.d0
    double precision, parameter :: fac3      = -7.5d0
    
    integer          :: i,j,k,ndo

    allocate(  dU(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,2))
    allocate(Unew(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,2))

    if (ng < 6) then
       print *,'NOT ENOUGH GHOST CELLS IN ADVANCE ',ng
       stop
    end if
 
    ! First call to get RHS = Lap(U)
    ndo = 4

    do k = lo(3)-ndo,hi(3)+ndo
       do j = lo(2)-ndo,hi(2)+ndo
          do i = lo(1)-ndo,hi(1)+ndo

             dU(i,j,k,1) = U(i,j,k,2)
  
             dU(i,j,k,2) = (  fac1 * ( U(i-2,j,k,1) + U(i+2,j,k,1) + &
                                       U(i,j-2,k,1) + U(i,j+2,k,1) + &
                                       U(i,j,k-2,1) + U(i,j,k+2,1) ) &
                            + fac2 * ( U(i-1,j,k,1) + U(i+1,j,k,1) + &
                                       U(i,j-1,k,1) + U(i,j+1,k,1) + &
                                       U(i,j,k-1,1) + U(i,j,k+1,1) ) &
                            + fac3 *   U(i,j,k,1)  ) / (dx*dx)
          end do
       end do
    end do

    ! First update
    ! This Unew lives at t^{n+1}
    Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) = & 
                 U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) &
         + dt * dU(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2)
 
    ! Second call to get RHS = Lap(Unew)
    ndo = 2

    do k = lo(3)-ndo,hi(3)+ndo
       do j = lo(2)-ndo,hi(2)+ndo
          do i = lo(1)-ndo,hi(1)+ndo

             dU(i,j,k,1) = Unew(i,j,k,2)

             dU(i,j,k,2) = (  fac1 * ( Unew(i-2,j,k,1) + Unew(i+2,j,k,1) + &
                                       Unew(i,j-2,k,1) + Unew(i,j+2,k,1) + &
                                       Unew(i,j,k-2,1) + Unew(i,j,k+2,1) ) &
                            + fac2 * ( Unew(i-1,j,k,1) + Unew(i+1,j,k,1) + &
                                       Unew(i,j-1,k,1) + Unew(i,j+1,k,1) + &
                                       Unew(i,j,k-1,1) + Unew(i,j,k+1,1) ) &
                            + fac3 *   Unew(i,j,k,1)  ) / (dx*dx)
          end do
       end do
    end do

    ! Second update
    ! This Unew lives at t^{n+1/2}
    Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) = & 
                .75d0 *    U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) &
              + .25d0 * Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) &
         + dt * .25d0 *   dU(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) 

    ! Third call to get RHS = Lap(Unew)
    ndo = 0

    do k = lo(3)-ndo,hi(3)+ndo
       do j = lo(2)-ndo,hi(2)+ndo
          do i = lo(1)-ndo,hi(1)+ndo

             dU(i,j,k,1) = Unew(i,j,k,2)

             dU(i,j,k,2) = (  fac1 * ( Unew(i-2,j,k,1) + Unew(i+2,j,k,1) + &
                                       Unew(i,j-2,k,1) + Unew(i,j+2,k,1) + &
                                       Unew(i,j,k-2,1) + Unew(i,j,k+2,1) ) &
                            + fac2 * ( Unew(i-1,j,k,1) + Unew(i+1,j,k,1) + &
                                       Unew(i,j-1,k,1) + Unew(i,j+1,k,1) + &
                                       Unew(i,j,k-1,1) + Unew(i,j,k+1,1) ) &
                            + fac3 *   Unew(i,j,k,1)  ) / (dx*dx)

          end do
       end do
    end do

    ! Third update
    ! This Unew lives at t^{n+1}
    U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2) =        &
                third     *    U(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2)   &
              + twothirds * Unew(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2)   &
         + dt * twothirds *   dU(lo(1)-ndo:hi(1)+ndo,lo(2)-ndo:hi(2)+ndo,lo(3)-ndo:hi(3)+ndo,1:2)

    deallocate(dU,Unew)

  end subroutine advance_3d

end module advance_module

