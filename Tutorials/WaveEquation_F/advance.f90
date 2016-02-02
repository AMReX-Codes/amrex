module advance_module

  implicit none

  private

  public :: advance

contains
  
  subroutine advance(data,dx,dt)

    use multifab_module

    type(multifab) , intent(inout) :: data
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: dt

    ! local variables
    integer :: lo(data%dim), hi(data%dim)
    integer :: dm, ng, i

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    ! set these here so we don't have to pass them into the subroutine
    dm = data%dim
    ng = data%ng

    ! here, nfabs() return the number of boxes local to the MPI task,
    ! so this loop only goes over the boxes local to a processor.
    do i=1,nfabs(data)
       dp => dataptr(data,i)
       lo = lwb(get_box(data,i))
       hi = upb(get_box(data,i))
       select case(dm)
       case (2)
          call advance_2d(dp(:,:,1,:), ng, lo, hi, dx, dt)
       case (3)
          call advance_3d(dp(:,:,:,:), ng, lo, hi, dx, dt)
       end select
    end do
    
    ! fill ghost cells
    ! this only fills periodic ghost cells and ghost cells for neighboring
    ! grids at the same level.  Physical boundary ghost cells are filled
    ! using multifab_physbc.  But this problem is periodic, so this
    ! call is sufficient.
    call multifab_fill_boundary(data)

  end subroutine advance

  subroutine advance_2d(U, ng, lo, hi, dx, dt)

    integer          :: lo(2), hi(2), ng
    double precision :: U(lo(1)-ng:,lo(2)-ng:,:)
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
    double precision :: U(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
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

