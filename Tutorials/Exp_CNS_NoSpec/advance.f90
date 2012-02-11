 module advance_module

  use bl_error_module
  use multifab_module

  implicit none

  private

  public :: advance

contains
  
  subroutine advance(U,dt)

    type(multifab),   intent(inout) :: U
    double precision, intent(in   ) :: dt

    double precision, parameter :: OneThird      = 1.d0/3.d0
    double precision, parameter :: TwoThirds     = 2.d0/3.d0
    double precision, parameter :: OneQuarter    = 1.d0/4.d0
    double precision, parameter :: ThreeQuarters = 3.d0/4.d0

    integer        :: lo(U%dim), hi(U%dim)
    integer        :: i, j, k, m, n, nc, ng
    type(layout)   :: la
    type(multifab) :: D, F, Unew

    double precision, pointer, dimension(:,:,:,:) :: up, dp, fp, unp

    nc = ncomp(U)
    ng = nghost(U)
    la = get_layout(U)
    !
    ! Sync ghost cells and periodic boundary cells in U prior to calculating D & F
    !
    call multifab_fill_boundary(U)

    call multifab_build(D,    la, nc, ng)
    call multifab_build(F,    la, nc, ng)
    call multifab_build(Unew, la, nc, ng)
    !
    ! Calculate D at time N.
    !
    do n=1,nboxes(D)
       if ( remote(D,n) ) cycle

       dp => dataptr(D,n)

       lo = lwb(get_box(D,n))
       hi = upb(get_box(D,n))
       !
       ! Use U.
       !
       dp = 0.0d0
    end do
    !
    ! Calculate F at time N.
    !
    do n=1,nboxes(F)
       if ( remote(F,n) ) cycle

       fp => dataptr(F,n)

       lo = lwb(get_box(F,n))
       hi = upb(get_box(F,n))
       !
       ! Use U.
       !
       fp = 0.0d0
    end do
    !
    ! Calculate U at time N+1/3.
    !
    do n=1,nboxes(U)
       if ( remote(U,n) ) cycle

       dp  => dataptr(D,   n)
       fp  => dataptr(F,   n)
       up  => dataptr(U,   n)
       unp => dataptr(Unew,n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       do m = 1, nc
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   unp(i,j,k,m) = up(i,j,k,m) + dt * (dp(i,j,k,m) + fp(i,j,k,m))
                end do
             end do
          end do
       end do
    end do
    !
    ! Sync U^1/3 prior to calculating D & F.
    !
    call multifab_fill_boundary(Unew)
    !
    ! Calculate D at time N+1/3.
    !
    do n=1,nboxes(D)
       if ( remote(D,n) ) cycle

       dp => dataptr(D,n)

       lo = lwb(get_box(D,n))
       hi = upb(get_box(D,n))
       !
       ! Use Unew.
       !
       dp = 0.0d0
    end do
    !
    ! Calculate F at time N+1/3.
    !
    do n=1,nboxes(F)
       if ( remote(F,n) ) cycle

       fp => dataptr(F,n)

       lo = lwb(get_box(F,n))
       hi = upb(get_box(F,n))
       !
       ! Use Unew.
       !
       fp = 0.0d0
    end do
    !
    ! Calculate U at time N+2/3.
    !
    do n=1,nboxes(U)
       if ( remote(U,n) ) cycle

       dp  => dataptr(D,   n)
       fp  => dataptr(F,   n)
       up  => dataptr(U,   n)
       unp => dataptr(Unew,n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       do m = 1, nc
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   unp(i,j,k,m) = ThreeQuarters * up(i,j,k,m) + &
                        OneQuarter * (unp(i,j,k,m) + dt * (dp(i,j,k,m) + fp(i,j,k,m)))
                end do
             end do
          end do
       end do
    end do
    !
    ! Sync U^2/3 prior to calculating D & F.
    !
    call multifab_fill_boundary(Unew)
    !
    ! Calculate D at time N+2/3.
    !
    do n=1,nboxes(D)
       if ( remote(D,n) ) cycle

       dp => dataptr(D,n)

       lo = lwb(get_box(D,n))
       hi = upb(get_box(D,n))
       !
       ! Use Unew.
       !
       dp = 0.0d0
    end do
    !
    ! Calculate F at time N+2/3.
    !
    do n=1,nboxes(F)
       if ( remote(F,n) ) cycle

       fp => dataptr(F,n)

       lo = lwb(get_box(F,n))
       hi = upb(get_box(F,n))
       !
       ! Use Unew.
       !
       fp = 0.0d0
    end do
    !
    ! Calculate U at time N+1.
    !
    do n=1,nboxes(U)
       if ( remote(U,n) ) cycle

       dp  => dataptr(D,   n)
       fp  => dataptr(F,   n)
       up  => dataptr(U,   n)
       unp => dataptr(Unew,n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       do m = 1, nc
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   up(i,j,k,m) = OneThird * up(i,j,k,m) + &
                        TwoThirds * (unp(i,j,k,m) + dt * (dp(i,j,k,m) + fp(i,j,k,m)))
                end do
             end do
          end do
       end do
    end do

    call destroy(Unew)
    call destroy(F)
    call destroy(D)

  end subroutine advance

end module advance_module

