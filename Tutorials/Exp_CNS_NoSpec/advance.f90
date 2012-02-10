module advance_module

  use bl_error_module
  use multifab_module

  implicit none

  private

  public :: advance

contains
  
  subroutine advance(U,dx,dt)

    type(multifab),   intent(inout) :: U
    double precision, intent(in   ) :: dx
    double precision, intent(in   ) :: dt

    double precision, parameter :: OneThird      = 1.d0/3.d0
    double precision, parameter :: TwoThirds     = 2.d0/3.d0
    double precision, parameter :: OneQuarter    = 1.d0/4.d0
    double precision, parameter :: ThreeQuarters = 3.d0/4.d0

    integer        :: lo(U%dim), hi(U%dim)
    integer        :: i, j, k, m, n, nc, ng
    type(layout)   :: la
    type(multifab) :: D, F, U13, U23

    double precision, pointer, dimension(:,:,:,:) ::   up, dp, fp, u13p, u23p

    nc = ncomp(U)
    ng = nghost(U)
    la = get_layout(U)
    !
    ! Sync up ghost cells and periodic boundary cells in U.
    !
    call multifab_fill_boundary(U)

    call multifab_build(D,  la,nc,ng)
    call multifab_build(F,  la,nc,ng)
    call multifab_build(U13,la,nc,ng)
    call multifab_build(U23,la,nc,ng)
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
    do n=1,nboxes(U13)
       if ( remote(U13,n) ) cycle

       dp   => dataptr(D,n)
       fp   => dataptr(F,n)
       up   => dataptr(U,n)
       u13p => dataptr(U13,n)

       lo = lwb(get_box(U13,n))
       hi = upb(get_box(U13,n))

       do m = 1, nc
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   u13p(i,j,k,m) = up(i,j,k,m) + dt * (dp(i,j,k,m) + fp(i,j,k,m))
                end do
             end do
          end do
       end do
    end do
    !
    ! Sync up ghost cells and periodic boundary cells in U13.
    !
    call multifab_fill_boundary(U13)
    !
    ! Calculate D at time N+1/3.
    !
    do n=1,nboxes(D)
       if ( remote(D,n) ) cycle

       dp => dataptr(D,n)

       lo = lwb(get_box(D,n))
       hi = upb(get_box(D,n))
       !
       ! Use U13.
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
       ! Use U13.
       !
       fp = 0.0d0
    end do
    !
    ! Calculate U at time N+2/3.
    !
    do n=1,nboxes(U23)
       if ( remote(U23,n) ) cycle

       dp   => dataptr(D,n)
       fp   => dataptr(F,n)
       up   => dataptr(U,n)
       u13p => dataptr(U13,n)
       u23p => dataptr(U23,n)

       lo = lwb(get_box(U23,n))
       hi = upb(get_box(U23,n))

       do m = 1, nc
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   u23p(i,j,k,m) = ThreeQuarters * up(i,j,k,m) + &
                        OneQuarter * (u13p(i,j,k,m) + dt * (dp(i,j,k,m) + fp(i,j,k,m)))
                end do
             end do
          end do
       end do
    end do
    !
    ! Sync up ghost cells and periodic boundary cells in U23.
    !
    call multifab_fill_boundary(U13)
    !
    ! Calculate D at time N+2/3.
    !
    do n=1,nboxes(D)
       if ( remote(D,n) ) cycle

       dp => dataptr(D,n)

       lo = lwb(get_box(D,n))
       hi = upb(get_box(D,n))
       !
       ! Use U23.
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
       ! Use U23.
       !
       fp = 0.0d0
    end do
    !
    ! Calculate U at time N+1.
    !
    do n=1,nboxes(U23)
       if ( remote(U23,n) ) cycle

       dp   => dataptr(D,n)
       fp   => dataptr(F,n)
       up   => dataptr(U,n)
       u23p => dataptr(U23,n)

       lo = lwb(get_box(U23,n))
       hi = upb(get_box(U23,n))

       do m = 1, nc
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   up(i,j,k,m) = OneThird * up(i,j,k,m) + &
                        TwoThirds * (u23p(i,j,k,m) + dt * (dp(i,j,k,m) + fp(i,j,k,m)))
                end do
             end do
          end do
       end do
    end do

    call destroy(U23)
    call destroy(U13)
    call destroy(F)
    call destroy(D)

  end subroutine advance

end module advance_module

