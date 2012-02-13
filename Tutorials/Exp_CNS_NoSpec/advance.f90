 module advance_module

  use bl_error_module
  use multifab_module

  implicit none

  private

  integer, parameter :: irho = 1
  integer, parameter :: imx  = 2
  integer, parameter :: imy  = 3
  integer, parameter :: imz  = 4
  integer, parameter :: iene = 5

  integer, parameter :: qrho  = 1
  integer, parameter :: qu    = 2
  integer, parameter :: qv    = 3
  integer, parameter :: qw    = 4
  integer, parameter :: qpres = 5
  integer, parameter :: qT    = 5

  public :: advance

contains
  
  subroutine advance (U,dt,dx)

    type(multifab),   intent(inout) :: U
    double precision, intent(in   ) :: dt

    double precision, parameter :: OneThird      = 1.d0/3.d0
    double precision, parameter :: TwoThirds     = 2.d0/3.d0
    double precision, parameter :: OneQuarter    = 1.d0/4.d0
    double precision, parameter :: ThreeQuarters = 3.d0/4.d0

    integer          :: lo(U%dim), hi(U%dim)
    integer          :: i, j, k, m, n, nc, ng
    double precision :: dx(U%dim), courno, courno_proc
    type(layout)     :: la
    type(multifab)   :: D, F, Unew, Q

    double precision, pointer, dimension(:,:,:,:) :: up, dp, fp, unp, qp

    nc = ncomp(U)
    ng = nghost(U)
    la = get_layout(U)
    !
    ! Sync U prior to calculating D & F.
    !
    call multifab_fill_boundary(U)

    call multifab_build(D,    la, nc, 0)
    call multifab_build(F,    la, nc, 0)
    call multifab_build(Q,    la, nc, ng+1)
    call multifab_build(Unew, la, nc, ng)
    !
    ! Calculate primitive variables based on U.
    !
    courno_proc = 1.0e-50

    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

       up => dataptr(U,n)
       qp => dataptr(Q,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim(lo,hi,up,qp,courno_proc,dx,ng)
    end do

    call parallel_reduce(courno, courno_proc, MPI_MAX)

    if ( parallel_IOProcessor() ) then
       print*, "courno = ", courno
    end if
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

       up => dataptr(U,n)
       qp => dataptr(Q,n)
       fp => dataptr(F,n)

       lo = lwb(get_box(F,n))
       hi = upb(get_box(F,n))

       call hypterm(lo,hi,ng,dx,up,qp,fp)
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
    ! Calculate primitive variables based on U^1/3.
    !
    courno_proc = 1.0e-50

    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

       up => dataptr(Unew,n)
       qp => dataptr(Q,   n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim(lo,hi,up,qp,courno_proc,dx,ng)
    end do

    !call parallel_reduce(courno, courno_proc, MPI_MAX)

    !if ( parallel_IOProcessor() ) then
    !   print*, "courno = ", courno
    !end if

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

       up => dataptr(Unew,n)
       qp => dataptr(Q,   n)
       fp => dataptr(F,   n)

       lo = lwb(get_box(F,n))
       hi = upb(get_box(F,n))

       call hypterm(lo,hi,ng,dx,up,qp,fp)
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
    ! Calculate primitive variables based on U^2/3.
    !
    courno_proc = 1.0e-50

    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

       up => dataptr(Unew,n)
       qp => dataptr(Q,   n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim(lo,hi,up,qp,courno_proc,dx,ng)
    end do

    !call parallel_reduce(courno, courno_proc, MPI_MAX)

    !if ( parallel_IOProcessor() ) then
    !   print*, "courno = ", courno
    !end if

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

       up => dataptr(Unew,n)
       qp => dataptr(Q,   n)
       fp => dataptr(F,   n)

       lo = lwb(get_box(F,n))
       hi = upb(get_box(F,n))

       call hypterm(lo,hi,ng,dx,up,qp,fp)
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
    call destroy(Q)
    call destroy(F)
    call destroy(D)

  end subroutine advance

  subroutine ctoprim (lo,hi,u,q,courno,dx,ng)

    integer,          intent(in   ) :: lo(3), hi(3), ng
    double precision, intent(in   ) :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(out  ) :: q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,6)
    double precision, intent(inout) :: courno

    integer          :: i, j, k
    double precision :: c, eint, courx, coury, courz, courmx, courmy, courmz

    double precision, parameter :: GAMMA = 1.4d0
    double precision, parameter :: CV    = 8.3333333333d6

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3)-ng,hi(3)+ng
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng

             q(i,j,k,1) = u(i,j,k,1)
             q(i,j,k,2) = u(i,j,k,2)/u(i,j,k,1)
             q(i,j,k,3) = u(i,j,k,3)/u(i,j,k,1)
             q(i,j,k,4) = u(i,j,k,4)/u(i,j,k,1)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! Get gamc, p, T, c, csml using q state
    !$OMP PARALLEL DO PRIVATE(i,j,k,pt_index)
    do k = lo(3)-ng,hi(3)+ng
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng

             eint =  u(i,j,k,5)/u(i,j,k,1)      & 
                  - 0.5d0*( q(i,j,k,2)**2 + q(i,j,k,3)**2 + q(i,j,k,4) **2)
             q(i,j,k,5) = (GAMMA-1.d0)*eint*u(i,j,k,1)
             q(i,j,k,6) = eint/CV

          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,courx,coury,courz) REDUCTION(max:courmx,courmy,courmz)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             c     = sqrt(GAMMA*q(i,j,k,5)/q(i,j,k,1))
             courx = ( c+abs(q(i,j,k,2)) ) * dx(1)
             coury = ( c+abs(q(i,j,k,3)) ) * dx(2)
             courz = ( c+abs(q(i,j,k,4)) ) * dx(3)

             courmx = max( courmx, courx )
             courmy = max( courmy, coury )
             courmz = max( courmz, courz )

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    !
    ! Compute running max of Courant number over grids.
    !
    courno = max( courmx, courmy, courmz , courno )

  end subroutine ctoprim

  subroutine hypterm (lo,hi,ng,dx,cons,q,flux)

    double precision, parameter :: ALP =  0.8d0
    double precision, parameter :: BET = -0.2d0
    double precision, parameter :: GAM =  4.d0/105.d0
    double precision, parameter :: DEL = -1.d0/280.d0

    ! inputs:  lo,hi,ng,cons,q
    ! cons --  rho, rho u , rho v, rho w, rho E
    ! cons --  rho, u ,  v,  w, p, T
    ! outputs: flux

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,5)
    double precision, intent(in ) ::    q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,6)
    double precision, intent(out) :: flux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)

    integer          :: i,j,k
    double precision :: unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = q(i+1,j,k,qu)
             unp2 = q(i+2,j,k,qu)
             unp3 = q(i+3,j,k,qu)
             unp4 = q(i+4,j,k,qu)

             unm1 = q(i-1,j,k,qu)
             unm2 = q(i-2,j,k,qu)
             unm3 = q(i-3,j,k,qu)
             unm4 = q(i-4,j,k,qu)

             flux(i,j,k,irho)= -(ALP*(cons(i+1,j,k,imx)-cons(i-1,j,k,imx)) &
                  + BET*(cons(i+2,j,k,imx)-cons(i-2,j,k,imx))              &
                  + GAM*(cons(i+3,j,k,imx)-cons(i-3,j,k,imx))              &
                  + DEL*(cons(i+4,j,k,imx)-cons(i-4,j,k,imx)))/dx(1)

             flux(i,j,k,imx)= -(ALP*(cons(i+1,j,k,imx)*unp1-cons(i-1,j,k,imx)*unm1 &
                  + (q(i+1,j,k,qpres)-q(i-1,j,k,qpres)))                           &
                  + BET*(cons(i+2,j,k,imx)*unp2-cons(i-2,j,k,imx)*unm2             &
                  + (q(i+2,j,k,qpres)-q(i-2,j,k,qpres)))                           &
                  + GAM*(cons(i+3,j,k,imx)*unp3-cons(i-3,j,k,imx)*unm3             &
                  + (q(i+3,j,k,qpres)-q(i-3,j,k,qpres)))                           &
                  + DEL*(cons(i+4,j,k,imx)*unp4-cons(i-4,j,k,imx)*unm4             &
                  + (q(i+4,j,k,qpres)-q(i-4,j,k,qpres)))) / dx(1)

             flux(i,j,k,imy)= -                                         &
                  (ALP*(cons(i+1,j,k,imy)*unp1-cons(i-1,j,k,imy)*unm1)  &
                  + BET*(cons(i+2,j,k,imy)*unp2-cons(i-2,j,k,imy)*unm2) &
                  + GAM*(cons(i+3,j,k,imy)*unp3-cons(i-3,j,k,imy)*unm3) &
                  + DEL*(cons(i+4,j,k,imy)*unp4-cons(i-4,j,k,imy)*unm4))/dx(1)

             flux(i,j,k,imz)= -                                         &
                  (ALP*(cons(i+1,j,k,imz)*unp1-cons(i-1,j,k,imz)*unm1)  &
                  + BET*(cons(i+2,j,k,imz)*unp2-cons(i-2,j,k,imz)*unm2) &
                  + GAM*(cons(i+3,j,k,imz)*unp3-cons(i-3,j,k,imz)*unm3) &
                  + DEL*(cons(i+4,j,k,imz)*unp4-cons(i-4,j,k,imz)*unm4))/dx(1)

             flux(i,j,k,iene)= -                                         &
                  (ALP*(cons(i+1,j,k,iene)*unp1-cons(i-1,j,k,iene)*unm1  &
                  + (q(i+1,j,k,qpres)*unp1-q(i-1,j,k,qpres)*unm1))       &
                  + BET*(cons(i+2,j,k,iene)*unp2-cons(i-2,j,k,iene)*unm2 &
                  + (q(i+2,j,k,qpres)*unp2-q(i-2,j,k,qpres)*unm2))       &
                  + GAM*(cons(i+3,j,k,iene)*unp3-cons(i-3,j,k,iene)*unm3 &
                  + (q(i+3,j,k,qpres)*unp3-q(i-3,j,k,qpres)*unm3))       &
                  + DEL*(cons(i+4,j,k,iene)*unp4-cons(i-4,j,k,iene)*unm4 &
                  + (q(i+4,j,k,qpres)*unp4-q(i-4,j,k,qpres)*unm4))) / dx(1) 

          enddo
       enddo
    enddo

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = q(i,j+1,k,qv)
             unp2 = q(i,j+2,k,qv)
             unp3 = q(i,j+3,k,qv)
             unp4 = q(i,j+4,k,qv)

             unm1 = q(i,j-1,k,qv)
             unm2 = q(i,j-2,k,qv)
             unm3 = q(i,j-3,k,qv)
             unm4 = q(i,j-4,k,qv)

             flux(i,j,k,irho)=flux(i,j,k,irho) -              &
                  (ALP*(cons(i,j+1,k,imy)-cons(i,j-1,k,imy))  &
                  + BET*(cons(i,j+2,k,imy)-cons(i,j-2,k,imy)) &
                  + GAM*(cons(i,j+3,k,imy)-cons(i,j-3,k,imy)) &
                  + DEL*(cons(i,j+4,k,imy)-cons(i,j-4,k,imy)))/dx(2)

             flux(i,j,k,imx)=flux(i,j,k,imx) -                          &
                  (ALP*(cons(i,j+1,k,imx)*unp1-cons(i,j-1,k,imx)*unm1)  &
                  + BET*(cons(i,j+2,k,imx)*unp2-cons(i,j-2,k,imx)*unm2) &
                  + GAM*(cons(i,j+3,k,imx)*unp3-cons(i,j-3,k,imx)*unm3) &
                  + DEL*(cons(i,j+4,k,imx)*unp4-cons(i,j-4,k,imx)*unm4))/dx(2)

             flux(i,j,k,imy)=flux(i,j,k,imy) -                         &
                  (ALP*(cons(i,j+1,k,imy)*unp1-cons(i,j-1,k,imy)*unm1  &
                  + (q(i,j+1,k,qpres)-q(i,j-1,k,qpres)))               &
                  + BET*(cons(i,j+2,k,imy)*unp2-cons(i,j-2,k,imy)*unm2 &
                  + (q(i,j+2,k,qpres)-q(i,j-2,k,qpres)))               &
                  + GAM*(cons(i,j+3,k,imy)*unp3-cons(i,j-3,k,imy)*unm3 &
                  + (q(i,j+3,k,qpres)-q(i,j-3,k,qpres)))               &
                  + DEL*(cons(i,j+4,k,imy)*unp4-cons(i,j-4,k,imy)*unm4 &
                  + (q(i,j+4,k,qpres)-q(i,j-4,k,qpres)))) / dx(2)

             flux(i,j,k,imz)=flux(i,j,k,imz) -                          &
                  (ALP*(cons(i,j+1,k,imz)*unp1-cons(i,j-1,k,imz)*unm1)  &
                  + BET*(cons(i,j+2,k,imz)*unp2-cons(i,j-2,k,imz)*unm2) &
                  + GAM*(cons(i,j+3,k,imz)*unp3-cons(i,j-3,k,imz)*unm3) &
                  + DEL*(cons(i,j+4,k,imz)*unp4-cons(i,j-4,k,imz)*unm4))/dx(2)

             flux(i,j,k,iene)=flux(i,j,k,iene) -                         &
                  (ALP*(cons(i,j+1,k,iene)*unp1-cons(i,j-1,k,iene)*unm1  &
                  + (q(i,j+1,k,qpres)*unp1-q(i,j-1,k,qpres)*unm1))       &
                  + BET*(cons(i,j+2,k,iene)*unp2-cons(i,j-2,k,iene)*unm2 &
                  + (q(i,j+2,k,qpres)*unp2-q(i,j-2,k,qpres)*unm2))       &
                  + GAM*(cons(i,j+3,k,iene)*unp3-cons(i,j-3,k,iene)*unm3 &
                  + (q(i,j+3,k,qpres)*unp3-q(i,j-3,k,qpres)*unm3))       &
                  + DEL*(cons(i,j+4,k,iene)*unp4-cons(i,j-4,k,iene)*unm4 &
                  + (q(i,j+4,k,qpres)*unp4-q(i,j-4,k,qpres)*unm4))) / dx(2)

          enddo
       enddo
    enddo

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = q(i,j,k+1,qw)
             unp2 = q(i,j,k+2,qw)
             unp3 = q(i,j,k+3,qw)
             unp4 = q(i,j,k+4,qw)

             unm1 = q(i,j,k-1,qw)
             unm2 = q(i,j,k-2,qw)
             unm3 = q(i,j,k-3,qw)
             unm4 = q(i,j,k-4,qw)

             flux(i,j,k,irho)=flux(i,j,k,irho) -              &
                  (ALP*(cons(i,j,k+1,imz)-cons(i,j,k-1,imz))  &
                  + BET*(cons(i,j,k+2,imz)-cons(i,j,k-2,imz)) &
                  + GAM*(cons(i,j,k+3,imz)-cons(i,j,k-3,imz)) &
                  + DEL*(cons(i,j,k+4,imz)-cons(i,j,k-4,imz)))/dx(3)

             flux(i,j,k,imx)=flux(i,j,k,imx) -                          &
                  (ALP*(cons(i,j,k+1,imx)*unp1-cons(i,j,k-1,imx)*unm1)  &
                  + BET*(cons(i,j,k+2,imx)*unp2-cons(i,j,k-2,imx)*unm2) &
                  + GAM*(cons(i,j,k+3,imx)*unp3-cons(i,j,k-3,imx)*unm3) &
                  + DEL*(cons(i,j,k+4,imx)*unp4-cons(i,j,k-4,imx)*unm4))/dx(3)

             flux(i,j,k,imy)=flux(i,j,k,imy) -                          &
                  (ALP*(cons(i,j,k+1,imy)*unp1-cons(i,j,k-1,imy)*unm1)  &
                  + BET*(cons(i,j,k+2,imy)*unp2-cons(i,j,k-2,imy)*unm2) &
                  + GAM*(cons(i,j,k+3,imy)*unp3-cons(i,j,k-3,imy)*unm3) &
                  + DEL*(cons(i,j,k+4,imy)*unp4-cons(i,j,k-4,imy)*unm4))/dx(3)

             flux(i,j,k,imz)=flux(i,j,k,imz) -                         &
                  (ALP*(cons(i,j,k+1,imz)*unp1-cons(i,j,k-1,imz)*unm1  &
                  + (q(i,j,k+1,qpres)-q(i,j,k-1,qpres)))               &
                  + BET*(cons(i,j,k+2,imz)*unp2-cons(i,j,k-2,imz)*unm2 &
                  + (q(i,j,k+2,qpres)-q(i,j,k-2,qpres)))               &
                  + GAM*(cons(i,j,k+3,imz)*unp3-cons(i,j,k-3,imz)*unm3 &
                  + (q(i,j,k+3,qpres)-q(i,j,k-3,qpres)))               &
                  + DEL*(cons(i,j,k+4,imz)*unp4-cons(i,j,k-4,imz)*unm4 &
                  + (q(i,j,k+4,qpres)-q(i,j,k-4,qpres)))) / dx(3)

             flux(i,j,k,iene)=flux(i,j,k,iene) -                         &
                  (ALP*(cons(i,j,k+1,iene)*unp1-cons(i,j,k-1,iene)*unm1  &
                  + (q(i,j,k+1,qpres)*unp1-q(i,j,k-1,qpres)*unm1))       &
                  + BET*(cons(i,j,k+2,iene)*unp2-cons(i,j,k-2,iene)*unm2 &
                  + (q(i,j,k+2,qpres)*unp2-q(i,j,k-2,qpres)*unm2))       &
                  + GAM*(cons(i,j,k+3,iene)*unp3-cons(i,j,k-3,iene)*unm3 &
                  + (q(i,j,k+3,qpres)*unp3-q(i,j,k-3,qpres)*unm3))       &
                  + DEL*(cons(i,j,k+4,iene)*unp4-cons(i,j,k-4,iene)*unm4 &
                  + (q(i,j,k+4,qpres)*unp4-q(i,j,k-4,qpres)*unm4))) / dx(3)

          enddo
       enddo
    enddo

  end subroutine hypterm

end module advance_module

