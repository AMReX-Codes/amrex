 module advance_module

  use bl_error_module
  use multifab_module
  use bl_prof_module

  implicit none

  private
  !
  ! These index constants are shared with the initial data routine.
  !
  integer, parameter, public :: irho = 1
  integer, parameter, public :: imx  = 2
  integer, parameter, public :: imy  = 3
  integer, parameter, public :: imz  = 4
  integer, parameter, public :: iene = 5

  integer, parameter :: qu    = 2
  integer, parameter :: qv    = 3
  integer, parameter :: qw    = 4
  integer, parameter :: qpres = 5

  double precision, parameter :: ALP =  0.8d0
  double precision, parameter :: BET = -0.2d0
  double precision, parameter :: GAM =  4.d0/105.d0
  double precision, parameter :: DEL = -1.d0/280.d0

  public :: advance

contains


  subroutine advance (U,dt,dx,cfl,eta,alam)

    use bl_prof_module

    type(multifab),   intent(inout) :: U
    double precision, intent(out  ) :: dt
    double precision, intent(in   ) :: dx(U%dim), cfl, eta, alam

    integer          :: lo(U%dim), hi(U%dim), i, j, k, m, n, nc, ng
    double precision :: courno, courno_proc
    type(layout)     :: la
    type(multifab)   :: D, F, Unew, Q

    double precision, pointer, dimension(:,:,:,:) :: up, dp, fp, unp, qp
    !
    ! Some arithmetic constants.
    !
    double precision, parameter :: OneThird      = 1.d0/3.d0
    double precision, parameter :: TwoThirds     = 2.d0/3.d0
    double precision, parameter :: OneQuarter    = 1.d0/4.d0
    double precision, parameter :: ThreeQuarters = 3.d0/4.d0

    type(bl_prof_timer), save :: bpt_advance

    call build(bpt_advance, "bpt_advance")

    nc = ncomp(U)
    ng = nghost(U)
    la = get_layout(U)
    !
    ! Sync U prior to calculating D & F.
    !

    call multifab_fill_boundary(U)

    call multifab_build(D,    la, nc,   0)
    call multifab_build(F,    la, nc,   0)
    call multifab_build(Q,    la, nc+1, ng)
    call multifab_build(Unew, la, nc,   ng)
    !
    ! Calculate primitive variables based on U.
    !
    ! Also calculate courno so we can set "dt".
    !
    courno_proc = 1.0d-50


    do n=1,nfabs(Q)

       up => dataptr(U,n)
       qp => dataptr(Q,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim(lo,hi,up,qp,dx,ng,courno=courno_proc)
    end do

    call parallel_reduce(courno, courno_proc, MPI_MAX)

    dt = cfl / courno

    if ( parallel_IOProcessor() ) then
       print*, "dt,courno", dt, courno
    end if
    !
    ! Calculate D at time N.
    !
    do n=1,nfabs(D)

       qp => dataptr(Q,n)
       dp => dataptr(D,n)

       lo = lwb(get_box(D,n))
       hi = upb(get_box(D,n))

       call diffterm(lo,hi,ng,dx,qp,dp,ETA,ALAM)
    end do

    !
    ! Calculate F at time N.
    !
    do n=1,nfabs(F)

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
    do n=1,nfabs(U)

       dp  => dataptr(D,   n)
       fp  => dataptr(F,   n)
       up  => dataptr(U,   n)
       unp => dataptr(Unew,n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       do m = 1, nc
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   unp(i,j,k,m) = up(i,j,k,m) + dt * (dp(i,j,k,m) + fp(i,j,k,m))
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do
    !
    ! Sync U^1/3 prior to calculating D & F.
    !
    call multifab_fill_boundary(Unew)
    !
    ! Calculate primitive variables based on U^1/3.
    !
    do n=1,nfabs(Q)

       up => dataptr(Unew,n)
       qp => dataptr(Q,   n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim(lo,hi,up,qp,dx,ng)
    end do
    !
    ! Calculate D at time N+1/3.
    !
    do n=1,nfabs(D)

       qp => dataptr(Q,n)
       dp => dataptr(D,n)

       lo = lwb(get_box(D,n))
       hi = upb(get_box(D,n))

       call diffterm(lo,hi,ng,dx,qp,dp,ETA,ALAM)
    end do
    !
    ! Calculate F at time N+1/3.
    !
    do n=1,nfabs(F)

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
    do n=1,nfabs(U)

       dp  => dataptr(D,   n)
       fp  => dataptr(F,   n)
       up  => dataptr(U,   n)
       unp => dataptr(Unew,n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       do m = 1, nc
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   unp(i,j,k,m) = ThreeQuarters * up(i,j,k,m) + &
                        OneQuarter * (unp(i,j,k,m) + dt * (dp(i,j,k,m) + fp(i,j,k,m)))
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do
    !
    ! Sync U^2/3 prior to calculating D & F.
    !
    call multifab_fill_boundary(Unew)
    !
    ! Calculate primitive variables based on U^2/3.
    !
    do n=1,nfabs(Q)

       up => dataptr(Unew,n)
       qp => dataptr(Q,   n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim(lo,hi,up,qp,dx,ng)
    end do
    !
    ! Calculate D at time N+2/3.
    !
    do n=1,nfabs(D)

       qp => dataptr(Q,n)
       dp => dataptr(D,n)

       lo = lwb(get_box(D,n))
       hi = upb(get_box(D,n))

       call diffterm(lo,hi,ng,dx,qp,dp,ETA,ALAM)
    end do
    !
    ! Calculate F at time N+2/3.
    !
    do n=1,nfabs(F)

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
    do n=1,nfabs(U)

       dp  => dataptr(D,   n)
       fp  => dataptr(F,   n)
       up  => dataptr(U,   n)
       unp => dataptr(Unew,n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       do m = 1, nc
          !$OMP PARALLEL DO PRIVATE(i,j,k)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   up(i,j,k,m) = OneThird * up(i,j,k,m) + &
                        TwoThirds * (unp(i,j,k,m) + dt * (dp(i,j,k,m) + fp(i,j,k,m)))
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do

    call destroy(Unew)
    call destroy(Q)
    call destroy(F)
    call destroy(D)

    call destroy(bpt_advance)

  end subroutine advance



  subroutine ctoprim (lo,hi,u,q,dx,ng,courno)

    use bl_prof_module

    integer,          intent(in ) :: lo(3), hi(3), ng
    double precision, intent(in ) :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)
    double precision, intent(in ) :: dx(3)
    double precision, intent(out) :: q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,6)

    double precision, intent(inout), optional :: courno

    integer          :: i, j, k
    double precision :: c, eint, courx, coury, courz, courmx, courmy, courmz, rhoinv
    double precision :: dx1inv, dx2inv, dx3inv, CVinv

    double precision, parameter :: GAMMA = 1.4d0
    double precision, parameter :: CV    = 8.3333333333d6

    type(bl_prof_timer), save :: bpt_ctoprim, bpt_ctoprim_loop1, bpt_ctoprim_loop2

    call build(bpt_ctoprim, "bpt_ctoprim")

    CVinv = 1.0d0 / CV

    call build(bpt_ctoprim_loop1, "bpt_ctoprim_loop1")
    !$OMP PARALLEL DO PRIVATE(i,j,k,eint,rhoinv)
    do k = lo(3)-ng,hi(3)+ng
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng

             rhoinv     = 1.0d0/u(i,j,k,1)
             q(i,j,k,1) = u(i,j,k,1)
             q(i,j,k,2) = u(i,j,k,2)*rhoinv
             q(i,j,k,3) = u(i,j,k,3)*rhoinv
             q(i,j,k,4) = u(i,j,k,4)*rhoinv

             eint = u(i,j,k,5)*rhoinv - 0.5d0*(q(i,j,k,2)**2 + q(i,j,k,3)**2 + q(i,j,k,4)**2)

             q(i,j,k,5) = (GAMMA-1.d0)*eint*u(i,j,k,1)
             q(i,j,k,6) = eint * CVinv

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_ctoprim_loop1)

    if ( present(courno) ) then

       courmx = -Huge(courmx)
       courmy = -Huge(courmy)
       courmz = -Huge(courmz)

       dx1inv = 1.0d0 / dx(1)
       dx2inv = 1.0d0 / dx(2)
       dx3inv = 1.0d0 / dx(3)

       call build(bpt_ctoprim_loop2, "bpt_ctoprim_loop2")
       !$OMP PARALLEL DO PRIVATE(i,j,k,c,courx,coury,courz) REDUCTION(max:courmx,courmy,courmz)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                c     = sqrt(GAMMA*q(i,j,k,5)/q(i,j,k,1))
                courx = ( c+abs(q(i,j,k,2)) ) * dx1inv
                coury = ( c+abs(q(i,j,k,3)) ) * dx2inv
                courz = ( c+abs(q(i,j,k,4)) ) * dx3inv

                courmx = max( courmx, courx )
                courmy = max( courmy, coury )
                courmz = max( courmz, courz )

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
       call destroy(bpt_ctoprim_loop2)
       !
       ! Compute running max of Courant number over grids.
       !
       courno = max( courmx, courmy, courmz , courno )

    end if

    call destroy(bpt_ctoprim)

  end subroutine ctoprim



  subroutine hypterm (lo,hi,ng,dx,cons,q,flux)

    use bl_prof_module

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,5)
    double precision, intent(in ) ::    q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,6)
    double precision, intent(out) :: flux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)

    integer          :: i,j,k
    double precision :: unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4
    double precision :: dxinv(3)

    type(bl_prof_timer), save :: bpt_hypterm, bpt_hypterm_loop1, bpt_hypterm_loop2, bpt_hypterm_loop3

    call build(bpt_hypterm, "bpt_hypterm")


    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do


    call build(bpt_hypterm_loop1, "bpt_hypterm_loop1")
    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
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

             flux(i,j,k,irho)= - &
                   (ALP*(cons(i+1,j,k,imx)-cons(i-1,j,k,imx)) &
                  + BET*(cons(i+2,j,k,imx)-cons(i-2,j,k,imx)) &
                  + GAM*(cons(i+3,j,k,imx)-cons(i-3,j,k,imx)) &
                  + DEL*(cons(i+4,j,k,imx)-cons(i-4,j,k,imx)))*dxinv(1)

             flux(i,j,k,imx)= - &
                   (ALP*(cons(i+1,j,k,imx)*unp1-cons(i-1,j,k,imx)*unm1 &
                  + (q(i+1,j,k,qpres)-q(i-1,j,k,qpres)))               &
                  + BET*(cons(i+2,j,k,imx)*unp2-cons(i-2,j,k,imx)*unm2 &
                  + (q(i+2,j,k,qpres)-q(i-2,j,k,qpres)))               &
                  + GAM*(cons(i+3,j,k,imx)*unp3-cons(i-3,j,k,imx)*unm3 &
                  + (q(i+3,j,k,qpres)-q(i-3,j,k,qpres)))               &
                  + DEL*(cons(i+4,j,k,imx)*unp4-cons(i-4,j,k,imx)*unm4 &
                  + (q(i+4,j,k,qpres)-q(i-4,j,k,qpres))))*dxinv(1)

             flux(i,j,k,imy)= - &
                   (ALP*(cons(i+1,j,k,imy)*unp1-cons(i-1,j,k,imy)*unm1) &
                  + BET*(cons(i+2,j,k,imy)*unp2-cons(i-2,j,k,imy)*unm2) &
                  + GAM*(cons(i+3,j,k,imy)*unp3-cons(i-3,j,k,imy)*unm3) &
                  + DEL*(cons(i+4,j,k,imy)*unp4-cons(i-4,j,k,imy)*unm4))*dxinv(1)

             flux(i,j,k,imz)= - &
                   (ALP*(cons(i+1,j,k,imz)*unp1-cons(i-1,j,k,imz)*unm1) &
                  + BET*(cons(i+2,j,k,imz)*unp2-cons(i-2,j,k,imz)*unm2) &
                  + GAM*(cons(i+3,j,k,imz)*unp3-cons(i-3,j,k,imz)*unm3) &
                  + DEL*(cons(i+4,j,k,imz)*unp4-cons(i-4,j,k,imz)*unm4))*dxinv(1)

             flux(i,j,k,iene)= - &
                   (ALP*(cons(i+1,j,k,iene)*unp1-cons(i-1,j,k,iene)*unm1 &
                  + (q(i+1,j,k,qpres)*unp1-q(i-1,j,k,qpres)*unm1))       &
                  + BET*(cons(i+2,j,k,iene)*unp2-cons(i-2,j,k,iene)*unm2 &
                  + (q(i+2,j,k,qpres)*unp2-q(i-2,j,k,qpres)*unm2))       &
                  + GAM*(cons(i+3,j,k,iene)*unp3-cons(i-3,j,k,iene)*unm3 &
                  + (q(i+3,j,k,qpres)*unp3-q(i-3,j,k,qpres)*unm3))       &
                  + DEL*(cons(i+4,j,k,iene)*unp4-cons(i-4,j,k,iene)*unm4 &
                  + (q(i+4,j,k,qpres)*unp4-q(i-4,j,k,qpres)*unm4)))*dxinv(1) 
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_hypterm_loop1)

    call build(bpt_hypterm_loop2, "bpt_hypterm_loop2")
    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
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

             flux(i,j,k,irho)=flux(i,j,k,irho) - &
                   (ALP*(cons(i,j+1,k,imy)-cons(i,j-1,k,imy)) &
                  + BET*(cons(i,j+2,k,imy)-cons(i,j-2,k,imy)) &
                  + GAM*(cons(i,j+3,k,imy)-cons(i,j-3,k,imy)) &
                  + DEL*(cons(i,j+4,k,imy)-cons(i,j-4,k,imy)))*dxinv(2)

             flux(i,j,k,imx)=flux(i,j,k,imx) - &
                   (ALP*(cons(i,j+1,k,imx)*unp1-cons(i,j-1,k,imx)*unm1) &
                  + BET*(cons(i,j+2,k,imx)*unp2-cons(i,j-2,k,imx)*unm2) &
                  + GAM*(cons(i,j+3,k,imx)*unp3-cons(i,j-3,k,imx)*unm3) &
                  + DEL*(cons(i,j+4,k,imx)*unp4-cons(i,j-4,k,imx)*unm4))*dxinv(2)

             flux(i,j,k,imy)=flux(i,j,k,imy) - &
                   (ALP*(cons(i,j+1,k,imy)*unp1-cons(i,j-1,k,imy)*unm1 &
                  + (q(i,j+1,k,qpres)-q(i,j-1,k,qpres)))               &
                  + BET*(cons(i,j+2,k,imy)*unp2-cons(i,j-2,k,imy)*unm2 &
                  + (q(i,j+2,k,qpres)-q(i,j-2,k,qpres)))               &
                  + GAM*(cons(i,j+3,k,imy)*unp3-cons(i,j-3,k,imy)*unm3 &
                  + (q(i,j+3,k,qpres)-q(i,j-3,k,qpres)))               &
                  + DEL*(cons(i,j+4,k,imy)*unp4-cons(i,j-4,k,imy)*unm4 &
                  + (q(i,j+4,k,qpres)-q(i,j-4,k,qpres))))*dxinv(2)

             flux(i,j,k,imz)=flux(i,j,k,imz) - &
                   (ALP*(cons(i,j+1,k,imz)*unp1-cons(i,j-1,k,imz)*unm1) &
                  + BET*(cons(i,j+2,k,imz)*unp2-cons(i,j-2,k,imz)*unm2) &
                  + GAM*(cons(i,j+3,k,imz)*unp3-cons(i,j-3,k,imz)*unm3) &
                  + DEL*(cons(i,j+4,k,imz)*unp4-cons(i,j-4,k,imz)*unm4))*dxinv(2)

             flux(i,j,k,iene)=flux(i,j,k,iene) - &
                   (ALP*(cons(i,j+1,k,iene)*unp1-cons(i,j-1,k,iene)*unm1 &
                  + (q(i,j+1,k,qpres)*unp1-q(i,j-1,k,qpres)*unm1))       &
                  + BET*(cons(i,j+2,k,iene)*unp2-cons(i,j-2,k,iene)*unm2 &
                  + (q(i,j+2,k,qpres)*unp2-q(i,j-2,k,qpres)*unm2))       &
                  + GAM*(cons(i,j+3,k,iene)*unp3-cons(i,j-3,k,iene)*unm3 &
                  + (q(i,j+3,k,qpres)*unp3-q(i,j-3,k,qpres)*unm3))       &
                  + DEL*(cons(i,j+4,k,iene)*unp4-cons(i,j-4,k,iene)*unm4 &
                  + (q(i,j+4,k,qpres)*unp4-q(i,j-4,k,qpres)*unm4)))*dxinv(2)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_hypterm_loop2)

    call build(bpt_hypterm_loop3, "bpt_hypterm_loop3")
    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
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

             flux(i,j,k,irho)=flux(i,j,k,irho) - &
                   (ALP*(cons(i,j,k+1,imz)-cons(i,j,k-1,imz)) &
                  + BET*(cons(i,j,k+2,imz)-cons(i,j,k-2,imz)) &
                  + GAM*(cons(i,j,k+3,imz)-cons(i,j,k-3,imz)) &
                  + DEL*(cons(i,j,k+4,imz)-cons(i,j,k-4,imz)))*dxinv(3)

             flux(i,j,k,imx)=flux(i,j,k,imx) - &
                   (ALP*(cons(i,j,k+1,imx)*unp1-cons(i,j,k-1,imx)*unm1) &
                  + BET*(cons(i,j,k+2,imx)*unp2-cons(i,j,k-2,imx)*unm2) &
                  + GAM*(cons(i,j,k+3,imx)*unp3-cons(i,j,k-3,imx)*unm3) &
                  + DEL*(cons(i,j,k+4,imx)*unp4-cons(i,j,k-4,imx)*unm4))*dxinv(3)

             flux(i,j,k,imy)=flux(i,j,k,imy) - &
                   (ALP*(cons(i,j,k+1,imy)*unp1-cons(i,j,k-1,imy)*unm1) &
                  + BET*(cons(i,j,k+2,imy)*unp2-cons(i,j,k-2,imy)*unm2) &
                  + GAM*(cons(i,j,k+3,imy)*unp3-cons(i,j,k-3,imy)*unm3) &
                  + DEL*(cons(i,j,k+4,imy)*unp4-cons(i,j,k-4,imy)*unm4))*dxinv(3)

             flux(i,j,k,imz)=flux(i,j,k,imz) - &
                   (ALP*(cons(i,j,k+1,imz)*unp1-cons(i,j,k-1,imz)*unm1 &
                  + (q(i,j,k+1,qpres)-q(i,j,k-1,qpres)))               &
                  + BET*(cons(i,j,k+2,imz)*unp2-cons(i,j,k-2,imz)*unm2 &
                  + (q(i,j,k+2,qpres)-q(i,j,k-2,qpres)))               &
                  + GAM*(cons(i,j,k+3,imz)*unp3-cons(i,j,k-3,imz)*unm3 &
                  + (q(i,j,k+3,qpres)-q(i,j,k-3,qpres)))               &
                  + DEL*(cons(i,j,k+4,imz)*unp4-cons(i,j,k-4,imz)*unm4 &
                  + (q(i,j,k+4,qpres)-q(i,j,k-4,qpres))))*dxinv(3)

             flux(i,j,k,iene)=flux(i,j,k,iene) - &
                   (ALP*(cons(i,j,k+1,iene)*unp1-cons(i,j,k-1,iene)*unm1 &
                  + (q(i,j,k+1,qpres)*unp1-q(i,j,k-1,qpres)*unm1))       &
                  + BET*(cons(i,j,k+2,iene)*unp2-cons(i,j,k-2,iene)*unm2 &
                  + (q(i,j,k+2,qpres)*unp2-q(i,j,k-2,qpres)*unm2))       &
                  + GAM*(cons(i,j,k+3,iene)*unp3-cons(i,j,k-3,iene)*unm3 &
                  + (q(i,j,k+3,qpres)*unp3-q(i,j,k-3,qpres)*unm3))       &
                  + DEL*(cons(i,j,k+4,iene)*unp4-cons(i,j,k-4,iene)*unm4 &
                  + (q(i,j,k+4,qpres)*unp4-q(i,j,k-4,qpres)*unm4)))*dxinv(3)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_hypterm_loop3)

    call destroy(bpt_hypterm)

  end subroutine hypterm



  subroutine diffterm (lo,hi,ng,dx,q,difflux,eta,alam)

    use bl_prof_module

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,6)
    double precision, intent(out) :: difflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)
    double precision, intent(in ) :: eta, alam

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz

    double precision :: dxinv(3)
    double precision :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
    double precision :: divu, uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz,txx,tyy,tzz
    double precision :: mechwork, uxy,uxz,vyz,wzx,wzy,vyx
    integer          :: i,j,k

    double precision, parameter :: OneThird   = 1.0d0/3.0d0
    double precision, parameter :: TwoThirds  = 2.0d0/3.0d0
    double precision, parameter :: FourThirds = 4.0d0/3.0d0

    double precision, parameter :: CENTER = -205.d0/72.d0
    double precision, parameter :: OFF1   =    8.d0/5.d0
    double precision, parameter :: OFF2   =   -0.2d0
    double precision, parameter :: OFF3   =    8.d0/315.d0
    double precision, parameter :: OFF4   =   -1.d0/560.d0

    type(bl_prof_timer), save :: bpt_diffterm, bpt_diffterm_loop123, bpt_diffterm_loop4
    type(bl_prof_timer), save :: bpt_diffterm_loop5, bpt_diffterm_loop6, bpt_diffterm_loop7


    call build(bpt_diffterm, "bpt_diffterm")

    allocate(ux(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(uy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(uz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    difflux(:,:,:,irho) = 0.0d0

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    call build(bpt_diffterm_loop123, "bpt_diffterm_loop123")
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1),hi(1)

             ux(i,j,k)= &
                   (ALP*(q(i+1,j,k,qu)-q(i-1,j,k,qu)) &
                  + BET*(q(i+2,j,k,qu)-q(i-2,j,k,qu)) &
                  + GAM*(q(i+3,j,k,qu)-q(i-3,j,k,qu)) &
                  + DEL*(q(i+4,j,k,qu)-q(i-4,j,k,qu)))*dxinv(1)

             vx(i,j,k)= &
                   (ALP*(q(i+1,j,k,qv)-q(i-1,j,k,qv)) &
                  + BET*(q(i+2,j,k,qv)-q(i-2,j,k,qv)) &
                  + GAM*(q(i+3,j,k,qv)-q(i-3,j,k,qv)) &
                  + DEL*(q(i+4,j,k,qv)-q(i-4,j,k,qv)))*dxinv(1)

             wx(i,j,k)= &
                   (ALP*(q(i+1,j,k,qw)-q(i-1,j,k,qw)) &
                  + BET*(q(i+2,j,k,qw)-q(i-2,j,k,qw)) &
                  + GAM*(q(i+3,j,k,qw)-q(i-3,j,k,qw)) &
                  + DEL*(q(i+4,j,k,qw)-q(i-4,j,k,qw)))*dxinv(1)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng

             uy(i,j,k)= &
                   (ALP*(q(i,j+1,k,qu)-q(i,j-1,k,qu)) &
                  + BET*(q(i,j+2,k,qu)-q(i,j-2,k,qu)) &
                  + GAM*(q(i,j+3,k,qu)-q(i,j-3,k,qu)) &
                  + DEL*(q(i,j+4,k,qu)-q(i,j-4,k,qu)))*dxinv(2)

             vy(i,j,k)= &
                   (ALP*(q(i,j+1,k,qv)-q(i,j-1,k,qv)) &
                  + BET*(q(i,j+2,k,qv)-q(i,j-2,k,qv)) &
                  + GAM*(q(i,j+3,k,qv)-q(i,j-3,k,qv)) &
                  + DEL*(q(i,j+4,k,qv)-q(i,j-4,k,qv)))*dxinv(2)

             wy(i,j,k)= &
                   (ALP*(q(i,j+1,k,qw)-q(i,j-1,k,qw)) &
                  + BET*(q(i,j+2,k,qw)-q(i,j-2,k,qw)) &
                  + GAM*(q(i,j+3,k,qw)-q(i,j-3,k,qw)) &
                  + DEL*(q(i,j+4,k,qw)-q(i,j-4,k,qw)))*dxinv(2)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             uz(i,j,k)= &
                   (ALP*(q(i,j,k+1,qu)-q(i,j,k-1,qu)) &
                  + BET*(q(i,j,k+2,qu)-q(i,j,k-2,qu)) &
                  + GAM*(q(i,j,k+3,qu)-q(i,j,k-3,qu)) &
                  + DEL*(q(i,j,k+4,qu)-q(i,j,k-4,qu)))*dxinv(3)

             vz(i,j,k)= &
                   (ALP*(q(i,j,k+1,qv)-q(i,j,k-1,qv)) &
                  + BET*(q(i,j,k+2,qv)-q(i,j,k-2,qv)) &
                  + GAM*(q(i,j,k+3,qv)-q(i,j,k-3,qv)) &
                  + DEL*(q(i,j,k+4,qv)-q(i,j,k-4,qv)))*dxinv(3)

             wz(i,j,k)= &
                   (ALP*(q(i,j,k+1,qw)-q(i,j,k-1,qw)) &
                  + BET*(q(i,j,k+2,qw)-q(i,j,k-2,qw)) &
                  + GAM*(q(i,j,k+3,qw)-q(i,j,k-3,qw)) &
                  + DEL*(q(i,j,k+4,qw)-q(i,j,k-4,qw)))*dxinv(3)
          enddo
       enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call destroy(bpt_diffterm_loop123)

    call build(bpt_diffterm_loop4, "bpt_diffterm_loop4")
    !$OMP PARALLEL DO PRIVATE(i,j,k,uxx,uyy,uzz,vyx,wzx)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             uxx = (CENTER*q(i,j,k,qu) &
                  + OFF1*(q(i+1,j,k,qu)+q(i-1,j,k,qu)) &
                  + OFF2*(q(i+2,j,k,qu)+q(i-2,j,k,qu)) &
                  + OFF3*(q(i+3,j,k,qu)+q(i-3,j,k,qu)) &
                  + OFF4*(q(i+4,j,k,qu)+q(i-4,j,k,qu)))*dxinv(1)**2

             uyy = (CENTER*q(i,j,k,qu) &
                  + OFF1*(q(i,j+1,k,qu)+q(i,j-1,k,qu)) &
                  + OFF2*(q(i,j+2,k,qu)+q(i,j-2,k,qu)) &
                  + OFF3*(q(i,j+3,k,qu)+q(i,j-3,k,qu)) &
                  + OFF4*(q(i,j+4,k,qu)+q(i,j-4,k,qu)))*dxinv(2)**2

             uzz = (CENTER*q(i,j,k,qu) &
                  + OFF1*(q(i,j,k+1,qu)+q(i,j,k-1,qu)) &
                  + OFF2*(q(i,j,k+2,qu)+q(i,j,k-2,qu)) &
                  + OFF3*(q(i,j,k+3,qu)+q(i,j,k-3,qu)) &
                  + OFF4*(q(i,j,k+4,qu)+q(i,j,k-4,qu)))*dxinv(3)**2

             vyx = (ALP*(vy(i+1,j,k)-vy(i-1,j,k)) &
                  + BET*(vy(i+2,j,k)-vy(i-2,j,k)) &
                  + GAM*(vy(i+3,j,k)-vy(i-3,j,k)) &
                  + DEL*(vy(i+4,j,k)-vy(i-4,j,k)))*dxinv(1)

             wzx = (ALP*(wz(i+1,j,k)-wz(i-1,j,k)) &
                  + BET*(wz(i+2,j,k)-wz(i-2,j,k)) &
                  + GAM*(wz(i+3,j,k)-wz(i-3,j,k)) &
                  + DEL*(wz(i+4,j,k)-wz(i-4,j,k)))*dxinv(1)

             difflux(i,j,k,imx) = eta*(FourThirds*uxx + uyy + uzz + OneThird*(vyx+wzx))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_diffterm_loop4)

    call build(bpt_diffterm_loop5, "bpt_diffterm_loop5")
    !$OMP PARALLEL DO PRIVATE(i,j,k,vxx,vyy,vzz,uxy,wzy)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             vxx = (CENTER*q(i,j,k,qv) &
                  + OFF1*(q(i+1,j,k,qv)+q(i-1,j,k,qv)) &
                  + OFF2*(q(i+2,j,k,qv)+q(i-2,j,k,qv)) &
                  + OFF3*(q(i+3,j,k,qv)+q(i-3,j,k,qv)) &
                  + OFF4*(q(i+4,j,k,qv)+q(i-4,j,k,qv)))*dxinv(1)**2

             vyy = (CENTER*q(i,j,k,qv) &
                  + OFF1*(q(i,j+1,k,qv)+q(i,j-1,k,qv)) &
                  + OFF2*(q(i,j+2,k,qv)+q(i,j-2,k,qv)) &
                  + OFF3*(q(i,j+3,k,qv)+q(i,j-3,k,qv)) &
                  + OFF4*(q(i,j+4,k,qv)+q(i,j-4,k,qv)))*dxinv(2)**2

             vzz = (CENTER*q(i,j,k,qv) &
                  + OFF1*(q(i,j,k+1,qv)+q(i,j,k-1,qv)) &
                  + OFF2*(q(i,j,k+2,qv)+q(i,j,k-2,qv)) &
                  + OFF3*(q(i,j,k+3,qv)+q(i,j,k-3,qv)) &
                  + OFF4*(q(i,j,k+4,qv)+q(i,j,k-4,qv)))*dxinv(3)**2

             uxy = (ALP*(ux(i,j+1,k)-ux(i,j-1,k)) &
                  + BET*(ux(i,j+2,k)-ux(i,j-2,k)) &
                  + GAM*(ux(i,j+3,k)-ux(i,j-3,k)) &
                  + DEL*(ux(i,j+4,k)-ux(i,j-4,k)))*dxinv(2)

             wzy = (ALP*(wz(i,j+1,k)-wz(i,j-1,k)) &
                  + BET*(wz(i,j+2,k)-wz(i,j-2,k)) &
                  + GAM*(wz(i,j+3,k)-wz(i,j-3,k)) &
                  + DEL*(wz(i,j+4,k)-wz(i,j-4,k)))*dxinv(2)

             difflux(i,j,k,imy) = eta*(vxx + FourThirds*vyy + vzz + OneThird*(uxy+wzy))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_diffterm_loop5)

    call build(bpt_diffterm_loop6, "bpt_diffterm_loop6")
    !$OMP PARALLEL DO PRIVATE(i,j,k,wxx,wyy,wzz,uxz,vyz)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             wxx = (CENTER*q(i,j,k,qw) &
                  + OFF1*(q(i+1,j,k,qw)+q(i-1,j,k,qw)) &
                  + OFF2*(q(i+2,j,k,qw)+q(i-2,j,k,qw)) &
                  + OFF3*(q(i+3,j,k,qw)+q(i-3,j,k,qw)) &
                  + OFF4*(q(i+4,j,k,qw)+q(i-4,j,k,qw)))*dxinv(1)**2

             wyy = (CENTER*q(i,j,k,qw) &
                  + OFF1*(q(i,j+1,k,qw)+q(i,j-1,k,qw)) &
                  + OFF2*(q(i,j+2,k,qw)+q(i,j-2,k,qw)) &
                  + OFF3*(q(i,j+3,k,qw)+q(i,j-3,k,qw)) &
                  + OFF4*(q(i,j+4,k,qw)+q(i,j-4,k,qw)))*dxinv(2)**2

             wzz = (CENTER*q(i,j,k,qw) &
                  + OFF1*(q(i,j,k+1,qw)+q(i,j,k-1,qw)) &
                  + OFF2*(q(i,j,k+2,qw)+q(i,j,k-2,qw)) &
                  + OFF3*(q(i,j,k+3,qw)+q(i,j,k-3,qw)) &
                  + OFF4*(q(i,j,k+4,qw)+q(i,j,k-4,qw)))*dxinv(3)**2

             uxz = (ALP*(ux(i,j,k+1)-ux(i,j,k-1)) &
                  + BET*(ux(i,j,k+2)-ux(i,j,k-2)) &
                  + GAM*(ux(i,j,k+3)-ux(i,j,k-3)) &
                  + DEL*(ux(i,j,k+4)-ux(i,j,k-4)))*dxinv(3)

             vyz = (ALP*(vy(i,j,k+1)-vy(i,j,k-1)) &
                  + BET*(vy(i,j,k+2)-vy(i,j,k-2)) &
                  + GAM*(vy(i,j,k+3)-vy(i,j,k-3)) &
                  + DEL*(vy(i,j,k+4)-vy(i,j,k-4)))*dxinv(3)

             difflux(i,j,k,imz) = eta*(wxx + wyy + FourThirds*wzz + OneThird*(uxz+vyz))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_diffterm_loop6)

    call build(bpt_diffterm_loop7, "bpt_diffterm_loop7")
    !$OMP PARALLEL DO PRIVATE(i,j,k,txx,tyy,tzz) &
    !$OMP PRIVATE(divu,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz,mechwork)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             txx = (CENTER*q(i,j,k,6) &
                  + OFF1*(q(i+1,j,k,6)+q(i-1,j,k,6)) &
                  + OFF2*(q(i+2,j,k,6)+q(i-2,j,k,6)) &
                  + OFF3*(q(i+3,j,k,6)+q(i-3,j,k,6)) &
                  + OFF4*(q(i+4,j,k,6)+q(i-4,j,k,6)))*dxinv(1)**2

             tyy = (CENTER*q(i,j,k,6) &
                  + OFF1*(q(i,j+1,k,6)+q(i,j-1,k,6)) &
                  + OFF2*(q(i,j+2,k,6)+q(i,j-2,k,6)) &
                  + OFF3*(q(i,j+3,k,6)+q(i,j-3,k,6)) &
                  + OFF4*(q(i,j+4,k,6)+q(i,j-4,k,6)))*dxinv(2)**2

             tzz = (CENTER*q(i,j,k,6) &
                  + OFF1*(q(i,j,k+1,6)+q(i,j,k-1,6)) &
                  + OFF2*(q(i,j,k+2,6)+q(i,j,k-2,6)) &
                  + OFF3*(q(i,j,k+3,6)+q(i,j,k-3,6)) &
                  + OFF4*(q(i,j,k+4,6)+q(i,j,k-4,6)))*dxinv(3)**2

             divu  = TwoThirds*(ux(i,j,k)+vy(i,j,k)+wz(i,j,k))
             tauxx = 2.d0*ux(i,j,k) - divu
             tauyy = 2.d0*vy(i,j,k) - divu
             tauzz = 2.d0*wz(i,j,k) - divu
             tauxy = uy(i,j,k)+vx(i,j,k)
             tauxz = uz(i,j,k)+wx(i,j,k)
             tauyz = vz(i,j,k)+wy(i,j,k)

             mechwork = tauxx*ux(i,j,k) + &
                        tauyy*vy(i,j,k) + &
                        tauzz*wz(i,j,k) + tauxy**2+tauxz**2+tauyz**2

             mechwork = eta*mechwork &
                  + difflux(i,j,k,imx)*q(i,j,k,qu) &
                  + difflux(i,j,k,imy)*q(i,j,k,qv) &
                  + difflux(i,j,k,imz)*q(i,j,k,qw)

             difflux(i,j,k,iene) = alam*(txx+tyy+tzz) + mechwork
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_diffterm_loop7)

    call destroy(bpt_diffterm)

  end subroutine diffterm

end module advance_module

