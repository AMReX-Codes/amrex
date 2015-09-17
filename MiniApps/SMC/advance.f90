module advance_module

  use bl_error_module
  use kernels_module
  use multifab_module
  use time_module
  use transport_properties
  use variables_module

  implicit none

  double precision, save, public :: wt_fillboundary = 0.d0
  double precision, save, public :: wt_ctoprim = 0.d0
  double precision, save, public :: wt_chemterm = 0.d0
  double precision, save, public :: wt_transprop = 0.d0
  double precision, save, public :: wt_hypdiff = 0.d0

  private
  public advance

contains

  subroutine advance(U, dt, courno, dx, istep)
    use smcdata_module
    type(multifab),    intent(inout) :: U
    double precision,  intent(inout) :: dt, courno
    double precision,  intent(in   ) :: dx(3)
    integer, intent(in) :: istep

    call multifab_setval(Unew, 0.d0, .true.)

    call dUdt(U, Uprime, Q, mu, xi, lam, Ddiag, dx, courno, istep)
    call set_dt(dt, courno, istep)
    call update_rk3(Zero,Unew, One,U, dt,Uprime)
    call reset_density(Unew)

    call dUdt(Unew, Uprime, Q, mu, xi, lam, Ddiag, dx)
    call update_rk3(OneQuarter, Unew, ThreeQuarters, U, OneQuarter*dt, Uprime)
    call reset_density(Unew)

    call dUdt(Unew, Uprime, Q, mu, xi, lam, Ddiag, dx)
    call update_rk3(OneThird, U, TwoThirds, Unew, TwoThirds*dt, Uprime)
    call reset_density(U)
  end subroutine advance


  !
  ! Compute new time-step size
  !
  subroutine set_dt(dt, courno, istep)

    use probin_module, only : cflfac, fixed_dt, init_shrink, max_dt, max_dt_growth, small_dt, stop_time

    double precision, intent(inout) :: dt
    double precision, intent(in   ) :: courno
    integer,          intent(in   ) :: istep

    double precision :: dtold

    if (fixed_dt > 0.d0) then

       dt = fixed_dt

       if (parallel_IOProcessor()) then
          print*, ""
          print*, "   Setting fixed dt =",dt
          print*, ""
       end if

    else

       dtold = dt
       dt    = cflfac / courno

       if (parallel_IOProcessor()) then
          print*, "   CFL: dt =", dt
       end if

       if (istep .eq. 1) then
          dt = dt * init_shrink
          if (parallel_IOProcessor()) then
             print*,'   Limited by init_shrink: dt =',dt
          end if
       else
          if (dt .gt. dtold * max_dt_growth) then
             dt = dtold * max_dt_growth
             if (parallel_IOProcessor()) then
                print*,'   Limited by dt_growth: dt =',dt
             end if
          end if
       end if

       if(dt .gt. max_dt) then
          if (parallel_IOProcessor()) then
             print*,'   Limited by max_dt: dt =',max_dt
          end if
          dt = max_dt
       end if

       if (dt < small_dt) then
          call bl_error("ERROR: timestep < small_dt")
       endif

       if (stop_time > 0.d0) then
          if (time + dt > stop_time) then
             dt = stop_time - time
             if (parallel_IOProcessor()) then
                print*, "   Limited by stop_time: dt =",dt
             end if
          end if
       end if
       
    end if

  end subroutine set_dt


  !
  ! Compute U1 = a U1 + b U2 + c Uprime.
  !
  subroutine update_rk3 (a,U1,b,U2,c,Uprime)

    type(multifab),   intent(in   ) :: U2, Uprime
    type(multifab),   intent(inout) :: U1
    double precision, intent(in   ) :: a, b, c

    integer :: lo(3), hi(3), i, j, k, m, n, nc
    type(mfiter) :: mfi
    type(box) :: tbx
    double precision, pointer, dimension(:,:,:,:) :: u1p, u2p, upp

    nc = ncomp(U1)

    !$omp parallel private(i,j,k,m,n,lo,hi,mfi,tbx,u1p,u2p,upp)

    call mfiter_build(mfi, U1, tiling=.true.)

    do while (next_tile(mfi,n))

       tbx = get_tilebox(mfi)
       lo = lwb(tbx)
       hi = upb(tbx)

       u1p => dataptr(U1,    n)
       u2p => dataptr(U2,    n)
       upp => dataptr(Uprime,n)

       do m = 1, nc
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   u1p(i,j,k,m) = a * u1p(i,j,k,m) + b * u2p(i,j,k,m) + c * upp(i,j,k,m)
                end do
             end do
          end do
       end do
    end do
    !$omp end parallel

  end subroutine update_rk3


  !
  ! Compute dU/dt given U.
  !
  ! The Courant number (courno) is also computed if passed.
  !
  subroutine dUdt (U, Uprime, Q, mu, xi, lam, Ddiag, dx, courno, istep)

    use probin_module, only: cfl_int, fixed_dt, do_tiling

    type(multifab),   intent(inout) :: U, Uprime, Q, mu, xi, lam, Ddiag
    double precision, intent(in   ) :: dx(3)
    integer,          intent(in   ), optional :: istep
    double precision, intent(inout), optional :: courno

    integer :: lo(3), hi(3)
    integer :: n, ng
    type(mfiter) :: mfi
    type(box) :: tbx

    logical :: update_courno
    double precision :: courno_proc

    integer :: qlo(4), qhi(4), uplo(4), uphi(4), ulo(4), uhi(4)
    double precision, pointer, dimension(:,:,:,:) :: up, qp, mup, xip, lamp, Ddp, upp

    double precision :: wt1, wt2

    update_courno = .false.
    if (present(courno) .and. present(istep) .and. fixed_dt.le.0.d0) then
       if (mod(istep,cfl_int).eq.1 .or. cfl_int.le.1) then
          update_courno = .true.
       end if
    end if

    wt1 = parallel_wtime()
    call multifab_fill_boundary(U)
    wt_fillboundary = wt_fillboundary + (parallel_wtime()-wt1)

    call multifab_setval(Uprime, 0.d0)

    !
    ! Calculate primitive variables based on U
    !
    wt1 = parallel_wtime()
    call ctoprim(U, Q)
    wt_ctoprim = wt_ctoprim + (parallel_wtime()-wt1)

    if (update_courno) then
       courno_proc = -1.d50
       call compute_courno(Q, dx, courno_proc)
    end if

    ! 
    ! chemistry
    !
    wt1 = parallel_wtime()
    do n=1,nfabs(Q)

       qp  => dataptr(Q,n)
       upp => dataptr(Uprime,n)

       qlo = lbound(qp)
       qhi = ubound(qp)
       uplo = lbound(upp)
       uphi = ubound(upp)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call chemterm_3d(lo,hi,qp,qlo(1:3),qhi(1:3),upp,uplo(1:3),uphi(1:3))
    end do
    wt2 = parallel_wtime()
    wt_chemterm = wt_chemterm + (wt2-wt1)

    !
    ! transport coefficients
    !
    call get_transport_properties(Q, mu, xi, lam, Ddiag)
    wt1 = parallel_wtime()
    wt_transprop = wt_transprop + (wt1-wt2)

    !
    ! Hyperbolic and Transport terms
    !
    if (do_tiling) then
       !$omp parallel private(n,lo,hi,up,ulo,uhi,upp,uplo,uphi) &
       !$omp private(qp,qlo,qhi,mup,xip,lamp,Ddp,mfi,tbx)
       
       call mfiter_build(mfi, Q, tiling=.true.)
       
       do while (next_tile(mfi,n))
          
          tbx = get_tilebox(mfi)
          lo = lwb(tbx)
          hi = upb(tbx)
          
          up => dataptr(U,n)
          upp=> dataptr(Uprime,n)
          qp => dataptr(Q,n)
          mup  => dataptr(mu   , n)
          xip  => dataptr(xi   , n)
          lamp => dataptr(lam  , n)
          Ddp  => dataptr(Ddiag, n)
          
          ulo = lbound(up)
          uhi = ubound(up)
          qlo = lbound(qp)
          qhi = ubound(qp)
          uplo = lbound(upp)
          uphi = ubound(upp)
          
          call narrow_diffterm_3d(lo,hi,dx,qp,qlo(1:3),qhi(1:3),upp,uplo(1:3),uphi(1:3), &
               mup,xip,lamp,Ddp)
          
          call hypterm_3d(lo,hi,dx,up,ulo(1:3),uhi(1:3),qp,qlo(1:3),qhi(1:3),&
               upp,uplo(1:3),uphi(1:3))
          
       end do
       !$omp end parallel
    else

       ng = nghost(U)

       do n=1,nfabs(Q)
          qp  => dataptr(Q,n)
          up  => dataptr(U,n)
          upp => dataptr(Uprime,n)
          
          mup  => dataptr(mu   , n)
          xip  => dataptr(xi   , n)
          lamp => dataptr(lam  , n)
          Ddp  => dataptr(Ddiag, n)
          
          lo = lwb(get_box(Q,n))
          hi = upb(get_box(Q,n))
          
          call orig_narrow_diffterm_3d(lo,hi,ng,dx,qp,upp,mup,xip,lamp,Ddp)
          
          call orig_hypterm_3d(lo,hi,ng,dx,up,qp,upp)
       end do
       
    end if

    wt2 = parallel_wtime()
    wt_hypdiff = wt_hypdiff + (wt2-wt1)

    if (update_courno) then
       call parallel_reduce(courno, courno_proc, MPI_MAX)
    end if

  end subroutine dUdt

  subroutine compute_courno(Q, dx, courno)
    type(multifab), intent(in) :: Q
    double precision, intent(in) :: dx(3)
    double precision, intent(inout) :: courno

    integer :: n, lo(3), hi(3), qlo(4), qhi(4)
    type(mfiter) :: mfi
    type(box) :: tbx
    double precision, pointer :: qp(:,:,:,:)

    !$omp parallel private(mfi, tbx, n, lo, hi, qlo, qhi, qp) &
    !$omp reduction(max:courno)

    call mfiter_build(mfi, Q, tiling=.true.)

    do while (next_tile(mfi,n))

       tbx = get_tilebox(mfi)
       lo = lwb(tbx)
       hi = upb(tbx)

       qp => dataptr(Q,n)
       qlo = lbound(qp)
       qhi = ubound(qp)

       call comp_courno_3d(lo,hi,dx,qp,qlo(1:3),qhi(1:3),courno)

    end do
    !$omp end parallel
  end subroutine compute_courno

end module advance_module
