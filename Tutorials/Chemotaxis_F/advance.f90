module advance_module
  
  use multifab_module
  use chemotaxis_kernels
  use dtypes_module
  use sdcquad_module
  use write_plotfile_module

  implicit none

contains
  
  subroutine advance_fe(q, dt, ctx)

    type(multifab),  intent(inout) :: q
    real(dp_t),      intent(in)    :: dt
    type(cht_ctx_t), intent(inout) :: ctx

    integer        :: nc
    type(multifab) :: f
    type(layout)   :: la

    nc = ncomp(q)
    la = get_layout(q)

    ! sync q and fill ghost cells
    call fill_boundary(q)

    ! build flux multifab
    call build(f, la, nc, 0)

    ! compute forward euler update
    call dqdt(q, ctx, f)
    call saxpy(q, dt, f)

    ! cleanup
    call destroy(f)

  end subroutine advance_fe


  subroutine advance_sdc(q, dt, ctx, sdc)

    type(multifab),  intent(inout) :: q
    real(dp_t),      intent(inout) :: dt
    type(cht_ctx_t), intent(inout) :: ctx
    type(sdcquad),   intent(in)    :: sdc

    integer          :: k, m, nc, ng
    double precision :: res, res_local
    type(layout)     :: la
    type(multifab)   :: qSDC(sdc%nnodes), fSDC(sdc%nnodes)

    nc = ncomp(q)
    ng = nghost(q)
    la = get_layout(q)

    ! set provisional solution
    do m = 1, sdc%nnodes
       call build(qSDC(m), la, nc, ng)
       call build(fSDC(m), la, nc, 0)
    end do

    call copy(qSDC(1), q)
    call dqdt(qSDC(1), ctx, fSDC(1))

    do m = 2, sdc%nnodes
       call copy(qSDC(m), qSDC(1))
       call copy(fSDC(m), fSDC(1))
    end do

    res = 0.0d0

    ! perform sdc iterations
    do k = 1, sdc%iters
       call sdc_sweep(qSDC, fSDC, dt, ctx, sdc)

       if (sdc%tol_residual > 0.d0) then
          res_local = sdc_residual(qSDC, fSDC, dt, sdc)
          call parallel_reduce(res, res_local, MPI_MAX)

          if (res < sdc%tol_residual) then
             exit
          end if
       end if
    end do

    print *, "SDC: iters:", k, "residual:", res
 
    call copy(q, qSDC(sdc%nnodes))
    
    do m = 1, sdc%nnodes
       call destroy(qSDC(m))
       call destroy(fSDC(m))
    end do
    
  end subroutine advance_sdc

  !
  ! Perform one SDC sweep.
  !
  subroutine sdc_sweep (qSDC,fSDC,dt,ctx,sdc)
    type(sdcquad),     intent(in   ) :: sdc
    type(multifab),    intent(inout) :: qSDC(sdc%nnodes), fSDC(sdc%nnodes)
    double precision,  intent(in   ) :: dt
    type(cht_ctx_t),   intent(in   ) :: ctx

    integer        :: m, n, nc
    type(multifab) :: S(sdc%nnodes-1)
    type(layout)   :: la

    double precision :: dtsdc(sdc%nnodes-1)

    la = get_layout(qSDC(1))
    nc = ncomp(qSDC(1))

    !
    ! Compute integrals (compact forward Euler)
    ! 
    do m = 1, sdc%nnodes-1
       call build(S(m), la, nc, 0)
       call setval(S(m), 0.0d0)
       do n = 1, sdc%nnodes
          call saxpy(S(m), sdc%smats(m,n,1), fSDC(n))
       end do
    end do

    !
    ! Perform sub-step correction
    !
    dtsdc = dt * (sdc%nodes(2:sdc%nnodes) - sdc%nodes(1:sdc%nnodes-1))
    do m = 1, sdc%nnodes-1

       ! U(m+1) = U(m) + dt dUdt(m) + dt S(m)

       call copy(qSDC(m+1), qSDC(m))
       call saxpy(qSDC(m+1), dtsdc(m), fSDC(m))
       call saxpy(qSDC(m+1), dt, S(m))

       call dqdt(qSDC(m+1), ctx, fSDC(m+1))

    end do


    do m = 1, sdc%nnodes-1
       call destroy(S(m))
    end do

  end subroutine sdc_sweep


  !
  ! Compute SDC residual.
  !
  function sdc_residual (qSDC,fSDC,dt,sdc) result(res)
    real(dp_t)                      :: res
    type(sdcquad),    intent(in   ) :: sdc
    type(multifab),   intent(inout) :: qSDC(sdc%nnodes), fSDC(sdc%nnodes)
    real(dp_t),       intent(in   ) :: dt      

    integer        :: m, n, nc
    type(multifab) :: R
    type(layout)   :: la

    la = get_layout(qSDC(1))
    nc = ncomp(qSDC(1))

    !
    ! Compute integral
    ! 
    call build(R, la, nc, 0)
    call copy(R, qSDC(1))

    do m = 1, sdc%nnodes-1
       do n = 1, sdc%nnodes
          call saxpy(R, dt*sdc%smat(m,n), fSDC(n))
       end do
    end do

    call saxpy(R, -1.0d0, qSDC(sdc%nnodes))
    
    res = norm_inf(R)

    call destroy(R)

  end function sdc_residual



  subroutine dqdt(q,ctx,f)
    type(multifab),   intent(inout) :: q, f
    type(cht_ctx_t),  intent(in)    :: ctx

    integer        :: n, ng, lo(2), hi(2)

    double precision, pointer, dimension(:,:,:,:) :: qp, fp 

    ng = nghost(q)

    ! sync q and fill ghost cells
    call fill_boundary(q)

    do n=1, nfabs(q)

       qp => dataptr(q,n)
       fp => dataptr(f,n)

       lo = lwb(get_box(q,n))
       hi = upb(get_box(q,n))

       fp = 0.0d0

       ! XXX: pass invdx here instead of dx

       call cell_motility           (fp(:,:,1,iu), qp(:,:,1,iu),               lo, hi, ng, ctx%dx, ctx%diff)
       call chemotactic_sensitivity (fp(:,:,1,iu), qp(:,:,1,iu), qp(:,:,1,iv), lo, hi, ng, ctx%dx, ctx%chi, ctx%alpha, ctx%gamma)
       call signal_diffusion        (fp(:,:,1,iv),               qp(:,:,1,iv), lo, hi, ng, ctx%dx)
       call signal_production       (fp(:,:,1,iv), qp(:,:,1,iu), qp(:,:,1,iv), lo, hi, ng, ctx%dx, ctx%phi)
       call signal_degradation      (fp(:,:,1,iv), qp(:,:,1,iu), qp(:,:,1,iv), lo, hi, ng, ctx%dx)

    end do

  end subroutine dqdt

end module advance_module

