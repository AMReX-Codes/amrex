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

    integer        :: n, nc, ng, lo(2), hi(2)
    type(multifab) :: f
    type(layout)   :: la

    double precision, pointer, dimension(:,:,:,:) :: qp, fp 

    nc = ncomp(q)
    ng = nghost(q)
    la = get_layout(q)

    ! sync q and fill ghost cells
    call fill_boundary(q)

    ! build flux multifab
    call build(f, la, nc, 0)

    do n=1, nboxes(q)
       if ( remote(q,n) ) cycle

       qp => dataptr(q,n)
       fp => dataptr(f,n)

       lo = lwb(get_box(q,n))
       hi = upb(get_box(q,n))

       fp = 0.0d0

       ! XXX: pass invdx here instead of dx

       call cell_motility           (fp(:,:,1,iu), qp(:,:,1,iu),               lo, hi, ng, ctx%dx, ctx%diff)
       call chemotactic_sensitivity (fp(:,:,1,iu), qp(:,:,1,iu), qp(:,:,1,iv), lo, hi, ng, ctx%dx, ctx%chi)
       call signal_diffusion        (fp(:,:,1,iv), qp(:,:,1,iv),               lo, hi, ng, ctx%dx)
       call signal_production       (fp(:,:,1,iv), qp(:,:,1,iu), qp(:,:,1,iv), lo, hi, ng, ctx%dx)
       call signal_degradation      (fp(:,:,1,iv), qp(:,:,1,iu), qp(:,:,1,iv), lo, hi, ng, ctx%dx)
    end do

    call saxpy(q, dt, f)

    ! cleanup
    call destroy(f)

  end subroutine advance_fe


  subroutine advance_sdc(q, dt, ctx, sdc)

    type(multifab),  intent(inout) :: q
    real(dp_t),      intent(in)    :: dt
    type(cht_ctx_t), intent(inout) :: ctx
    type(sdcquad),   intent(in)    :: sdc

    
  end subroutine advance_sdc

end module advance_module

