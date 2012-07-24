module advance_module
  
  use multifab_module
  use dtypes_module
  use sdcquad_module
  use write_plotfile_module

  implicit none

  ! Throughout this module we use the fourth-order stencils from
  ! Kadioglu, Klein, and Minion 2008.
  !
  ! The edge-average of a derivative normal to the edge is given by,
  ! eg,
  !
  !   (\phi_x)_{i+1/2,j} = [ 1, -15, 15, -1 ] / 12 h
  !
  ! for offsets [ -1, 0, 1, 2 ] so that
  !  
  !   (\phi_x)_{i+1/2,j} - (\phi_x)_{i-1/2,j} = [ 1, -14, 0, 14, -1 ] / 12 h
  !
  ! for centered offsets.

  double precision, parameter :: dxm2 =   1.0d0/12
  double precision, parameter :: dxm1 = -14.0d0/12
  double precision, parameter :: dxp1 =  14.0d0/12
  double precision, parameter :: dxp2 =  -1.0d0/12


contains
  
  subroutine advance_fe(q, dt, ctx)

    type(multifab),  intent(inout) :: q
    real(dp_t),      intent(in)    :: dt
    type(cht_ctx_t), intent(inout) :: ctx

    integer        :: n, lo(2), hi(2)
    type(multifab) :: f

    double precision, pointer, dimension(:,:,:,:) :: qp, fp 

    ! sync q and fill ghost cells
    call fill_boundary(q)

    ! build flux multifab
    call build(f, ctx%la, ctx%nc, 0)

    do n=1, nboxes(q)
       if ( remote(q,n) ) cycle

       qp => dataptr(q,n)
       fp => dataptr(f,n)

       lo = lwb(get_box(q,n))
       hi = upb(get_box(q,n))

       call cell_diffusive_flux_2d(qp, lo, hi, ctx%ng, ctx, fp)

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


  subroutine cell_diffusive_flux_2d(qp, lo, hi, ng, ctx, fp)
    integer,          intent(in)    :: lo(2), hi(2), ng
    type(cht_ctx_t),  intent(inout) :: ctx
    double precision, intent(in)    :: qp(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,2)
    double precision, intent(out)   :: fp(lo(1):hi(1),lo(2):hi(2),2)

    integer :: i, j 

    ! diffusive flux: edge difference of: grad u - chi u grad v

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! \nabla u
          fp(i,j,1) = ( &
               + dxm2 * (qp(i-2,j,1) + qp(i,j-2,1)) &
               + dxm1 * (qp(i-1,j,1) + qp(i,j-1,1)) &
               + dxp1 * (qp(i+1,j,1) + qp(i,j+1,1)) &
               + dxp2 * (qp(i+2,j,1) + qp(i,j+2,1)) &
               ) / ctx%dx
          
       end do
    end do
  end subroutine cell_diffusive_flux_2d


end module advance_module

