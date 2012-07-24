module chemotaxis_module

  use boxlib
  use multifab_module

  use dtypes_module
  use advance_module
  use sdcquad_module
  use write_plotfile_module

  implicit none

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cht_main(ctx, opts, sdc)
    type(cht_ctx_t),  intent(inout) :: ctx
    type(cht_opts_t), intent(inout) :: opts
    type(sdcquad),    intent(in)    :: sdc

    integer,    allocatable :: lo(:), hi(:)
    real(dp_t), allocatable :: prob_lo(:), prob_hi(:)
    logical,    allocatable :: is_periodic(:)
    
    integer        :: n
    real(dp_t)     :: dt, time
  
    type(box)      :: bx
    type(boxarray) :: ba
    type(layout)   :: la
    type(multifab) :: q

    ! allocate
    allocate(lo(ctx%dim),hi(ctx%dim))
    allocate(prob_lo(ctx%dim),prob_hi(ctx%dim))
    allocate(is_periodic(ctx%dim))

    ! physical problem is a box on (-1,-1) to (1,1), periodic on all sides
    prob_lo     = -1.d0
    prob_hi     =  1.d0
    is_periodic = .true.

    ! create a box from (0,0) to (n_cell-1,n_cell-1)
    lo = 0
    hi = ctx%n_cell-1
    bx = make_box(lo,hi)

    ctx%dx = (prob_hi(1)-prob_lo(1)) / ctx%n_cell

    ! build layout
    call build(ba,bx)
    call build(la,ba,bx,pmask=is_periodic)

    ctx%la = la

    ! build q (solution) multifab
    call build(q,la,ctx%nc,ctx%ng)
    
    ! set initial condition and write plot 0
    call cht_initial(q,ctx)
    call write_plotfile(la,q,0,ctx%dx,0.0d0,prob_lo,prob_hi)

    ! run
    select case(opts%method)
    ! case("rk")
    !    call advance_rk(q,dt,ctx)
    case("fe")
       time = 0.0d0
       do n = 1, 100
          call advance_fe(q,dt,ctx)
          time = time + dt
       end do
       call write_plotfile(la,q,n,ctx%dx,time,prob_lo,prob_hi)

    case("sdc")
       call advance_sdc(q,dt,ctx,sdc)
    end select

    ! destroy/deallocate
    call destroy(ba)
    call destroy(la)
    call destroy(q)

    deallocate(lo,hi,prob_lo,prob_hi,is_periodic)

  end subroutine cht_main

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cht_initial(q,ctx)
    type(multifab),  intent(inout) :: q
    type(cht_ctx_t), intent(in)    :: ctx

    integer        :: n, i, j, lo(2), hi(2)

    double precision :: x, y
    double precision, pointer, dimension(:,:,:,:) :: qp

    do n=1, nboxes(q)
       if ( remote(q,n) ) cycle

       qp => dataptr(q,n)

       lo = lwb(get_box(q,n))
       hi = upb(get_box(q,n))

       do j = lo(2), hi(2)
          y = -1.0d0 + j * ctx%dx
          do i = lo(1), hi(2)
             x = -1.0d0 + i * ctx%dx

             qp(i,j,1,1) = 1.0d0 + sin(x)
          end do
       end do

    end do

  end subroutine cht_initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module chemotaxis_module
