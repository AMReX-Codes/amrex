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

  subroutine chemotaxis(ctx, opts, sdc)
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

    dt = 0.01d0

    ! create a box from (0,0) to (n_cell-1,n_cell-1)
    lo = 0
    hi = ctx%n_cell-1
    bx = make_box(lo,hi)

    ctx%dx    = (prob_hi(1)-prob_lo(1)) / ctx%n_cell
    ctx%invdx = 1.0d0 / ctx%dx

    ! build layout
    call build(ba,bx)
    call build(la,ba,bx,pmask=is_periodic)

    ! build q (solution) multifab
    call build(q,la,ctx%nc,ctx%ng)
    
    ! set initial condition and write plot 0
    call initial(q,ctx)
    call write_plotfile(la,q,0,ctx%dx,0.0d0,prob_lo,prob_hi)

    ! run
    time = 0.0d0
    do n = 1, opts%nsteps
       select case(opts%method)
       case("fe")
          call advance_fe(q,dt,ctx)
       case("sdc")
          call advance_sdc(q,dt,ctx,sdc)
       end select

       time = time + dt
       if (mod(n, opts%plot_int) .eq. 0 .or. n .eq. opts%nsteps) then
          call write_plotfile(la,q,n,ctx%dx,time,prob_lo,prob_hi)
       end if
    end do

    ! destroy/deallocate
    call destroy(ba)
    call destroy(la)
    call destroy(q)

    deallocate(lo,hi,prob_lo,prob_hi,is_periodic)

  end subroutine chemotaxis

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initial(q,ctx)
    type(multifab),  intent(inout) :: q
    type(cht_ctx_t), intent(in)    :: ctx

    integer        :: n, i, j, lo(2), hi(2)

    double precision :: x, y, d
    double precision, pointer, dimension(:,:,:,:) :: qp

    ! double precision, parameter :: pi = 3.141592653589793d0

    do n=1, nboxes(q)
       if ( remote(q,n) ) cycle

       qp => dataptr(q,n)

       lo = lwb(get_box(q,n))
       hi = upb(get_box(q,n))

       do j = lo(2), hi(2)
          y = -1.0d0 + j * ctx%dx
          do i = lo(1), hi(2)
             x = -1.0d0 + i * ctx%dx

             qp(i,j,1,iu) = 1.0d0

             d = (x)**2 + (y)**2
             qp(i,j,1,iv) = 1.0d0 + 0.1d0 * dexp(-10.0d0 * d)
          end do
       end do

    end do

  end subroutine initial

end module chemotaxis_module
