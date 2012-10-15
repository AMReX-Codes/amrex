module chemotaxis_module

  use boxlib
  use multifab_module
  use mt19937_module

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

    integer        :: lo(dim), hi(dim), n
    logical        :: is_periodic(dim)
    real(dp_t)     :: dt, time 
  
    type(box)      :: bx
    type(boxarray) :: ba
    type(layout)   :: la
    type(multifab) :: q

    ! physical problem is a box on (-1,-1) to (1,1), periodic on all sides
    ctx%prob_lo     = -10.0d0
    ctx%prob_hi     =  10.0d0
    is_periodic = .true.

    ! create a box from (0,0) to (n_cell-1,n_cell-1)
    lo = 0
    hi = ctx%n_cell-1
    bx = make_box(lo,hi)

    ctx%dx    = (ctx%prob_hi(1) - ctx%prob_lo(1)) / ctx%n_cell
    ctx%invdx = 1.0d0 / ctx%dx

    ! build layout
    call build(ba,bx)
    call build(la,ba,bx,pmask=is_periodic)

    ! build q (solution) multifab
    call build(q,la,ctx%nc,ctx%ng)
    
    ! set initial condition and write plot 0
    call initial(q,ctx)
    call write_plotfile(la,q,0,ctx%dx,0.0d0,ctx%prob_lo,ctx%prob_hi)

    ! run
    dt = ctx%dt
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
          print *, "OUTPUT: step:", n, "time:", time
          call write_plotfile(la,q,n,ctx%dx,time,ctx%prob_lo,ctx%prob_hi)
       end if
    end do

    ! destroy/deallocate
    call destroy(ba)
    call destroy(la)
    call destroy(q)

  end subroutine chemotaxis

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initial(q,ctx)
    type(multifab),  intent(inout) :: q
    type(cht_ctx_t), intent(in)    :: ctx

    integer        :: n, p, i, j, lo(2), hi(2)

    double precision :: x, y, px, py, d, a
    double precision, pointer, dimension(:,:,:,:) :: qp

    ! double precision, parameter :: pi  = 3.141592653589793d0
    ! double precision, parameter :: pad = 4.0d0
    double precision, parameter :: pad = 2.0d0

    call init_genrand(368)

    do n=1, nfabs(q)

       qp => dataptr(q,n)

       lo = lwb(get_box(q,n))
       hi = upb(get_box(q,n))

       qp(:,:,1,iu) = 1.0d0
       qp(:,:,1,iv) = 0.5d0

       do p = 1, 800
          call mt_random_number(px)
          call mt_random_number(py)
          call mt_random_number(a)
          a = -1.0d0 + 2.0d0 * a

          px = ctx%prob_lo(1) + pad + px * (ctx%prob_hi(1) - ctx%prob_lo(1) - 2*pad)
          py = ctx%prob_lo(2) + pad + py * (ctx%prob_hi(2) - ctx%prob_lo(2) - 2*pad)

          do j = lo(2), hi(2)
             y = ctx%prob_lo(2) + j * ctx%dx
             do i = lo(1), hi(2)
                x = ctx%prob_lo(1) + i * ctx%dx
          
                d = (x - px)**2 + (y - py)**2
                qp(i,j,1,iv) = qp(i,j,1,iv) + 0.01d0 * a * dexp(-4.0d0*d) 
             end do
          end do
       end do

    end do
  end subroutine initial

end module chemotaxis_module
