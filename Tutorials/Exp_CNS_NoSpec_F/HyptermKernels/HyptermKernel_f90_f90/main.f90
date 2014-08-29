program main

  use boxlib
  use parallel
  use multifab_module
  use advance_module
  use fabio_module

  implicit none
  integer, parameter :: DM = 3
  integer, parameter :: NG  = 4
  integer, parameter :: NC  = 5

!--------------------------------------------------------------
    integer                   :: ngrids,ib,j,k
    double precision, pointer :: dp(:,:,:,:)
    double precision :: xloc,yloc,zloc,rholoc,eloc,uvel,vvel,wvel,scale
!--------------------------------------------------------------
  integer            :: nsteps, n_cell, max_grid_size
  integer            :: un, farg, narg
  logical            :: need_inputs_file, found_inputs_file
  character(len=128) :: inputs_file_name
  integer            :: i, lo(DM), hi(DM), istep
  double precision   :: prob_lo(DM), prob_hi(DM), cfl, eta, alam
  double precision   :: dx(DM), dt, time, start_time, end_time
  logical            :: is_periodic(DM)
  type(box)          :: bx
  type(boxarray)     :: ba
  type(layout)       :: la
  type(multifab)     :: U
  double precision   :: sleep_start_time, sleep_end_time

  real(kind=dp_t) :: starts, ends

  call boxlib_initialize()

  sleep_start_time = parallel_wtime()
  call Sleep(1)
  sleep_end_time = parallel_wtime()
  if ( parallel_IOProcessor() ) then
     write(6,42),"Sleep time (s) =",sleep_end_time-sleep_start_time
  end if

  nsteps        = 1
  n_cell        = 64
  max_grid_size = 64
  cfl           = 0.5d0
  eta           = 1.8d-4
  alam          = 1.5d2

  prob_lo     = -2.3d0
  prob_hi     =  2.3d0
  is_periodic = .true.
  !
  lo = 0
  hi = n_cell-1
  bx = make_box(lo,hi)

  do i = 1,DM
     dx(i) = (prob_hi(i)-prob_lo(i)) / n_cell
     write(6,42),"dx = ", dx(i)
  end do

  call boxarray_build_bx(ba,bx)

  call boxarray_maxsize(ba,max_grid_size)

  call layout_build_ba(la,ba,boxarray_bbox(ba),pmask=is_periodic)

  call destroy(ba)

  call multifab_build(U,la,NC,NG)

!------------------------------------------------- initialize U
    ngrids = U%ng

    do ib=1,nfabs(U)
      dp => dataptr(U,ib)
      lo = lwb(get_box(U,ib))
      hi = upb(get_box(U,ib))
      scale = 1.0d0

      do k=lo(3)-ng,hi(3)+ng
         zloc = dfloat(k)*dx(3)/scale
         do j=lo(2)-ng,hi(2)+ng
            yloc = dfloat(j)*dx(2)/scale
            do i=lo(1)-ng,hi(1)+ng
               xloc = dfloat(i)*dx(1)/scale

               uvel   = 1.1d4*sin(1*xloc)*sin(2*yloc)*sin(3*zloc)
               vvel   = 1.0d4*sin(2*xloc)*sin(4*yloc)*sin(1*zloc)
               wvel   = 1.2d4*sin(3*xloc)*cos(2*yloc)*sin(2*zloc)
               rholoc = 1.0d-3 + 1.0d-5*sin(1*xloc)*cos(2*yloc)*cos(3*zloc)
               eloc   = 2.5d9  + 1.0d-3*sin(2*xloc)*cos(2*yloc)*sin(2*zloc)

               dp(i,j,k,irho) = rholoc
               dp(i,j,k,imx)  = rholoc*uvel
               dp(i,j,k,imy)  = rholoc*vvel
               dp(i,j,k,imz)  = rholoc*wvel
               dp(i,j,k,iene) = rholoc*(eloc + (uvel**2+vvel**2+wvel**2)/2)

            enddo
         enddo
      enddo
    end do

  call fabio_multifab_write_d(U, ".", "mfUInit", nOutfiles = 1)
!--------------------------------------------------------------

  istep = 0
  time  = 0.d0

  start_time = parallel_wtime()

  do istep=1,nsteps

     call advance(U,dt,dx,cfl,eta,alam)

     time = time + dt

  end do

  end_time = parallel_wtime()

  call destroy(U)
  call destroy(la)


  if ( parallel_IOProcessor() ) then
     write(6,42),"Run time (s) =",end_time-start_time
     write(6,42),"Run time per iteration (s) =",(end_time-start_time) / nsteps
     print *,"-----------------"
  end if

42  format(a,f12.8)

  call boxlib_finalize()

end program main
