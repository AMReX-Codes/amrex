subroutine orig_pingpong(sync)
  use fboxlib_mpi
  implicit none
  integer :: sync
  integer, parameter :: Ntimes=3, Npower=15, Nbase=1024
  integer, parameter :: dp_t = kind(1.0d0)
  integer :: istatus(MPI_STATUS_SIZE)
  integer :: ierr, npes, ipe, ICOMM

  real(kind=dp_t) :: sec, secmin, secmax, secavg, tar,      &
       ratemin, ratemax, rateavg, fmegabytes

  integer :: Nbytes, Nw
  real(kind=dp_t) :: tc0, tc1

  integer :: r1, s1

  real(kind=dp_t), allocatable, dimension(:) :: A, B
  integer, allocatable, dimension(:) :: nsize

  integer :: n,i,k
  logical :: lrs

  ! ****************************************************************
  ICOMM = MPI_COMM_WORLD     ! Use sorter parameter-- this if Fortran.

  allocate( nsize(0:Npower) )

  do i = 0, Npower    ! Make array of word sizes for send/recv arrays.
     nsize(i) = (Nbase*2**i)/8
  end do

  N = nsize(Npower)
  if( (Nbase .ne. nsize(0)*8) .or. (N .eq. 0) ) stop
  if( (Ntimes .lt. 3) .or. (Ntimes .gt. 99)   ) stop
  allocate( A(N), B(N) )

  call MPI_Initialized(lrs, ierr)

  if ( .not. lrs) call MPI_Init(ierr)
  call MPI_Comm_size(ICOMM, npes, ierr)
  call MPI_Comm_rank(ICOMM, ipe, ierr)

  a=ipe
  b=ipe

  tar = 0
  do i = 1, 100
     tc0 = MPI_Wtime()
     tc1 = MPI_Wtime()
     tar = (tc1-tc0) + tar
  end do
  tar = tar/100          ! Used to subtract calling time of MPI_Wtime().

  do k = 0, Npower     !Loop of number of runs.

     Nw = nsize(k)
     ! Initial time info.
     secavg = 0.0; secmin= HUGE(secmin); secmax=-HUGE(secmax) 

     do i = 1,Ntimes            !Loop number of runs.

        tc0 = MPI_Wtime()
        select case ( sync )
        case (0)
           if ( ipe == 0 ) then
              call MPI_send(a,Nw,MPI_DOUBLE_PRECISION,1,1,ICOMM,ierr)
              call MPI_recv(b,Nw,MPI_DOUBLE_PRECISION,1,0,ICOMM,istatus,ierr)
           else
              call MPI_recv(b,Nw,MPI_DOUBLE_PRECISION,0,1,ICOMM,istatus,ierr)
              call MPI_send(a,Nw,MPI_DOUBLE_PRECISION,0,0,ICOMM,ierr)
           end if
        case (1)
           if ( ipe .eq. 0 ) then
              call MPI_isend(a,Nw,MPI_DOUBLE_PRECISION,1,1,ICOMM,s1,ierr)
              call MPI_irecv(b,Nw,MPI_DOUBLE_PRECISION,1,0,ICOMM,r1,ierr)
           else
              call MPI_irecv(b,Nw,MPI_DOUBLE_PRECISION,0,1,ICOMM,r1,ierr)
              call MPI_isend(a,Nw,MPI_DOUBLE_PRECISION,0,0,ICOMM,s1,ierr)
           end if
           call MPI_Wait(r1,istatus,ierr)
           call MPI_Wait(s1,istatus,ierr)
        case (2)
           if ( ipe .eq. 0 ) then
              call MPI_irecv(b,Nw,MPI_DOUBLE_PRECISION,1,0,ICOMM,r1,ierr)
              call MPI_send(a,Nw,MPI_DOUBLE_PRECISION,1,1,ICOMM,ierr)
           else
              call MPI_irecv(b,Nw,MPI_DOUBLE_PRECISION,0,1,ICOMM,r1,ierr)
              call MPI_send(a,Nw,MPI_DOUBLE_PRECISION,0,0,ICOMM,ierr)
           end if
           call MPI_Wait(r1,istatus,ierr)
        case (3)
           if ( ipe .eq. 0 ) then
              call MPI_isend(a,Nw,MPI_DOUBLE_PRECISION,1,1,ICOMM,s1,ierr)
              call MPI_recv(b,Nw,MPI_DOUBLE_PRECISION,1,0,ICOMM,istatus,ierr)
           else
              call MPI_isend(a,Nw,MPI_DOUBLE_PRECISION,0,0,ICOMM,s1,ierr)
              call MPI_recv(b,Nw,MPI_DOUBLE_PRECISION,0,1,ICOMM,istatus,ierr)
           end if
           call MPI_Wait(s1,istatus,ierr)
        case default
           stop 'what 2'
        end select

        tc1 = MPI_Wtime()
        sec   = (tc1-tc0) - tar

        if ( ipe .eq. 0 ) then
           if ( any(b(1:Nw) /= 1) ) then
              print *, 'count = ', count(b /= 1)
              stop 'what 0'
           end if
        else
           if ( any(b(1:Nw) /= 0) ) then
              print *, 'count = ', count(b /= 0)
              stop 'what 1'
           end if
        end if

        secmin = min(sec, secmin)                ! Calc. Min, Max, & Avg.
        secmax = max(sec, secmax)
        secavg = secavg + sec

     end do

     secavg = secavg/real(Ntimes)
     Nbytes = Nw*8
     ! Determine rates.
     fmegabytes = real(Nbytes,kind=dp_t)/real(1024*1024,kind=dp_t)
     ratemax = real(2*fmegabytes)/secmin      ! 2 for round trip
     ratemin = real(2*fmegabytes)/secmax
     rateavg = real(2*fmegabytes)/secavg

     ! Report results.
     if ( ipe.eq.0 ) then
        if ( k.eq.0 ) then
           write(*,'("RUNS  BYTES",14x,"RATE(MB/s)",25x,"TIME(sec)")')
           write(*,'(15x,"MIN         AVG         MAX", 7x,&
                &       "MAX         AVG         MIN ")' )
        end if
        write(*,'(i2,1x,i8,                        &
             &             F9.4,"   ", F9.4,"   ", F9.4, 3x,&
             &             F9.6,"   ", F9.6,"   ", F9.6)' ) &
             &         Ntimes, Nbytes, ratemin, rateavg, ratemax, &
             &         secmax, secavg, secmin
     end if
  end do

  if ( .not. lrs ) call MPI_Finalize(ierr)

end subroutine orig_pingpong
subroutine my_pingpong(sync)
  use parallel
  implicit none
  integer :: sync
  integer, parameter :: Ntimes=3, Npower=15, Nbase=1024
  integer :: npes, ipe

  real(kind=dp_t) :: sec, secmin, secmax, secavg, tar,      &
       ratemin, ratemax, rateavg, fmegabytes

  integer :: Nbytes, Nw
  real(kind=dp_t) :: tc0, tc1

  integer :: r1, s1

  real(kind=dp_t), allocatable, dimension(:) :: A, B
  integer, allocatable, dimension(:) :: nsize

  integer :: n,i,k
  logical :: lrs


  ! ****************************************************************

  lrs = parallel_initialized()

  allocate( nsize(0:Npower) )

  do i = 0, Npower    ! Make array of word sizes for send/recv arrays.
     nsize(i) = (Nbase*2**i)/8
  end do

  N = nsize(Npower)
  if( (Nbase .ne. nsize(0)*8) .or. (N .eq. 0) ) stop
  if( (Ntimes .lt. 3) .or. (Ntimes .gt. 99)   ) stop
  allocate( A(N), B(N) )

! call parallel_initialize()
  if ( .not. lrs ) call parallel_initialize()
  npes = parallel_nprocs()
  ipe = parallel_myproc()

  a=ipe
  b=ipe

  tar = 0
  do i = 1, 100
     tc0 = parallel_wtime()
     tc1 = parallel_wtime()
     tar = (tc1-tc0) + tar
  end do
  tar = tar/100         ! Used to subtract calling time of MPI_Wtime().
  print *, ipe, ': tar = ', tar, 'wtick = ', parallel_wtick()

  do k = 0, Npower              ! Loop of number of runs.

     Nw = nsize(k)
     ! Initial time info.
     secavg = 0.0; secmin= HUGE(secmin); secmax=-HUGE(secmax)

     do i = 1,Ntimes            ! Loop number of runs.

        tc0 = parallel_wtime()
        select case ( sync )
        case (0)
           if ( ipe == 0 ) then
              call parallel_send_dv(a,Nw,1,1)
              call parallel_recv_dv(b,Nw,1,0)
           else
              call parallel_recv_dv(b,Nw,0,1)
              call parallel_send_dv(a,Nw,0,0)
           end if
        case (1)
           if ( ipe == 0 ) then
              s1 = parallel_isend_dv(a,Nw,1,1)
              r1 = parallel_irecv_dv(b,Nw,1,0)
           else
              r1 = parallel_irecv_dv(b,Nw,0,1)
              s1 = parallel_isend_dv(a,Nw,0,0)
           end if

           call parallel_wait(r1)
           call parallel_wait(s1)
        case (2)
           if ( ipe == 0 ) then
              r1 = parallel_irecv_dv(b,Nw,1,0)
              call parallel_send_dv(a,Nw,1,1)
           else
              r1 = parallel_irecv_dv(b,Nw,0,1)
              call parallel_send_dv(a,Nw,0,0)
           end if
           call parallel_wait(r1)
        case (3)
           if ( ipe == 0 ) then
              s1 = parallel_isend_dv(a,Nw,1,1)
              call parallel_recv_dv(b,Nw,1,0)
           else
              s1 = parallel_isend_dv(a,Nw,0,0)
              call parallel_recv_dv(b,Nw,0,1)
           end if
           call parallel_wait(s1)
        case default
           stop 'what 2'
        end select
        tc1 = parallel_wtime()
        sec   = (tc1-tc0) - tar

        if ( parallel_ioprocessor() ) then
           if ( any(b(1:Nw) /= 1) ) then
              print *, 'count = ', count(b /= 1)
              stop 'what 0'
           end if
        else
           if ( any(b(1:Nw) /= 0) ) then
              print *, 'count = ', count(b /= 0)
              stop 'what 1'
           end if
        end if

        secmin = min(sec, secmin) ! Calc. Min, Max, & Avg.
        secmax = max(sec, secmax)
        secavg = secavg + sec

     end do

     secavg = secavg/real(Ntimes)
     Nbytes = Nw*8
     ! Determine rates.
     fmegabytes = real(Nbytes,kind=dp_t)/real(1024*1024,kind=dp_t)
     ratemax = real(2*fmegabytes)/secmin ! 2 for round trip
     ratemin = real(2*fmegabytes)/secmax
     rateavg = real(2*fmegabytes)/secavg

     ! Report results.
     if ( parallel_ioprocessor() ) then
        if ( k.eq.0 ) then
           write(*,'("RUNS  BYTES",14x,"RATE(MB/s)",25x,"TIME(sec)")')
           write(*,'(15x,"MIN         AVG         MAX", 7x,&
                &       "MAX         AVG         MIN ")' )
        end if
        write(*,'(i2,1x,i8,                        &
             &             F9.4,"   ", F9.4,"   ", F9.4, 3x,&
             &             F9.6,"   ", F9.6,"   ", F9.6)' ) &
             &         Ntimes, Nbytes, ratemin, rateavg, ratemax, &
             &         secmax, secavg, secmin
     end if
  end do

  if ( .not. lrs) call parallel_finalize()

end subroutine my_pingpong

subroutine t_pingpong
  integer  :: sync = 3
  call my_pingpong(sync)
! call orig_pingpong(sync)
end subroutine t_pingpong
