subroutine t_particle

  use fabio_module
  use particle_module
  use multifab_module

  integer i, id
  type(particle_container) :: v, tv

  type(box)         :: bx,bx2
  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla
  type(boxarray)    :: ba
  double precision  :: dx(1,MAX_SPACEDIM), dx2(2,MAX_SPACEDIM)
  double precision  :: problo(MAX_SPACEDIM)
  double precision  :: probhi(MAX_SPACEDIM), time
  logical           :: pmask(MAX_SPACEDIM)
  character(len=256):: check_file_name
  character(len=5)  :: check_index

  type(multifab), allocatable :: mf(:)

  id = 1

  pmask = .true.

  problo = -1.0d0
  probhi = +1.0d0

  call particle_setverbose(.true.)

  call particle_container_setdebugging(.true.)
  !
  ! Uncomment one of the following lines.
  !
  ! Due to having some static storage in the particles code
  ! you can only run 2D or 3D but not both in the same executable.
  !
!  call test(dm=2, npart=10000)

  call test(dm=3, npart=1000)

contains

  subroutine test(dm,npart)

    integer, intent(in) :: dm
    integer, intent(in) :: npart

    type(particle), pointer :: pp(:)
    !
    ! Let's build a single level mla
    !
    call bl_assert(dm == 2 .or. dm == 3, 'only 2 or 3-D supported')

    if (parallel_IOProcessor()) then
       print*, 'Entered test() with dm = ', dm, ', npart = ', npart
    end if

    call build(v)

    if (dm == 2) then
       bx = make_box( (/0,0/), (/63,63/) )
       call build(ba,bx)
       call boxarray_maxsize(ba,32)
    else
       bx = make_box( (/0,0,0/), (/31,31,31/) )
       call build(ba,bx)
       call boxarray_maxsize(ba,16)
    end if

    if (parallel_IOProcessor()) call print(ba)

    do i = 1, dm
       dx(1,i) = (probhi(i) - problo(i)) / extent(bx,i)
    end do

    call build(mba, ba, bx)

    call destroy(ba)

    if (parallel_IOProcessor()) print*, 'pmask: ', pmask(1:dm)

    call build(mla, mba, pmask(1:dm))

    call destroy(mba)

    allocate(mf(1))

    call build(mf(1), mla%la(1), 5, 1)

    call setval(mf(1), 1.0d0)

    call init_random(v,npart,17971,mla,dx,problo,probhi)
    !
    ! A quick test of dataptr()
    !
    pp => dataptr(v)

    call bl_assert(size(pp) == capacity(v), 'dataptr does not appear to be correct')

    if (parallel_IOProcessor()) then
       print*, ''
       print*, 'size(v): ', size(v)
       print*, ''
       !     call print(v, 'after init_random')
       print*, ''
    end if

    call parallel_barrier()
    !
    ! Let's move the particles a bit.
    !
    time = 1.0d0

    do i = 1,10
       if (parallel_IOProcessor()) then
          print*, i, 'Calling move_random(one-level mla) ...'
          call flush(6)
       end if
       call move_random(v,mla,dx,problo,probhi)

       write(unit=check_index,fmt='(i5.5)') i
       check_file_name = 'chk' // check_index

       if ( parallel_IOProcessor() ) then
          call fabio_mkdir(check_file_name)
       end if
       call parallel_barrier()

       if (parallel_IOProcessor()) then
          print*, 'check_file_name: ', check_file_name
       end if

       call particle_container_checkpoint(v,check_file_name,mla)

       call build(tv)

       call particle_container_restart(tv,check_file_name,mla,dx,problo)

       call bl_assert(size(tv) == size(v), 'v and tv are NOT the same size')

       call destroy(tv)

       call timestamp(v, 'timestamp_onelev', mf, (/1,3,5/), (/"a", "b", "c"/), time)

       time = time + 0.1d0
    end do

    if (parallel_IOProcessor()) then
       print*, ''
       print*, 'size(v): ', size(v)
       print*, ''
       !     call print(v, 'after move_random')
       print*, ''
       call flush(6)
    end if

    call parallel_barrier()

!    call bl_error('Got Here')

    call destroy(mf(1))

    deallocate(mf)

    call destroy(mla)
    !
    ! Now let's try for a multi-level mla
    !
    call clear(v)

    call build(v)

    call destroy(ba)

    call build(mba,2,dm)

    mba%rr(1,1:dm) = 2

    call build(ba,bx)

    if (dm == 2) then
       call boxarray_maxsize(ba,32)
    else
       call boxarray_maxsize(ba,16)
    end if

    if (parallel_IOProcessor()) call print(ba, 'level 1')

    do i = 1, dm
       dx2(1,i) = (probhi(i) - problo(i)) / extent(bx,i)
       dx2(2,i) = dx2(1,i) / 2.0d0
    end do

    call copy(mba%bas(1),ba)
    mba%pd(1) = bx

    if (parallel_IOProcessor()) call print(mba%pd(1), 'pd(1)')

    call destroy(ba)

    if (dm == 2) then
       bx2 = make_box( (/32,32/), (/95,95/) )
       call build(ba,bx2)
       call boxarray_maxsize(ba,32)
    else
       bx2 = make_box( (/16,16,16/), (/47,47,47/) )
       call build(ba,bx2)
       call boxarray_maxsize(ba,16)
    end if

    if (parallel_IOProcessor()) call print(ba, 'level 2')

    call copy(mba%bas(2),ba)
    mba%pd(2) = refine(mba%pd(1),2)

    if (parallel_IOProcessor()) call print(mba%pd(2), 'pd(2)')

    call destroy(ba)

    call build(mla, mba, pmask(1:dm))

    call destroy(mba)

    allocate(mf(2))

    call build(mf(1), mla%la(1), 6, 2)
    call build(mf(2), mla%la(2), 6, 2)

    call setval(mf(1), 1.0d0)
    call setval(mf(2), 2.0d0)

    call init_random(v,npart,171717171,mla,dx2,problo,probhi)

    if (parallel_IOProcessor()) then
       print*, ''
       print*, 'size(v): ', size(v)
       print*, ''
       !     call print(v, 'after init_random using 2-level mla')
       print*, ''
       call flush(6)
    end if
    !
    ! Let's move the particles a bit.
    !
    do i = 1,10
       if (parallel_IOProcessor()) then
          print*, i, 'Calling move_random(two-level mla) ...'
          call flush(6)
       end if
       call move_random(v,mla,dx2,problo,probhi)

       write(unit=check_index,fmt='(i5.5)') i
       check_file_name = 'chk' // check_index

       if ( parallel_IOProcessor() ) then
          call fabio_mkdir(check_file_name)
       end if
       call parallel_barrier()

       if (parallel_IOProcessor()) then
          print*, 'check_file_name: ', check_file_name
       end if

       call particle_container_checkpoint(v,check_file_name,mla)

       call build(tv)

       call particle_container_restart(tv,check_file_name,mla,dx2,problo)

       call bl_assert(size(tv) == size(v), 'v and tv are NOT the same size')

       call destroy(tv)

       call timestamp(v, 'timestamp_twolev', mf, (/1,3,5/), (/"a", "b", "c"/), time)

       time = time + 0.1d0
    end do

    call destroy(mf(1))
    call destroy(mf(2))

    deallocate(mf)

    call destroy(mla)

    call destroy(v)

  end subroutine test

end subroutine t_particle
