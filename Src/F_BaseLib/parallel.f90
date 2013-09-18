
      module fboxlib_mpi
      implicit none
      include 'mpif.h'
      end module fboxlib_mpi

!! MPI wrappers
!! These wrappers are used so that a non-MPI version can coexist
!! with the MPI version.
!! Not all of MPI is supported.  Only sufficient to execute
!! fBoxLib calls.

module parallel

  ! Assumption:
  ! 1) The user has not replaced the default error handler which
  !    is MPI_ERRORS_ARE_FATAL

  use bl_types; use fboxlib_mpi

  implicit none

  integer, parameter :: parallel_root = 0
  integer, parameter, private :: io_processor_node = parallel_root
  integer, private :: m_nprocs = -1
  integer, private :: m_myproc = -1
  integer, private :: m_comm   = -1
  integer, private :: m_thread_support_level = MPI_THREAD_SINGLE

  interface parallel_reduce
     module procedure parallel_reduce_d
     module procedure parallel_reduce_r
     module procedure parallel_reduce_i
     module procedure parallel_reduce_l
     module procedure parallel_reduce_dv
     module procedure parallel_reduce_rv
     module procedure parallel_reduce_iv
     module procedure parallel_reduce_lv
  end interface parallel_reduce

!   interface parallel_isend
!      module procedure parallel_isend_dv
!      module procedure parallel_isend_rv
!      module procedure parallel_isend_iv
!      module procedure parallel_isend_lv
!   end interface parallel_isend

!   interface parallel_irecv
!      module procedure parallel_irecv_dv
!      module procedure parallel_irecv_rv
!      module procedure parallel_irecv_iv
!      module procedure parallel_irecv_lv
!   end interface parallel_irecv

  interface parallel_send
     module procedure parallel_send_d1
     module procedure parallel_send_d2
     module procedure parallel_send_d3
     module procedure parallel_send_d4
     module procedure parallel_send_d5
     module procedure parallel_send_d6
     module procedure parallel_send_d7
     module procedure parallel_send_i1
     module procedure parallel_send_i2
     module procedure parallel_send_i3
     module procedure parallel_send_i4
     module procedure parallel_send_i5
     module procedure parallel_send_i6
     module procedure parallel_send_i7
     module procedure parallel_send_l1
     module procedure parallel_send_l2
     module procedure parallel_send_l3
     module procedure parallel_send_l4
     module procedure parallel_send_l5
     module procedure parallel_send_l6
     module procedure parallel_send_l7
     module procedure parallel_send_z1
     module procedure parallel_send_z2
     module procedure parallel_send_z3
     module procedure parallel_send_z4
     module procedure parallel_send_z5
     module procedure parallel_send_z6
     module procedure parallel_send_z7
  end interface parallel_send

  interface parallel_recv
     module procedure parallel_recv_d1
     module procedure parallel_recv_d2
     module procedure parallel_recv_d3
     module procedure parallel_recv_d4
     module procedure parallel_recv_d5
     module procedure parallel_recv_d6
     module procedure parallel_recv_d7
     module procedure parallel_recv_i1
     module procedure parallel_recv_i2
     module procedure parallel_recv_i3
     module procedure parallel_recv_i4
     module procedure parallel_recv_i5
     module procedure parallel_recv_i6
     module procedure parallel_recv_i7
     module procedure parallel_recv_l1
     module procedure parallel_recv_l2
     module procedure parallel_recv_l3
     module procedure parallel_recv_l4
     module procedure parallel_recv_l5
     module procedure parallel_recv_l6
     module procedure parallel_recv_l7
     module procedure parallel_recv_z1
     module procedure parallel_recv_z2
     module procedure parallel_recv_z3
     module procedure parallel_recv_z4
     module procedure parallel_recv_z5
     module procedure parallel_recv_z6
     module procedure parallel_recv_z7
  end interface parallel_recv

  interface parallel_wait
     module procedure parallel_wait_one
     module procedure parallel_wait_vec
     module procedure parallel_wait_vec_vec
     module procedure parallel_wait_vec_vec_vec
  end interface parallel_wait

  interface parallel_test
     module procedure parallel_test_one
     module procedure parallel_test_vec
     module procedure parallel_test_vec_vec
     module procedure parallel_test_vec_vec_vec
  end interface parallel_test

  interface parallel_bcast
     module procedure parallel_bcast_d
     module procedure parallel_bcast_r
     module procedure parallel_bcast_i
     module procedure parallel_bcast_l
     module procedure parallel_bcast_c
     module procedure parallel_bcast_z
     module procedure parallel_bcast_dv
     module procedure parallel_bcast_rv
     module procedure parallel_bcast_iv
     module procedure parallel_bcast_lv
     module procedure parallel_bcast_cv
     module procedure parallel_bcast_zv
  end interface parallel_bcast

  interface parallel_scatter
     module procedure parallel_scatter_dv
     module procedure parallel_scatter_rv
     module procedure parallel_scatter_iv
     module procedure parallel_scatter_lv
     module procedure parallel_scatter_cv
     module procedure parallel_scatter_zv
  end interface parallel_scatter

  interface parallel_gather
     !
     ! Gather fixed size blocks to specified processor
     !
     module procedure parallel_gather_d
     module procedure parallel_gather_r
     module procedure parallel_gather_i
     module procedure parallel_gather_l
     module procedure parallel_gather_c
     module procedure parallel_gather_z
     !
     ! Gather variable sized blocks to specified processor
     !
     module procedure parallel_gather_dv
     module procedure parallel_gather_rv
     module procedure parallel_gather_iv
     module procedure parallel_gather_lv
     module procedure parallel_gather_cv
     module procedure parallel_gather_zv
  end interface parallel_gather

  interface parallel_allgather
     module procedure parallel_allgather_dv
     module procedure parallel_allgather_rv
     module procedure parallel_allgather_iv
     module procedure parallel_allgather_lv
     module procedure parallel_allgather_cv
     module procedure parallel_allgather_zv
  end interface parallel_allgather

  interface parallel_alltoall
     module procedure parallel_alltoall_d
     module procedure parallel_alltoall_dv
     module procedure parallel_alltoall_i
     module procedure parallel_alltoall_iv
     module procedure parallel_alltoall_l
     module procedure parallel_alltoall_lv
  end interface

contains

  function parallel_q() result(r)
    logical :: r
    r = .TRUE.
  end function parallel_q

  function parallel_initialized() result(r)
    logical :: r
    integer :: ierr
    external MPI_Initialized
    call MPI_Initialized(r, ierr)
  end function parallel_initialized

  subroutine parallel_initialize(comm, thread_support_level)
    integer, intent(in), optional :: comm, thread_support_level
    integer ierr
    logical flag
    external MPI_Init, MPI_Comm_Dup, MPI_Comm_Size, MPI_Comm_Rank
    call MPI_Initialized(flag, ierr)
    if ( .not. flag ) then
       if (present(thread_support_level)) then
          call MPI_Init_thread(thread_support_level, m_thread_support_level, ierr)
       else
          call MPI_Init(ierr)
       end if
    end if
    if ( present(comm) ) then
       call MPI_Comm_Dup(comm, m_comm, ierr)
    else
       call MPI_Comm_Dup(MPI_COMM_WORLD, m_comm, ierr)
    endif
    call MPI_Comm_Size(m_comm, m_nprocs, ierr)
    call MPI_Comm_Rank(m_comm, m_myproc, ierr)
    call parallel_barrier()
  end subroutine parallel_initialize
  subroutine parallel_finalize(do_finalize_MPI)
    logical, intent(in), optional :: do_finalize_MPI
    integer ierr
    external MPI_Comm_Free, MPI_Finalize
    !call MPI_Comm_Free(m_comm, ierr)  !Note: This is *supposed* to be the right way to do this, but it crashes on Linux.  comment out leads to small mem leak
    m_comm = MPI_COMM_WORLD
    if (present(do_finalize_MPI) ) then
       if (do_finalize_MPI) call MPI_Finalize(ierr)
    else
       call MPI_Finalize(ierr)
    endif
    
  end subroutine parallel_finalize

  subroutine parallel_abort(str)
    character(len=*), optional :: str
    external MPI_Abort
    integer :: ierr
    if ( parallel_IOProcessor() ) then
       if ( present(str) ) then
          print*, 'parallel_abort(): ', str
       else
          print*, 'parallel_abort() !!!'
       end if
    end if
    call MPI_Abort(m_comm, -1, ierr)
  end subroutine parallel_abort

  subroutine parallel_set_comm(comm)
    integer, intent(in) :: comm
    integer :: ierr
    m_comm = comm
    call MPI_Comm_Size(m_comm, m_nprocs, ierr)
    call MPI_Comm_Rank(m_comm, m_myproc, ierr)
    call parallel_barrier()
  end subroutine parallel_set_comm
  !
  ! Returns the null communicator.
  !
  function parallel_null_communicator() result(comm)
    integer :: comm
    comm = MPI_COMM_NULL
  end function parallel_null_communicator
  !
  ! Create a new communicator from the set of unique proc IDs in procs.
  ! The proc IDs in procs do not need to be unique.  They just need to
  ! constitute a subset of the proc IDs in MPI_COM_WORLD.  One possible
  ! way to call this would be with the processor map from a layout.
  ! Note that this will return MPI_COMM_NULL to those MPI procs
  ! that are not in the new communicator.
  !
  function parallel_create_communicator(procs) result(comm)
    use sort_i_module
    use vector_i_module

    integer              :: comm
    integer, intent(in)  :: procs(:)

    integer              :: i, world_group, this_group, ierr
    integer, allocatable :: ranks(:)
    type(vector_i)       :: v
    !
    ! Make sure all possible procs are in the current MPI_COMM_WORLD.
    !
    do i = 1, size(procs)
       if ( (procs(i) .lt. 0) .or. (procs(i) .ge. parallel_nprocs()) ) then
          if ( parallel_IOProcessor() ) then
             print*, 'procs must be in range: [0,nprocs)'
             call flush(6)
          end if
          call parallel_abort()
       end if
    end do

    allocate(ranks(size(procs)))

    ranks = procs
    !
    ! Sort & remove duplicates from "rank".
    !
    call sort(ranks)
    call build(v)
    call reserve(v,parallel_nprocs())
    call push_back(v,ranks(1))
    do i = 2, size(ranks)
       if ( ranks(i) .ne. back(v) ) call push_back(v,ranks(i))
    end do
    !
    ! Build a duplicate of the MPI_COMM_WORLD group.
    !
    call MPI_Comm_group(MPI_COMM_WORLD, world_group, ierr)

    call MPI_group_incl(world_group, size(v), dataptr(v), this_group, ierr)

    call destroy(v)

    deallocate(ranks)
    !
    ! This sets comm to MPI_COMM_NULL on those ranks not in this_group.
    !
    call MPI_Comm_create(MPI_COMM_WORLD, this_group, comm, ierr)

    call MPI_Group_free(this_group,  ierr)
    call MPI_Group_free(world_group, ierr)

  end function parallel_create_communicator
  
  subroutine parallel_free_communicator(comm)
    integer :: comm, ierr
    if ( comm .ne. MPI_COMM_NULL ) call MPI_Comm_free(comm, ierr)
  end subroutine parallel_free_communicator

  pure function parallel_communicator() result(r)
    integer :: r
    r = m_comm
  end function parallel_communicator
  pure function parallel_nprocs() result(r)
    integer r
    r = m_nprocs
  end function parallel_nprocs
  pure function parallel_myproc() result(r)
    integer r
    r = m_myproc
  end function parallel_myproc
  function parallel_IOProcessor(comm) result(r)
    logical :: r
    integer, intent(in), optional :: comm
    integer :: rank, ierr
    if ( present(comm) ) then
       if ( comm .eq. MPI_COMM_NULL ) then
          r = .false.
       else
          call MPI_Comm_rank(comm, rank, ierr)
          r = (rank == io_processor_node)
       end if
    else
       r = (m_myproc == io_processor_node)
    end if
  end function parallel_IOProcessor
  pure function parallel_IOProcessorNode() result(r)
    integer :: r
    r = io_processor_node
  end function parallel_IOProcessorNode
  pure function parallel_thread_support_level() result(r)
    integer :: r
    r = m_thread_support_level
  end function parallel_thread_support_level
  function parallel_wtime() result(r)
    real(kind=dp_t) :: r
    r = MPI_Wtime()
  end function parallel_wtime
  function parallel_wtick() result(r)
    real(kind=dp_t) :: r
    r = MPI_Wtick()
  end function parallel_wtick

  ! REDUCE:
  subroutine parallel_reduce_d(r, a, op, proc, comm)
    real(kind=dp_t), intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc
    integer, intent(in), optional :: comm
    real(kind=dp_t), intent(out) :: r
    integer ierr
    external MPI_Reduce, MPI_AllReduce
    integer :: l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    if ( present(proc) ) then
       CALL MPI_Reduce(a, r, 1, MPI_DOUBLE_PRECISION, op, proc, l_comm, ierr) 
    else
       call MPI_AllReduce(a, r, 1, MPI_DOUBLE_PRECISION, op, l_comm, ierr)
    end if
  end subroutine parallel_reduce_d
  subroutine parallel_reduce_r(r, a, op, proc, comm)
    real(kind=sp_t), intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc
    integer, intent(in), optional :: comm
    real(kind=sp_t), intent(out) :: r
    integer ierr
    external MPI_Reduce, MPI_AllReduce
    integer :: l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    if ( present(proc) ) then
       CALL MPI_Reduce(a, r, 1, MPI_REAL, op, proc, l_comm, ierr) 
    else
       call MPI_AllReduce(a, r, 1, MPI_REAL, op, l_comm, ierr)
    end if
  end subroutine parallel_reduce_r
  subroutine parallel_reduce_i(r, a, op, proc, comm)
    integer, intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc
    integer, intent(in), optional :: comm
    integer, intent(out) :: r
    integer ierr
    external MPI_Reduce, MPI_AllReduce
    integer :: l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    if ( present(proc) ) then
       CALL MPI_Reduce(a, r, 1, MPI_INTEGER, op, proc, l_comm, ierr) 
    else
       call MPI_AllReduce(a, r, 1, MPI_INTEGER, op, l_comm, ierr)
    end if
  end subroutine parallel_reduce_i
  subroutine parallel_reduce_l(r, a, op, proc, comm)
    logical, intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc
    integer, intent(in), optional :: comm
    logical, intent(out) :: r
    integer ierr
    external MPI_Reduce, MPI_AllReduce
    integer :: l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    if ( present(proc) ) then
       CALL MPI_Reduce(a, r, 1, MPI_LOGICAL, op, proc, l_comm, ierr) 
    else
       call MPI_AllReduce(a, r, 1, MPI_LOGICAL, op, l_comm, ierr)
    end if
  end subroutine parallel_reduce_l

  ! REDUCE_V:
  subroutine parallel_reduce_dv(r, a, op, proc, comm)
    real(kind=dp_t), intent(in) :: a(:)
    integer, intent(in) :: op
    integer, intent(in), optional :: proc
    integer, intent(in), optional :: comm
    real(kind=dp_t), intent(out) :: r(:)
    integer ierr
    external MPI_Reduce, MPI_AllReduce
    integer :: l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    if ( present(proc) ) then
       CALL MPI_Reduce(a, r, size(a), MPI_DOUBLE_PRECISION, op, proc, l_comm, ierr) 
    else
       call MPI_AllReduce(a, r, size(a), MPI_DOUBLE_PRECISION, op, l_comm, ierr)
    end if
  end subroutine parallel_reduce_dv
  subroutine parallel_reduce_rv(r, a, op, proc, comm)
    real(kind=sp_t), intent(in) :: a(:)
    integer, intent(in) :: op
    integer, intent(in), optional :: proc
    integer, intent(in), optional :: comm
    real(kind=sp_t), intent(out) :: r(:)
    integer ierr
    external MPI_Reduce, MPI_AllReduce
    integer :: l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    if ( present(proc) ) then
       CALL MPI_Reduce(a, r, size(a), MPI_REAL, op, proc, l_comm, ierr) 
    else
       call MPI_AllReduce(a, r, size(a), MPI_REAL, op, l_comm, ierr)
    end if
  end subroutine parallel_reduce_rv
  subroutine parallel_reduce_iv(r, a, op, proc, comm)
    integer, intent(in) :: a(:)
    integer, intent(in) :: op
    integer, intent(in), optional :: proc
    integer, intent(in), optional :: comm
    integer, intent(out) :: r(:)
    integer ierr
    external MPI_Reduce, MPI_AllReduce
    integer :: l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    if ( present(proc) ) then
       CALL MPI_Reduce(a, r, size(a), MPI_INTEGER, op, proc, l_comm, ierr) 
    else
       call MPI_AllReduce(a, r, size(a), MPI_INTEGER, op, l_comm, ierr)
    end if
  end subroutine parallel_reduce_iv
  subroutine parallel_reduce_lv(r, a, op, proc, comm)
    logical, intent(in) :: a(:)
    integer, intent(in) :: op
    integer, intent(in), optional :: proc
    integer, intent(in), optional :: comm
    logical, intent(out) :: r(:)
    integer ierr
    external MPI_Reduce, MPI_AllReduce
    integer :: l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    if ( present(proc) ) then
       CALL MPI_Reduce(a, r, size(a), MPI_LOGICAL, op, proc, l_comm, ierr) 
    else
       call MPI_AllReduce(a, r, size(a), MPI_LOGICAL, op, l_comm, ierr)
    end if
  end subroutine parallel_reduce_lv



  ! REDUCETION operator creation:
  function parallel_op_create_d(fcn, commute) result(r)
    integer :: r
    logical, intent(in), optional :: commute
    interface
       subroutine fcn(o, i, n, t)
         use bl_types, only : dp_t
         integer, intent(in) :: n
         real(kind=dp_t) :: o(n)
         real(kind=dp_t) :: i(n)
         integer, intent(in) :: t
       end subroutine fcn
    end interface
    integer :: ierr
    logical :: l_commute
    external MPI_Op_Create
    l_commute = .TRUE.; if ( present(commute) ) l_commute = commute
    CALL MPI_Op_Create(fcn, l_commute, r, ierr)
  end function parallel_op_create_d
  subroutine parallel_op_free(op)
    integer, intent(in) :: op
    integer :: ierr
    external MPI_Op_Free
    CALL MPI_Op_Free(op, ierr)
  end subroutine parallel_op_free


  ! ISEND:
  ! non-blocking send.
  function parallel_isend_dv(a, n, proc, tag, comm) result(r)
    integer :: r
    real(kind=dp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_ISend
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_ISend(a, n, MPI_DOUBLE_PRECISION, proc, tag, l_comm, r, ierr)
  end function parallel_isend_dv
  function parallel_isend_rv(a, n, proc, tag, comm) result(r)
    integer :: r
    real(kind=sp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_ISend
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_ISend(a, n, MPI_REAL, proc, tag, l_comm, r, ierr)
  end function parallel_isend_rv
  function parallel_isend_iv(a, n, proc, tag, comm) result(r)
    integer :: r
    integer, intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_ISend
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_ISend(a, n, MPI_INTEGER, proc, tag, l_comm, r, ierr)
  end function parallel_isend_iv
  function parallel_isend_lv(a, n, proc, tag, comm) result(r)
    integer :: r
    logical, intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_ISend
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_ISend(a, n, MPI_LOGICAL, proc, tag, l_comm, r, ierr)
  end function parallel_isend_lv
  function parallel_isend_zv(a, n, proc, tag, comm) result(r)
    integer :: r
    complex(dp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_ISend
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_ISend(a, n, MPI_DOUBLE_COMPLEX, proc, tag, l_comm, r, ierr)
  end function parallel_isend_zv

  ! IRECV:
  ! non-blocking receive.
  function parallel_irecv_dv(a, n, proc, tag, comm) result(r)
    integer :: r
    real(kind=dp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_IRecv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_IRecv(a, n, MPI_DOUBLE_PRECISION, proc, tag, l_comm, r, ierr)
  end function parallel_irecv_dv
  function parallel_irecv_rv(a, n, proc, tag, comm) result(r)
    integer :: r
    real(kind=sp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_IRecv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_IRecv(a, n, MPI_REAL, proc, tag, l_comm, r, ierr)
  end function parallel_irecv_rv
  function parallel_irecv_iv(a, n, proc, tag, comm) result(r)
    integer :: r
    integer, intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_IRecv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_IRecv(a, n, MPI_INTEGER, proc, tag, l_comm, r, ierr)
  end function parallel_irecv_iv
  function parallel_irecv_lv(a, n, proc, tag, comm) result(r)
    integer :: r
    logical, intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_IRecv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_IRecv(a, n, MPI_LOGICAL, proc, tag, l_comm, r, ierr)
  end function parallel_irecv_lv
  function parallel_irecv_zv(a, n, proc, tag, comm) result(r)
    integer :: r
    complex(dp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_IRecv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_IRecv(a, n, MPI_DOUBLE_COMPLEX, proc, tag, l_comm, r, ierr)
  end function parallel_irecv_zv



  ! WAIT:
  ! completes the isend/iwait calls
  subroutine parallel_wait_one(req, status)
    integer, intent(inout)  :: req
    integer, intent(out), optional :: status(MPI_STATUS_SIZE)
    integer :: ierr, lstatus(MPI_STATUS_SIZE)
    external MPI_Wait
    call MPI_Wait(req, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_wait_one

  subroutine parallel_wait_vec(req, index, status)
    integer, intent(inout)  :: req(:)
    integer, intent(out) :: index
    integer, intent(out), optional :: status(MPI_STATUS_SIZE)
    integer :: ierr, lstatus(MPI_STATUS_SIZE)
    external MPI_WaitAny
    call MPI_WaitAny(size(req), req, index, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_wait_vec

  subroutine parallel_wait_vec_vec(req, status)
    integer, intent(inout)  :: req(:)
    integer, intent(out), optional :: status(:,:)
    integer :: ierr, lstatus(MPI_STATUS_SIZE,size(req))
    external MPI_WaitAll
    call MPI_WaitAll(size(req), req, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_wait_vec_vec

  subroutine parallel_wait_vec_vec_vec(r, req, indicies, status)
    integer, intent(out) :: r
    integer, intent(inout)  :: req(:)
    integer, intent(out) :: indicies(:)
    integer, intent(out), optional :: status(:,:)
    integer :: ierr, lstatus(MPI_STATUS_SIZE,size(req))
    external MPI_WaitSome
    call MPI_WaitSome(size(req), req, r, indicies, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_wait_vec_vec_vec


  ! TEST:
  ! check for completion of non-blocking send:
  function parallel_test_one(req, status) result(r)
    integer, intent(inout)  :: req
    logical :: r
    integer, intent(out), optional :: status(MPI_STATUS_SIZE)
    integer :: ierr, lstatus(MPI_STATUS_SIZE)
    external MPI_Test
    call MPI_Test(req, r, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end function parallel_test_one

  function parallel_test_vec(req, index, status) result(r)
    integer, intent(inout)  :: req(:)
    integer, intent(out) :: Index
    logical :: r
    integer, intent(out), optional :: status(MPI_STATUS_SIZE)
    integer :: ierr, lstatus(MPI_STATUS_SIZE)
    external MPI_TestAny
    call MPI_TestAny(size(req), req, index, r, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end function parallel_test_vec

  function parallel_test_vec_vec(req, status) result(r)
    integer, intent(inout)  :: req(:)
    logical :: r
    integer, intent(out), optional :: status(:,:)
    integer :: ierr, lstatus(MPI_STATUS_SIZE,size(req))
    external MPI_TestAll
    call MPI_TestAll(size(req), req, r, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end function parallel_test_vec_vec

  function parallel_test_vec_vec_vec(req, indicies, status) result(r)
    integer, intent(inout)  :: req(:)
    integer :: r
    integer, intent(out) :: indicies(:)
    integer, intent(out), optional :: status(:,:)
    integer :: ierr, lstatus(MPI_STATUS_SIZE,size(req))
    external MPI_TestSome
    call MPI_TestSome(size(req), req, r, indicies, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end function parallel_test_vec_vec_vec



  ! SEND:
  ! Blocking Send.
  subroutine parallel_send_d1(a, proc, tag, comm)
    real(kind=dp_t), intent(in) :: a(:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_dv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_d1
  subroutine parallel_send_d2(a, proc, tag, comm)
    real(kind=dp_t), intent(in) :: a(:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_dv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_d2
  subroutine parallel_send_d3(a, proc, tag, comm)
    real(kind=dp_t), intent(in) :: a(:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_dv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_d3
  subroutine parallel_send_d4(a, proc, tag, comm)
    real(kind=dp_t), intent(in) :: a(:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_dv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_d4
  subroutine parallel_send_d5(a, proc, tag, comm)
    real(kind=dp_t), intent(in) :: a(:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_dv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_d5
  subroutine parallel_send_d6(a, proc, tag, comm)
    real(kind=dp_t), intent(in) :: a(:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_dv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_d6
  subroutine parallel_send_d7(a, proc, tag, comm)
    real(kind=dp_t), intent(in) :: a(:,:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_dv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_d7
  subroutine parallel_send_dv(a, n, proc, tag, comm)
    real(kind=dp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_Send
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Send(a, n, MPI_DOUBLE_PRECISION, proc, tag, l_comm, ierr)
  end subroutine parallel_send_dv
  subroutine parallel_send_rv(a, n, proc, tag, comm)
    real(kind=sp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_Send
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Send(a, n, MPI_REAL, proc, tag, l_comm, ierr)
  end subroutine parallel_send_rv
  subroutine parallel_send_i1(a, proc, tag, comm)
    integer, intent(in) :: a(:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_iv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_i1
  subroutine parallel_send_i2(a, proc, tag, comm)
    integer, intent(in) :: a(:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_iv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_i2
  subroutine parallel_send_i3(a, proc, tag, comm)
    integer, intent(in) :: a(:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_iv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_i3
  subroutine parallel_send_i4(a, proc, tag, comm)
    integer, intent(in) :: a(:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_iv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_i4
  subroutine parallel_send_i5(a, proc, tag, comm)
    integer, intent(in) :: a(:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_iv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_i5
  subroutine parallel_send_i6(a, proc, tag, comm)
    integer, intent(in) :: a(:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_iv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_i6
  subroutine parallel_send_i7(a, proc, tag, comm)
    integer, intent(in) :: a(:,:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_iv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_i7
  subroutine parallel_send_iv(a, n, proc, tag, comm)
    integer, intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_Send
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Send(a, n, MPI_INTEGER, proc, tag, l_comm, ierr)
  end subroutine parallel_send_iv
  subroutine parallel_send_l1(a, proc, tag, comm)
    logical, intent(in) :: a(:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_lv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_l1
  subroutine parallel_send_l2(a, proc, tag, comm)
    logical, intent(in) :: a(:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_lv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_l2
  subroutine parallel_send_l3(a, proc, tag, comm)
    logical, intent(in) :: a(:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_lv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_l3
  subroutine parallel_send_l4(a, proc, tag, comm)
    logical, intent(in) :: a(:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_lv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_l4
  subroutine parallel_send_l5(a, proc, tag, comm)
    logical, intent(in) :: a(:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_lv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_l5
  subroutine parallel_send_l6(a, proc, tag, comm)
    logical, intent(in) :: a(:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_lv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_l6
  subroutine parallel_send_l7(a, proc, tag, comm)
    logical, intent(in) :: a(:,:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_lv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_l7
  subroutine parallel_send_lv(a, n, proc, tag, comm)
    logical, intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_Send
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Send(a, n, MPI_LOGICAL, proc, tag, l_comm, ierr)
  end subroutine parallel_send_lv
  subroutine parallel_send_z1(a, proc, tag, comm)
    complex(kind=dp_t), intent(in) :: a(:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_zv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_z1
  subroutine parallel_send_z2(a, proc, tag, comm)
    complex(kind=dp_t), intent(in) :: a(:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_zv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_z2
  subroutine parallel_send_z3(a, proc, tag, comm)
    complex(kind=dp_t), intent(in) :: a(:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_zv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_z3
  subroutine parallel_send_z4(a, proc, tag, comm)
    complex(kind=dp_t), intent(in) :: a(:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_zv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_z4
  subroutine parallel_send_z5(a, proc, tag, comm)
    complex(kind=dp_t), intent(in) :: a(:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_zv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_z5
  subroutine parallel_send_z6(a, proc, tag, comm)
    complex(kind=dp_t), intent(in) :: a(:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_zv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_z6
  subroutine parallel_send_z7(a, proc, tag, comm)
    complex(kind=dp_t), intent(in) :: a(:,:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_zv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_z7
  subroutine parallel_send_zv(a, n, proc, tag, comm)
    complex(kind=dp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_Send
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Send(a, n, MPI_DOUBLE_COMPLEX, proc, tag, l_comm, ierr)
  end subroutine parallel_send_zv



  ! RECV:
  ! Blocking Receive.
  subroutine parallel_recv_d1(a, proc, tag, comm, status)
    real(kind=dp_t), intent(out) :: a(:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_dv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_d1
  subroutine parallel_recv_d2(a, proc, tag, comm, status)
    real(kind=dp_t), intent(out) :: a(:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_dv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_d2
  subroutine parallel_recv_d3(a, proc, tag, comm, status)
    real(kind=dp_t), intent(out) :: a(:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_dv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_d3
  subroutine parallel_recv_d4(a, proc, tag, comm, status)
    real(kind=dp_t), intent(out) :: a(:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_dv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_d4
  subroutine parallel_recv_d5(a, proc, tag, comm, status)
    real(kind=dp_t), intent(out) :: a(:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_dv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_d5
  subroutine parallel_recv_d6(a, proc, tag, comm, status)
    real(kind=dp_t), intent(out) :: a(:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_dv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_d6
  subroutine parallel_recv_d7(a, proc, tag, comm, status)
    real(kind=dp_t), intent(out) :: a(:,:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_dv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_d7
  subroutine parallel_recv_dv(a, n, proc, tag, comm, status)
    real(kind=dp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    external MPI_Recv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Recv(a, n, MPI_DOUBLE_PRECISION, proc, tag, l_comm, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_dv
  subroutine parallel_recv_rv(a, n, proc, tag, comm, status)
    real(kind=sp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    external MPI_Recv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Recv(a, n, MPI_REAL, proc, tag, l_comm, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_rv
  subroutine parallel_recv_i1(a, proc, tag, comm, status)
    integer, intent(out) :: a(:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_iv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_i1
  subroutine parallel_recv_i2(a, proc, tag, comm, status)
    integer, intent(out) :: a(:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_iv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_i2
  subroutine parallel_recv_i3(a, proc, tag, comm, status)
    integer, intent(out) :: a(:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_iv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_i3
  subroutine parallel_recv_i4(a, proc, tag, comm, status)
    integer, intent(out) :: a(:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_iv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_i4
  subroutine parallel_recv_i5(a, proc, tag, comm, status)
    integer, intent(out) :: a(:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_iv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_i5
  subroutine parallel_recv_i6(a, proc, tag, comm, status)
    integer, intent(out) :: a(:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_iv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_i6
  subroutine parallel_recv_i7(a, proc, tag, comm, status)
    integer, intent(out) :: a(:,:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_iv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_i7
  subroutine parallel_recv_iv(a, n, proc, tag, comm, status)
    integer, intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    external MPI_Recv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Recv(a, n, MPI_INTEGER, proc, tag, l_comm, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_iv
  subroutine parallel_recv_l1(a, proc, tag, comm, status)
    logical, intent(out) :: a(:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_lv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_l1
  subroutine parallel_recv_l2(a, proc, tag, comm, status)
    logical, intent(out) :: a(:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_lv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_l2
  subroutine parallel_recv_l3(a, proc, tag, comm, status)
    logical, intent(out) :: a(:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_lv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_l3
  subroutine parallel_recv_l4(a, proc, tag, comm, status)
    logical, intent(out) :: a(:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_lv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_l4
  subroutine parallel_recv_l5(a, proc, tag, comm, status)
    logical, intent(out) :: a(:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_lv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_l5
  subroutine parallel_recv_l6(a, proc, tag, comm, status)
    logical, intent(out) :: a(:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_lv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_l6
  subroutine parallel_recv_l7(a, proc, tag, comm, status)
    logical, intent(out) :: a(:,:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_lv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_l7
  subroutine parallel_recv_lv(a, n, proc, tag, comm, status)
    logical, intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    external MPI_Recv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Recv(a, n, MPI_LOGICAL, proc, tag, l_comm, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_lv
  subroutine parallel_recv_z1(a, proc, tag, comm, status)
    complex(kind=dp_t), intent(out) :: a(:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_zv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_z1
  subroutine parallel_recv_z2(a, proc, tag, comm, status)
    complex(kind=dp_t), intent(out) :: a(:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_zv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_z2
  subroutine parallel_recv_z3(a, proc, tag, comm, status)
    complex(kind=dp_t), intent(out) :: a(:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_zv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_z3
  subroutine parallel_recv_z4(a, proc, tag, comm, status)
    complex(kind=dp_t), intent(out) :: a(:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_zv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_z4
  subroutine parallel_recv_z5(a, proc, tag, comm, status)
    complex(kind=dp_t), intent(out) :: a(:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_zv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_z5
  subroutine parallel_recv_z6(a, proc, tag, comm, status)
    complex(kind=dp_t), intent(out) :: a(:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_zv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_z6
  subroutine parallel_recv_z7(a, proc, tag, comm, status)
    complex(kind=dp_t), intent(out) :: a(:,:,:,:,:,:,:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_zv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_z7
  subroutine parallel_recv_zv(a, n, proc, tag, comm, status)
    complex(kind=dp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    external MPI_Recv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Recv(a, n, MPI_DOUBLE_COMPLEX, proc, tag, l_comm, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_zv

  ! AlltoAll
  subroutine parallel_alltoall_d(a, b, n, comm)
    real(kind=dp_t), intent(in) :: b(*)
    real(kind=dp_t), intent(inout) :: a(*)
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    integer, intent(in) :: n
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    call MPI_Alltoall(b, n, MPI_DOUBLE_PRECISION, a, n, MPI_DOUBLE_PRECISION, l_comm, ierr)
  end subroutine parallel_alltoall_d
  ! AlltoAll
  subroutine parallel_alltoall_dv(a, ac, ad, b, bc, bd, comm)
    real(kind=dp_t), intent(in) :: b(*)
    real(kind=dp_t), intent(inout) :: a(*)
    integer, intent(in) :: ac(*), ad(*), bc(*), bd(*)
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    call MPI_Alltoallv(b, bc, bd, MPI_DOUBLE_PRECISION, a, ac, ad, MPI_DOUBLE_PRECISION, l_comm, ierr)
  end subroutine parallel_alltoall_dv

  subroutine parallel_alltoall_i(a, b, n, comm)
    integer, intent(in) :: b(*)
    integer, intent(inout) :: a(*)
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    integer, intent(in) :: n
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    call MPI_Alltoall(b, n, MPI_INTEGER, a, n, MPI_INTEGER, l_comm, ierr)
  end subroutine parallel_alltoall_i
  ! AlltoAll
  subroutine parallel_alltoall_iv(a, ac, ad, b, bc, bd, comm)
    integer, intent(in) :: b(*)
    integer, intent(inout) :: a(*)
    integer, intent(in) :: ac(*), ad(*), bc(*), bd(*)
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    call MPI_Alltoallv(b, bc, bd, MPI_INTEGER, a, ac, ad, MPI_INTEGER, l_comm, ierr)
  end subroutine parallel_alltoall_iv

  ! AlltoAll
  subroutine parallel_alltoall_l(a, b, n, comm)
    logical, intent(in) :: b(*)
    logical, intent(inout) :: a(*)
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    integer, intent(in) :: n
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    call MPI_Alltoall(b, n, MPI_LOGICAL, a, n, MPI_LOGICAL, l_comm, ierr)
  end subroutine parallel_alltoall_l
  ! AlltoAll
  subroutine parallel_alltoall_lv(a, ac, ad, b, bc, bd, comm)
    logical, intent(in) :: b(*)
    logical, intent(inout) :: a(*)
    integer, intent(in) :: ac(*), ad(*), bc(*), bd(*)
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    call MPI_Alltoallv(b, bc, bd, MPI_LOGICAL, a, ac, ad, MPI_LOGICAL, l_comm, ierr)
  end subroutine parallel_alltoall_lv

  ! Broadcast
  subroutine parallel_bcast_d(a, root, comm)
    real(kind=dp_t), intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node; 
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, l_root, l_comm, ierr)
  end subroutine parallel_bcast_d
  subroutine parallel_bcast_r(a, root, comm)
    real(kind=sp_t), intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, 1, MPI_REAL, l_root, l_comm, ierr)
  end subroutine parallel_bcast_r
  subroutine parallel_bcast_i(a, root, comm)
    integer, intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, 1, MPI_INTEGER, l_root, l_comm, ierr)
  end subroutine parallel_bcast_i
  subroutine parallel_bcast_l(a, root, comm)
    logical, intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, 1, MPI_LOGICAL, l_root, l_comm, ierr)
  end subroutine parallel_bcast_l
  subroutine parallel_bcast_c(a, root, comm)
    complex(kind=sp_t), intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, 1, MPI_COMPLEX, l_root, l_comm, ierr)
  end subroutine parallel_bcast_c
  subroutine parallel_bcast_z(a, root, comm)
    complex(kind=dp_t), intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, 1, MPI_DOUBLE_COMPLEX, l_root, l_comm, ierr)
  end subroutine parallel_bcast_z
  ! vector versions
  subroutine parallel_bcast_dv(a, root, comm)
    real(kind=dp_t), intent(in) :: a(:)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node; 
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, size(a), MPI_DOUBLE_PRECISION, l_root, l_comm, ierr)
  end subroutine parallel_bcast_dv
  subroutine parallel_bcast_rv(a, root, comm)
    real(kind=sp_t), intent(in) :: a(:)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, size(a), MPI_REAL, l_root, l_comm, ierr)
  end subroutine parallel_bcast_rv
  subroutine parallel_bcast_iv(a, root, comm)
    integer, intent(in) :: a(:)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, size(a), MPI_INTEGER, l_root, l_comm, ierr)
  end subroutine parallel_bcast_iv
  subroutine parallel_bcast_lv(a, root, comm)
    logical, intent(in) :: a(:)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, size(a), MPI_LOGICAL, l_root, l_comm, ierr)
  end subroutine parallel_bcast_lv
  subroutine parallel_bcast_cv(a, root, comm)
    complex(kind=sp_t), intent(in) :: a(:)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, size(a), MPI_COMPLEX, l_root, l_comm, ierr)
  end subroutine parallel_bcast_cv
  subroutine parallel_bcast_zv(a, root, comm)
    complex(kind=dp_t), intent(in) :: a(:)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, size(a), MPI_DOUBLE_COMPLEX, l_root, l_comm, ierr)
  end subroutine parallel_bcast_zv


  ! Barrier:
  subroutine parallel_barrier(comm)
    integer, intent(in), optional :: comm
    integer :: ierr, l_comm
    external MPI_Barrier
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Barrier(l_comm, ierr)
  end subroutine parallel_barrier


  !
  ! Gather fixed size blocks to specified processor
  !
  subroutine parallel_gather_d(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    real(kind=dp_t), intent(in) :: snd(*)
    real(kind=dp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Gather
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gather(snd, n, MPI_DOUBLE_PRECISION, &
         rcv, n, MPI_DOUBLE_PRECISION, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_d
  subroutine parallel_gather_r(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    real(kind=sp_t), intent(in) :: snd(*)
    real(kind=sp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Gather
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gather(snd, n, MPI_REAL, &
         rcv, n, MPI_REAL, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_r
  subroutine parallel_gather_i(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    integer, intent(in) :: snd(*)
    integer, intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Gather
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gather(snd, n, MPI_INTEGER, &
         rcv, n, MPI_INTEGER, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_i
  subroutine parallel_gather_l(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    logical, intent(in) :: snd(*)
    logical, intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Gather
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gather(snd, n, MPI_LOGICAL, &
         rcv, n, MPI_LOGICAL, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_l
  subroutine parallel_gather_c(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    complex(kind=sp_t), intent(in) :: snd(*)
    complex(kind=sp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Gather
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gather(snd, n, MPI_COMPLEX, &
         rcv, n, MPI_COMPLEX, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_c
  subroutine parallel_gather_z(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    complex(kind=dp_t), intent(in) :: snd(*)
    complex(kind=dp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Gather
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gather(snd, n, MPI_DOUBLE_COMPLEX, &
         rcv, n, MPI_DOUBLE_COMPLEX, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_z
  !
  ! Gather variable sized blocks to specified processor
  !
  subroutine parallel_gather_dv(snd, n, rcv, rcvc, rcvd, root, comm)
    real(kind=dp_t), intent(in) :: snd(*)
    integer, intent(in) :: n
    real(kind=dp_t), intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_root, l_comm
    external MPI_Gatherv
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    call MPI_Gatherv(snd, n, MPI_DOUBLE_PRECISION, &
         rcv, rcvc, rcvd, MPI_DOUBLE_PRECISION, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_dv
  subroutine parallel_gather_rv(snd, n, rcv, rcvc, rcvd, root, comm)
    real(kind=sp_t), intent(in) :: snd(*)
    integer, intent(in) :: n
    real(kind=sp_t), intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_root, l_comm
    external MPI_Gatherv
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    call MPI_Gatherv(snd, n, MPI_REAL, &
         rcv, rcvc, rcvd, MPI_REAL, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_rv
  subroutine parallel_gather_iv(snd, n, rcv, rcvc, rcvd, root, comm)
    integer, intent(in) :: snd(*)
    integer, intent(in) :: n
    integer, intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_root, l_comm
    external MPI_Gatherv
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gatherv(snd, n, MPI_INTEGER, &
         rcv, rcvc, rcvd, MPI_INTEGER, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_iv
  subroutine parallel_gather_lv(snd, n, rcv, rcvc, rcvd, root, comm)
    logical, intent(in) :: snd(*)
    integer, intent(in) :: n
    logical, intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_root, l_comm
    external MPI_Gatherv
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gatherv(snd, n, MPI_LOGICAL, &
         rcv, rcvc, rcvd, MPI_LOGICAL, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_lv
  subroutine parallel_gather_cv(snd, n, rcv, rcvc, rcvd, root, comm)
    complex(sp_t), intent(in) :: snd(*)
    integer, intent(in) :: n
    complex(sp_t), intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_root, l_comm
    external MPI_Gatherv
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gatherv(snd, n, MPI_COMPLEX, &
         rcv, rcvc, rcvd, MPI_COMPLEX, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_cv
  subroutine parallel_gather_zv(snd, n, rcv, rcvc, rcvd, root, comm)
    complex(dp_t), intent(in) :: snd(*)
    integer, intent(in) :: n
    complex(dp_t), intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_root, l_comm
    external MPI_Gatherv
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gatherv(snd, n, MPI_DOUBLE_COMPLEX, &
         rcv, rcvc, rcvd, MPI_DOUBLE_COMPLEX, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_zv


  ! Allgather:
  subroutine parallel_allgather_dv(snd, rcv, n, comm)
    integer, intent(in) :: n
    real(kind=dp_t), intent(in) :: snd(*)
    real(kind=dp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: ierr, l_comm
    external MPI_Allgather
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Allgather(snd, n, MPI_DOUBLE_PRECISION, rcv, n, MPI_DOUBLE_PRECISION, l_comm, ierr)
  end subroutine parallel_allgather_dv
  subroutine parallel_allgather_rv(snd, rcv, n, comm)
    integer, intent(in) :: n
    real(kind=sp_t), intent(in) :: snd(*)
    real(kind=sp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: ierr, l_comm
    external MPI_Allgather
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Allgather(snd, n, MPI_REAL, rcv, n, MPI_REAL, l_comm, ierr)
  end subroutine parallel_allgather_rv
  subroutine parallel_allgather_iv(snd, rcv, n, comm)
    integer, intent(in) :: n
    integer, intent(in) :: snd(*)
    integer, intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: ierr, l_comm
    external MPI_Allgather
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Allgather(snd, n, MPI_INTEGER, rcv, n, MPI_INTEGER, l_comm, ierr)
  end subroutine parallel_allgather_iv
  subroutine parallel_allgather_lv(snd, rcv, n, comm)
    integer, intent(in) :: n
    logical, intent(in) :: snd(*)
    logical, intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: ierr, l_comm
    external MPI_Allgather
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Allgather(snd, n, MPI_LOGICAL, rcv, n, MPI_LOGICAL, l_comm, ierr)
  end subroutine parallel_allgather_lv
  subroutine parallel_allgather_cv(snd, rcv, n, comm)
    integer, intent(in) :: n
    complex(kind=sp_t), intent(in) :: snd(*)
    complex(kind=sp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: ierr, l_comm
    external MPI_Allgather
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Allgather(snd, n, MPI_COMPLEX, rcv, n, MPI_COMPLEX, l_comm, ierr)
  end subroutine parallel_allgather_cv
  subroutine parallel_allgather_zv(snd, rcv, n, comm)
    integer, intent(in) :: n
    complex(kind=dp_t), intent(in) :: snd(*)
    complex(kind=dp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: ierr, l_comm
    external MPI_Allgather
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Allgather(snd, n, MPI_DOUBLE_COMPLEX, rcv, n, MPI_DOUBLE_COMPLEX, l_comm, ierr)
  end subroutine parallel_allgather_zv

  ! Scatter:
  subroutine parallel_scatter_dv(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    real(kind=dp_t), intent(in) :: snd(*)
    real(kind=dp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Scatter
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Scatter(snd, n, MPI_DOUBLE_PRECISION, &
         rcv, n, MPI_DOUBLE_PRECISION, &
         l_root, l_comm, ierr)
  end subroutine parallel_scatter_dv
  subroutine parallel_scatter_rv(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    real(kind=sp_t), intent(in) :: snd(*)
    real(kind=sp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Scatter
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Scatter(snd, n, MPI_REAL, &
         rcv, n, MPI_REAL, &
         l_root, l_comm, ierr)
  end subroutine parallel_scatter_rv
  subroutine parallel_scatter_iv(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    integer, intent(in) :: snd(*)
    integer, intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Scatter
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Scatter(snd, n, MPI_INTEGER, &
         rcv, n, MPI_INTEGER, &
         l_root, l_comm, ierr)
  end subroutine parallel_scatter_iv
  subroutine parallel_scatter_lv(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    logical, intent(in) :: snd(*)
    logical, intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Scatter
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Scatter(snd, n, MPI_LOGICAL, &
         rcv, n, MPI_LOGICAL, &
         l_root, l_comm, ierr)
  end subroutine parallel_scatter_lv
  subroutine parallel_scatter_cv(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    complex(kind=sp_t), intent(in) :: snd(*)
    complex(kind=sp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Scatter
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Scatter(snd, n, MPI_COMPLEX, &
         rcv, n, MPI_COMPLEX, &
         l_root, l_comm, ierr)
  end subroutine parallel_scatter_cv
  subroutine parallel_scatter_zv(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    complex(kind=dp_t), intent(in) :: snd(*)
    complex(kind=dp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Scatter
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Scatter(snd, n, MPI_DOUBLE_COMPLEX, &
         rcv, n, MPI_DOUBLE_COMPLEX, &
         l_root, l_comm, ierr)
  end subroutine parallel_scatter_zv

end module parallel
