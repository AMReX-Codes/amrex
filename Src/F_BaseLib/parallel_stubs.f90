! Parallel wrappers
!! These wrappers are used so that a non-MPI version can coexist
!! with the MPI version.

module parallel

  use bl_types

  implicit none


  ! Some selected values based on MPICH/1.2.5 unix implementation

  integer, parameter :: MPI_COMM_WORLD = 0

  integer, parameter :: MPI_UNDEFINED   = -32766
  integer, parameter :: MPI_STATUS_SIZE = 4
  integer, parameter :: MPI_SOURCE      = 2
  integer, parameter :: MPI_TAG         = 3
  integer, parameter :: MPI_ERROR       = 4

  integer, parameter :: MPI_MAX       = 100
  integer, parameter :: MPI_MIN       = 101
  integer, parameter :: MPI_SUM       = 102
  integer, parameter :: MPI_PROD      = 103
  integer, parameter :: MPI_LAND      = 104
  integer, parameter :: MPI_BAND      = 105
  integer, parameter :: MPI_LOR       = 106
  integer, parameter :: MPI_BOR       = 107
  integer, parameter :: MPI_LXOR      = 108
  integer, parameter :: MPI_BXOR      = 109
  integer, parameter :: MPI_MINLOC    = 110
  integer, parameter :: MPI_MAXLOC    = 111

  integer, parameter :: MPI_ANY_SOURCE = -2
  integer, parameter :: MPI_ANY_TAG    = -1

  integer, parameter :: MPI_THREAD_SINGLE     = 0
  integer, parameter :: MPI_THREAD_FUNNELED   = 1
  integer, parameter :: MPI_THREAD_SERIALIZED = 2
  integer, parameter :: MPI_THREAD_MULTIPLE   = 3

  integer, parameter, private :: io_processor_node = 0
  integer, private :: m_nprocs = 1
  integer, private :: m_myproc = 0
  integer, private :: m_comm   = -1
  integer, private :: m_thread_support_level = 0

  ! interface communicator
  !    module procedure parallel_communicator
  ! end interface communicator
  ! interface nprocs
  !   module procedure parallel_nprocs
  ! end interface nprocs
  ! interface myproc
  !    module procedure parallel_myproc
  ! end interface myproc
  ! interface IOProcessor
  !    module procedure parallel_IOProcessor
  ! end interface IOProcessor
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
     module procedure parallel_send_l1
     module procedure parallel_send_l2
     module procedure parallel_send_l3
     module procedure parallel_send_l4
     module procedure parallel_send_l5
     module procedure parallel_send_l6
     module procedure parallel_send_l7
     module procedure parallel_send_i1
     module procedure parallel_send_i2
     module procedure parallel_send_i3
     module procedure parallel_send_i4
     module procedure parallel_send_i5
     module procedure parallel_send_i6
     module procedure parallel_send_i7
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
     module procedure parallel_recv_l1
     module procedure parallel_recv_l2
     module procedure parallel_recv_l3
     module procedure parallel_recv_l4
     module procedure parallel_recv_l5
     module procedure parallel_recv_l6
     module procedure parallel_recv_l7
     module procedure parallel_recv_i1
     module procedure parallel_recv_i2
     module procedure parallel_recv_i3
     module procedure parallel_recv_i4
     module procedure parallel_recv_i5
     module procedure parallel_recv_i6
     module procedure parallel_recv_i7
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
  end interface parallel_wait

  interface parallel_test
     module procedure parallel_test_one
     module procedure parallel_test_vec
     module procedure parallel_test_vec_vec
  end interface parallel_test

  interface parallel_bcast
     module procedure parallel_bcast_d
     module procedure parallel_bcast_r
     module procedure parallel_bcast_i
     module procedure parallel_bcast_l
     module procedure parallel_bcast_dv
     module procedure parallel_bcast_rv
     module procedure parallel_bcast_iv
     module procedure parallel_bcast_lv
  end interface parallel_bcast

  interface parallel_scatter
     module procedure parallel_scatter_dv
     module procedure parallel_scatter_rv
     module procedure parallel_scatter_iv
     module procedure parallel_scatter_lv
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

  interface
     subroutine sys_abort()
     end subroutine sys_abort
  end interface

  logical, private :: g_init = .False.

contains

  function parallel_initialized() result(r)
    logical :: r
    r = g_init
  end function parallel_initialized

  subroutine parallel_initialize(comm, thread_support_level)
    integer, intent(in), optional :: comm, thread_support_level
    g_init = .True.
  end subroutine parallel_initialize

  subroutine parallel_finalize(do_finalize_MPI)
    logical, intent(in), optional :: do_finalize_MPI
    g_init = .False.
  end subroutine parallel_finalize

  function parallel_q() result(r)
    logical :: r
    r = .FALSE.
  end function parallel_q

  subroutine parallel_abort(str)
    character(len=*), optional :: str
    if ( present(str) ) then
       print*, 'parallel_abort(): ', str
    else
       print*, 'parallel_abort() !!!'
    end if
    call sys_abort()
  end subroutine parallel_abort

  subroutine parallel_set_comm(comm)
    integer, intent(in) :: comm
  end subroutine parallel_set_comm

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
  pure function parallel_IOProcessor(comm) result(r)
    logical :: r
    integer, intent(in), optional :: comm
    r = (parallel_myproc() == io_processor_node)
  end function parallel_IOProcessor
  pure function parallel_IOProcessorNode() result(r)
    integer :: r
    r = io_processor_node
  end function parallel_IOProcessorNode

  pure function parallel_null_communicator() result(r)
    integer :: r
    !
    ! Just gotta return something that isn't m_comm.
    !
    r = (m_comm - 1)
  end function parallel_null_communicator

  function parallel_create_communicator(procs) result(comm)
    integer              :: comm
    integer, intent(in)  :: procs(:)
    !
    ! We're running in serial.  Just return m_comm
    !
    comm = m_comm
  end function parallel_create_communicator

  subroutine parallel_free_communicator(comm)
    integer, intent(in) :: comm
  end subroutine parallel_free_communicator

  pure function parallel_thread_support_level() result(r)
    integer :: r
    r = m_thread_support_level
  end function parallel_thread_support_level
  function parallel_wtime() result(r)
    real(kind=dp_t) :: r
    interface
       subroutine wall_second(s)
         use bl_types
         real(kind=dp_t), intent(out) :: s
       end subroutine wall_second
    end interface
    call wall_second(r)
  end function parallel_wtime
  function parallel_wtick() result(r)
    real(kind=dp_t) :: r
    interface
       subroutine wall_second_tick(s)
         use bl_types
         real(kind=dp_t), intent(out) :: s
       end subroutine wall_second_tick
    end interface
    call wall_second_tick(r)
  end function parallel_wtick

  subroutine parallel_reduce_d(r, a, op, proc, comm)
    real(kind=dp_t), intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc, comm
    real(kind=dp_t), intent(out) :: r
    r = a
  end subroutine parallel_reduce_d
  subroutine parallel_reduce_r(r, a, op, proc, comm)
    real(kind=sp_t), intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc, comm
    real(kind=sp_t), intent(out) :: r
    r = a
  end subroutine parallel_reduce_r
  subroutine parallel_reduce_i(r, a, op, proc, comm)
    integer, intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc, comm
    integer, intent(out) :: r
    r = a
  end subroutine parallel_reduce_i
  subroutine parallel_reduce_l(r, a, op, proc, comm)
    logical, intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc, comm
    logical, intent(out) :: r
    r = a
  end subroutine parallel_reduce_l
  subroutine parallel_reduce_dv(r, a, op, proc, comm)
    real(kind=dp_t), intent(in) :: a(:)
    integer, intent(in) :: op
    integer, intent(in), optional :: proc, comm
    real(kind=dp_t), intent(out) :: r(:)
    r = a
  end subroutine parallel_reduce_dv
  subroutine parallel_reduce_rv(r, a, op, proc, comm)
    real(kind=sp_t), intent(in) :: a(:)
    integer, intent(in) :: op
    integer, intent(in), optional :: proc, comm
    real(kind=sp_t), intent(out) :: r(:)
    r = a
  end subroutine parallel_reduce_rv
  subroutine parallel_reduce_iv(r, a, op, proc, comm)
    integer, intent(in) :: a(:)
    integer, intent(in) :: op
    integer, intent(in), optional :: proc, comm
    integer, intent(out) :: r(:)
    r = a
  end subroutine parallel_reduce_iv
  subroutine parallel_reduce_lv(r, a, op, proc, comm)
    logical, intent(in) :: a(:)
    integer, intent(in) :: op
    integer, intent(in), optional :: proc, comm
    logical, intent(out) :: r(:)
    r = a
  end subroutine parallel_reduce_lv

  function parallel_isend_dv(a, n, proc, tag, comm) result(r)
    integer :: r
    real(kind=dp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_DV')
    r = -1
  end function parallel_isend_dv
  function parallel_isend_rv(a, n, proc, tag, comm) result(r)
    integer :: r
    real(kind=sp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_RV')
    r = -1
  end function parallel_isend_rv
  function parallel_isend_iv(a, n, proc, tag, comm) result(r)
    integer :: r
    integer, intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_RV')
    r = -1
  end function parallel_isend_iv
  function parallel_isend_lv(a, n, proc, tag, comm) result(r)
    integer :: r
    logical, intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_RV')
    r = -1
  end function parallel_isend_lv
  function parallel_isend_zv(a, n, proc, tag, comm) result(r)
    integer :: r
    complex(kind=dp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_DV')
    r = -1
  end function parallel_isend_zv



  function parallel_irecv_dv(a, n, proc, tag, comm) result(r)
    integer :: r
    real(kind=dp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_RV')
    a(1:n) = HUGE(a)
    r = -1
  end function parallel_irecv_dv
  function parallel_irecv_rv(a, n, proc, tag, comm) result(r)
    integer :: r
    real(kind=sp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_RV')
    a(1:n) = HUGE(a)
    r = -1
  end function parallel_irecv_rv
  function parallel_irecv_iv(a, n, proc, tag, comm) result(r)
    integer :: r
    integer, intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_RV')
    a(1:n) = HUGE(a)
    r = -1
  end function parallel_irecv_iv
  function parallel_irecv_lv(a, n, proc, tag, comm) result(r)
    integer :: r
    logical, intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_RV')
    a(1:n) = .FALSE.
    r = -1
  end function parallel_irecv_lv
  function parallel_irecv_zv(a, n, proc, tag, comm) result(r)
    integer :: r
    complex(dp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_RV')
    a(1:n) = Huge(1.0_dp_t)
    r = -1
  end function parallel_irecv_zv

  subroutine parallel_wait_one(req, status)
    integer, intent(in)  :: req
    integer, intent(out), optional :: status(MPI_STATUS_SIZE)
    integer :: lstatus(MPI_STATUS_SIZE)
    lstatus = 0
    call parallel_abort('PARALLEL_WAIT')
    if ( present(status) ) status = lstatus
  end subroutine parallel_wait_one

  subroutine parallel_wait_vec(req, index, status)
    integer, intent(in)  :: req(:)
    integer, intent(out) :: index
    integer, intent(out), optional :: status(MPI_STATUS_SIZE)
    integer :: lstatus(MPI_STATUS_SIZE)
    lstatus = 0
    if ( size(req) > 0 ) call parallel_abort('PARALLEL_WAIT')
    if ( present(status) ) status = lstatus
    index = MPI_UNDEFINED
  end subroutine parallel_wait_vec

  subroutine parallel_wait_vec_vec(req, status)
    integer, intent(in)  :: req(:)
    integer, intent(out), optional :: status(:,:)
    integer :: lstatus(MPI_STATUS_SIZE,size(req))
    lstatus = 0
    if ( size(req) > 0 ) call parallel_abort('PARALLEL_WAIT')
    if ( present(status) ) status = lstatus
  end subroutine parallel_wait_vec_vec



  function parallel_test_one(req, status) result(r)
    integer, intent(in)  :: req
    logical :: r
    integer, intent(out), optional :: status(MPI_STATUS_SIZE)
    integer :: lstatus(MPI_STATUS_SIZE)
    lstatus = 0
    call parallel_abort('PARALLEL_TEST')
    if ( present(status) ) status = lstatus
    r = .FALSE.
  end function parallel_test_one

  function parallel_test_vec(req, index, status) result(r)
    integer, intent(in)  :: req(:)
    integer, intent(out) :: Index
    logical :: r
    integer, intent(out), optional :: status(MPI_STATUS_SIZE)
    integer :: lstatus(MPI_STATUS_SIZE)
    lstatus = 0
    if ( size(req) > 0 ) call parallel_abort('PARALLEL_TEST')
    if ( present(status) ) status = lstatus
    index = MPI_UNDEFINED
    r = .FALSE.
  end function parallel_test_vec

  function parallel_test_vec_vec(req, status) result(r)
    integer, intent(in)  :: req(:)
    logical :: r
    integer, intent(out), optional :: status(:,:)
    integer :: lstatus(MPI_STATUS_SIZE,size(req))
    lstatus = 0
    if ( size(req) > 0 ) call parallel_abort('PARALLEL_TEST')
    if ( present(status) ) status = lstatus
    r = .FALSE.
  end function parallel_test_vec_vec



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
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_SEND')
  end subroutine parallel_send_dv
  subroutine parallel_send_rv(a, n, proc, tag, comm)
    real(kind=sp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_SEND')
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
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_SEND')
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
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_SEND')
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
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_SEND')
  end subroutine parallel_send_zv

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
    integer l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    lstatus = 0
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_RECV')
    a(1:n) = HUGE(a)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_dv
  subroutine parallel_recv_rv(a, n, proc, tag, comm, status)
    real(kind=sp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    lstatus = 0
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_RECV')
    a(1:n) = HUGE(a)
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
    integer l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    lstatus = 0
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_RECV')
    a(1:n) = HUGE(a)
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
    integer l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    lstatus = 0
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_RECV')
    a(1:n) = .FALSE.
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
    integer l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    lstatus = 0
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_RECV')
    a(1:n) = HUGE(1.0_dp_t)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_zv

  ! AlltoAll
  subroutine parallel_alltoall_d(a, b, n, comm)
    real(kind=dp_t), intent(in) :: b(*)
    real(kind=dp_t), intent(inout) :: a(*)
    integer, intent(in), optional :: comm
    integer l_comm
    integer, intent(in) :: n
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    a(1:n) = b(1:n)
  end subroutine parallel_alltoall_d
  ! AlltoAll
  subroutine parallel_alltoall_dv(a, ac, ad, b, bc, bd, comm)
    real(kind=dp_t), intent(in) :: b(*)
    real(kind=dp_t), intent(inout) :: a(*)
    integer, intent(in) :: ac(*), ad(*), bc(*), bd(*)
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    a(ad(1):ad(1)+ac(1)-1) = b(bd(1):bd(1)+bc(1)-1)
  end subroutine parallel_alltoall_dv

  subroutine parallel_alltoall_i(a, b, n, comm)
    integer, intent(in) :: b(*)
    integer, intent(inout) :: a(*)
    integer, intent(in), optional :: comm
    integer l_comm
    integer, intent(in) :: n
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    a(1:n) = b(1:n)
  end subroutine parallel_alltoall_i
  ! AlltoAll
  subroutine parallel_alltoall_iv(a, ac, ad, b, bc, bd, comm)
    integer, intent(in) :: b(*)
    integer, intent(inout) :: a(*)
    integer, intent(in) :: ac(*), ad(*), bc(*), bd(*)
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    a(ad(1):ad(1)+ac(1)-1) = b(bd(1):bd(1)+bc(1)-1)
  end subroutine parallel_alltoall_iv

  subroutine parallel_alltoall_l(a, b, n, comm)
    logical, intent(in) :: b(*)
    logical, intent(inout) :: a(*)
    integer, intent(in), optional :: comm
    integer l_comm
    integer, intent(in) :: n
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    a(1:n) = b(1:n)
  end subroutine parallel_alltoall_l
  ! AlltoAll
  subroutine parallel_alltoall_lv(a, ac, ad, b, bc, bd, comm)
    logical, intent(in) :: b(*)
    logical, intent(inout) :: a(*)
    integer, intent(in) :: ac(*), ad(*), bc(*), bd(*)
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    a(ad(1):ad(1)+ac(1)-1) = b(bd(1):bd(1)+bc(1)-1)
  end subroutine parallel_alltoall_lv

  subroutine parallel_bcast_d(a, root, comm)
    real(kind=dp_t), intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_comm, l_root
    l_root  = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
  end subroutine parallel_bcast_d
  subroutine parallel_bcast_r(a, root, comm)
    real(kind=sp_t), intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_comm, l_root
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
  end subroutine parallel_bcast_r
  subroutine parallel_bcast_i(a, root, comm)
    integer, intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_comm, l_root
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
  end subroutine parallel_bcast_i
  subroutine parallel_bcast_l(a, root, comm)
    logical, intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_comm, l_root
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
  end subroutine parallel_bcast_l
  ! vector versions
  subroutine parallel_bcast_dv(a, root, comm)
    real(kind=dp_t), intent(in) :: a(:)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_comm, l_root
    l_root  = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
  end subroutine parallel_bcast_dv
  subroutine parallel_bcast_rv(a, root, comm)
    real(kind=sp_t), intent(in) :: a(:)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_comm, l_root
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
  end subroutine parallel_bcast_rv
  subroutine parallel_bcast_iv(a, root, comm)
    integer, intent(in) :: a(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_comm, l_root
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
  end subroutine parallel_bcast_iv
  subroutine parallel_bcast_lv(a, root, comm)
    logical, intent(in) :: a(:)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_comm, l_root
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
  end subroutine parallel_bcast_lv



  subroutine parallel_barrier(comm)
    integer, intent(in), optional :: comm
    integer :: l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
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
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_d
  subroutine parallel_gather_r(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    real(kind=sp_t), intent(in) :: snd(*)
    real(kind=sp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_r
  subroutine parallel_gather_i(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    integer, intent(in) :: snd(*)
    integer, intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_i
  subroutine parallel_gather_l(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    logical, intent(in) :: snd(*)
    logical, intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_l
  subroutine parallel_gather_c(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    complex(sp_t), intent(in) :: snd(*)
    complex(sp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_c
  subroutine parallel_gather_z(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    complex(dp_t), intent(in) :: snd(*)
    complex(dp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
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
    integer l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_dv
  subroutine parallel_gather_rv(snd, n, rcv, rcvc, rcvd, root, comm)
    real(kind=sp_t), intent(in) :: snd(*)
    integer, intent(in) :: n
    real(kind=sp_t), intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_rv
  subroutine parallel_gather_iv(snd, n, rcv, rcvc, rcvd, root, comm)
    integer, intent(in) :: snd(*)
    integer, intent(in) :: n
    integer, intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_iv
  subroutine parallel_gather_lv(snd, n, rcv, rcvc, rcvd, root, comm)
    logical, intent(in) :: snd(*)
    integer, intent(in) :: n
    logical, intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_lv
  subroutine parallel_gather_cv(snd, n, rcv, rcvc, rcvd, root, comm)
    complex(sp_t), intent(in) :: snd(*)
    integer, intent(in) :: n
    complex(sp_t), intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_cv
  subroutine parallel_gather_zv(snd, n, rcv, rcvc, rcvd, root, comm)
    complex(dp_t), intent(in) :: snd(*)
    integer, intent(in) :: n
    complex(dp_t), intent(inout) :: rcv(*)
    integer, intent(in) :: rcvc(*), rcvd(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_zv


  subroutine parallel_allgather_dv(snd, rcv, n, comm)
    integer, intent(in) :: n
    real(kind=dp_t), intent(in) :: snd(*)
    real(kind=dp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_allgather_dv
  subroutine parallel_allgather_rv(snd, rcv, n, comm)
    integer, intent(in) :: n
    real(kind=sp_t), intent(in) :: snd(*)
    real(kind=sp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_allgather_rv
  subroutine parallel_allgather_iv(snd, rcv, n, comm)
    integer, intent(in) :: n
    integer, intent(in) :: snd(*)
    integer, intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_allgather_iv
  subroutine parallel_allgather_lv(snd, rcv, n, comm)
    integer, intent(in) :: n
    logical, intent(in) :: snd(*)
    logical, intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_allgather_lv
  subroutine parallel_allgather_zv(snd, rcv, n, comm)
    integer, intent(in) :: n
    complex(dp_t), intent(in) :: snd(*)
    complex(dp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_allgather_zv


  subroutine parallel_scatter_dv(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    real(kind=dp_t), intent(in) :: snd(*)
    real(kind=dp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_scatter_dv
  subroutine parallel_scatter_rv(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    real(kind=sp_t), intent(in) :: snd(*)
    real(kind=sp_t), intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_scatter_rv
  subroutine parallel_scatter_iv(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    integer, intent(in) :: snd(*)
    integer, intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_scatter_iv
  subroutine parallel_scatter_lv(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    logical, intent(in) :: snd(*)
    logical, intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_scatter_lv

end module parallel
