module omp_module

  implicit none

contains

  subroutine omp_set_num_threads(np)
    integer np
    if ( np /= 1 ) &
         print *, 'OMP_SET_NUM_THREADS(',np, '):omp_stub:WARN: no threads'
  end subroutine omp_set_num_threads
  integer function omp_get_num_threads()
    omp_get_num_threads = 1
  end function omp_get_num_threads
  integer function omp_get_max_threads()
    omp_get_max_threads = 1
  end function omp_get_max_threads
  integer function omp_get_thread_num()
    omp_get_thread_num = 0
  end function omp_get_thread_num
  integer function omp_get_num_procs()
    omp_get_num_procs = 1
  end function omp_get_num_procs
  subroutine omp_set_dynamic(flag)
    logical flag
    if ( .not. flag ) &
         print *, 'OMP_SET_DYNAMIC(',flag,':omp_stub:WARN: no threads'
  end subroutine omp_set_dynamic
  logical function omp_get_dynamic()
    omp_get_dynamic = .false.
  end function omp_get_dynamic
  logical function omp_in_parallel()
    omp_in_parallel = .false.
  end function omp_in_parallel
  subroutine omp_set_nested(flag)
    logical flag
    if ( .not. flag ) &
         print *, 'OMP_SET_NESTED(',flag,'):omp_stub:WARN: no threads'
  end subroutine omp_set_nested
  logical function omp_get_nested()
    omp_get_nested = .false.
  end function omp_get_nested
  subroutine omp_init_lock(lock)
    integer lock
    lock = -1
  end subroutine omp_init_lock
  subroutine omp_destroy_lock(lock)
    integer lock
    lock = 0
  end subroutine omp_destroy_lock
  subroutine omp_set_lock(lock)
    integer lock
    if(lock .eq. 0) then
       stop 'OMP_SET_LOCK:ERROR: lock not initialized'
    else if(lock .eq. 1) then
       stop 'OMP_SET_LOCK:ERROR: deadlock in using lock variable'
    else
       lock = 1
    end if
  end subroutine omp_set_lock
  subroutine omp_unset_lock(lock)
    integer lock
    if(lock .eq. 0) then
       stop 'OMP_UNSET_LOCK:ERROR: lock not initialized'
    else if(lock .eq. 1) then
       lock = -1
    else
       stop 'OMP_UNSET_LOCK:ERROR: lock not set'
    end if
  end subroutine omp_unset_lock
  logical function omp_test_lock(lock)
    integer lock
    if (lock .eq. -1) then
       lock = 1
       omp_test_lock = .true.
    else if(lock .eq. 1) then
       omp_test_lock = .false.
    else
       stop 'OMP_TEST_LOCK:ERROR: lock not initialized'
    end if
  end function omp_test_lock

end module omp_module
