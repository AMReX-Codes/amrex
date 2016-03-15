
module memprof_module

  use iso_c_binding

  implicit none

contains

  function fab_numdoubles() bind(c, name='memprof_fab_numdoubles')
    use fab_module
    integer(kind=c_long) :: fab_numdoubles
    type(mem_stats) :: ms
    ms = fab_mem_stats()
    fab_numdoubles = ms%num_alloc - ms%num_dealloc
  end function fab_numdoubles

  function fab_numdoubles_hwm() bind(c, name='memprof_fab_numdoubles_hwm')
    use fab_module
    integer(kind=c_long) :: fab_numdoubles_hwm
    fab_numdoubles_hwm = fab_get_high_water_mark()
  end function fab_numdoubles_hwm

  function boxarray_bytes() bind(c, name='memprof_boxarray_bytes')
    use boxarray_module
    integer(kind=c_long) :: boxarray_bytes
    type(mem_stats) :: ms
    ms = boxarray_mem_stats()
    ! 7 integers/box * 4 bytes/integer = 28 bytes/box
    boxarray_bytes = (ms%num_alloc - ms%num_dealloc) * 28
  end function boxarray_bytes

  function boxarray_bytes_hwm() bind(c, name='memprof_boxarray_bytes_hwm')
    use boxarray_module
    integer(kind=c_long) :: boxarray_bytes_hwm
    ! 7 integers/box * 4 bytes/integer = 28 bytes/box
    boxarray_bytes_hwm = boxarray_get_high_water_mark() * 28
  end function boxarray_bytes_hwm

end module memprof_module
