
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

end module memprof_module
