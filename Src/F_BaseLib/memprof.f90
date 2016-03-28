
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
    boxarray_bytes = (ms%num_alloc - ms%num_dealloc) * sizeof_box;
  end function boxarray_bytes

  function boxarray_bytes_hwm() bind(c, name='memprof_boxarray_bytes_hwm')
    use boxarray_module
    integer(kind=c_long) :: boxarray_bytes_hwm
    ! 7 integers/box * 4 bytes/integer = 28 bytes/box
    boxarray_bytes_hwm = boxarray_get_high_water_mark() * 28
  end function boxarray_bytes_hwm

  function tilearray_bytes() bind(c, name='memprof_tilearray_bytes')
    use layout_module, only : tilearray_total_bytes
    integer(kind=c_long) :: tilearray_bytes
    tilearray_bytes = 0_c_long
    !$omp parallel reduction(+:tilearray_bytes)
    tilearray_bytes = tilearray_bytes + tilearray_total_bytes
    !$omp end parallel
  end function tilearray_bytes

  function tilearray_bytes_hwm() bind(c, name='memprof_tilearray_bytes_hwm')
    use layout_module, only : tilearray_total_bytes_hwm
    integer(kind=c_long) :: tilearray_bytes_hwm
    tilearray_bytes_hwm = 0_c_long
    !$omp parallel reduction(+:tilearray_bytes_hwm)
    tilearray_bytes_hwm = tilearray_total_bytes_hwm
    !$omp end parallel
  end function tilearray_bytes_hwm

  function boxhash_bytes() bind(c, name='memprof_boxhash_bytes')
    use layout_module, only : boxhash_total_bytes
    integer(kind=c_long) :: boxhash_bytes
    boxhash_bytes = boxhash_total_bytes
  end function boxhash_bytes

  function boxhash_bytes_hwm() bind(c, name='memprof_boxhash_bytes_hwm')
    use layout_module, only : boxhash_total_bytes_hwm
    integer(kind=c_long) :: boxhash_bytes_hwm
    boxhash_bytes_hwm = boxhash_total_bytes_hwm
  end function boxhash_bytes_hwm

  function boxassoc_bytes() bind(c, name='memprof_boxassoc_bytes')
    use layout_module, only : boxassoc_total_bytes
    integer(kind=c_long) :: boxassoc_bytes
    boxassoc_bytes = boxassoc_total_bytes
  end function boxassoc_bytes

  function boxassoc_bytes_hwm() bind(c, name='memprof_boxassoc_bytes_hwm')
    use layout_module, only : boxassoc_total_bytes_hwm
    integer(kind=c_long) :: boxassoc_bytes_hwm
    boxassoc_bytes_hwm = boxassoc_total_bytes_hwm
  end function boxassoc_bytes_hwm

  function fgassoc_bytes() bind(c, name='memprof_fgassoc_bytes')
    use layout_module, only : fgassoc_total_bytes
    integer(kind=c_long) :: fgassoc_bytes
    fgassoc_bytes = fgassoc_total_bytes
  end function fgassoc_bytes

  function fgassoc_bytes_hwm() bind(c, name='memprof_fgassoc_bytes_hwm')
    use layout_module, only : fgassoc_total_bytes_hwm
    integer(kind=c_long) :: fgassoc_bytes_hwm
    fgassoc_bytes_hwm = fgassoc_total_bytes_hwm
  end function fgassoc_bytes_hwm

  function syncassoc_bytes() bind(c, name='memprof_syncassoc_bytes')
    use layout_module, only : syncassoc_total_bytes
    integer(kind=c_long) :: syncassoc_bytes
    syncassoc_bytes = syncassoc_total_bytes
  end function syncassoc_bytes

  function syncassoc_bytes_hwm() bind(c, name='memprof_syncassoc_bytes_hwm')
    use layout_module, only : syncassoc_total_bytes_hwm
    integer(kind=c_long) :: syncassoc_bytes_hwm
    syncassoc_bytes_hwm = syncassoc_total_bytes_hwm
  end function syncassoc_bytes_hwm

  function copyassoc_bytes() bind(c, name='memprof_copyassoc_bytes')
    use layout_module, only : copyassoc_total_bytes
    integer(kind=c_long) :: copyassoc_bytes
    copyassoc_bytes = copyassoc_total_bytes
  end function copyassoc_bytes

  function copyassoc_bytes_hwm() bind(c, name='memprof_copyassoc_bytes_hwm')
    use layout_module, only : copyassoc_total_bytes_hwm
    integer(kind=c_long) :: copyassoc_bytes_hwm
    copyassoc_bytes_hwm = copyassoc_total_bytes_hwm
  end function copyassoc_bytes_hwm

  function fluxassoc_bytes() bind(c, name='memprof_fluxassoc_bytes')
    use layout_module, only : fluxassoc_total_bytes
    integer(kind=c_long) :: fluxassoc_bytes
    fluxassoc_bytes = fluxassoc_total_bytes
  end function fluxassoc_bytes

  function fluxassoc_bytes_hwm() bind(c, name='memprof_fluxassoc_bytes_hwm')
    use layout_module, only : fluxassoc_total_bytes_hwm
    integer(kind=c_long) :: fluxassoc_bytes_hwm
    fluxassoc_bytes_hwm = fluxassoc_total_bytes_hwm
  end function fluxassoc_bytes_hwm

end module memprof_module
