module bl_mem_stat_module
  
  use bl_types
  use parallel

  implicit none

  type mem_stats
     integer(kind =ll_t) :: num_alloc = 0_ll_t
     integer(kind =ll_t) :: num_dealloc = 0_ll_t
     integer(kind =ll_t) :: cnt_alloc = 0_ll_t
     integer(kind =ll_t) :: cnt_dealloc = 0_ll_t
  end type mem_stats

  interface print
     module procedure mem_stats_print
  end interface

  interface mem_stats_alloc
     module procedure mem_stats_alloc_ll
     module procedure mem_stats_alloc_i
     module procedure mem_stats_alloc_c
  end interface

  interface mem_stats_dealloc
     module procedure mem_stats_dealloc_ll
     module procedure mem_stats_dealloc_i
     module procedure mem_stats_dealloc_c
  end interface

contains

  subroutine mem_stats_print(ms, str, unit, advance, total)
    use bl_IO_module
    use bl_string_module
    type(mem_stats), intent(in) :: ms
    character(len=*), intent(in), optional :: str
    character(len=*), intent(in), optional :: advance
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: total
    integer :: un
    character(len=3) :: adv
    logical :: ltotal
    ltotal = .true.; if ( present(total) ) ltotal = total
    if ( parallel_IOProcessor() ) then
       un = unit_stdout(unit)
       adv = unit_advance(advance)
       if ( present(str) ) then
          write(unit = un, fmt = '(A,": ")', advance = 'no') str
       end if
       if ( ltotal ) then
          write(unit=un, fmt='(4i15)', advance = 'no') &
               ms%num_alloc,                           &
               ms%num_alloc-ms%num_dealloc,            &
               ms%cnt_alloc,                           &
               ms%cnt_alloc-ms%cnt_dealloc
       else
          write(unit=un, fmt='(2i15)', advance = 'no') &
               ms%num_alloc-ms%num_dealloc,            &
               ms%cnt_alloc-ms%cnt_dealloc
       end if
       if ( eq_i(adv, "YES") ) write(unit=un, fmt='()')
    end if
  end subroutine mem_stats_print

  subroutine mem_stats_alloc_ll(ms, vol)
    type(mem_stats), intent(inout) :: ms
    integer(kind=ll_t) :: vol
    ms%num_alloc = ms%num_alloc + vol
    ms%cnt_alloc = ms%cnt_alloc + 1
  end subroutine mem_stats_alloc_ll

  subroutine mem_stats_alloc_i(ms, vol)
    type(mem_stats), intent(inout) :: ms
    integer :: vol
    ms%num_alloc = ms%num_alloc + int(vol,kind=ll_t)
    ms%cnt_alloc = ms%cnt_alloc + 1
  end subroutine mem_stats_alloc_i

  subroutine mem_stats_alloc_c(ms)
    type(mem_stats), intent(inout) :: ms
    ms%cnt_alloc = ms%cnt_alloc + 1
  end subroutine mem_stats_alloc_c

  subroutine mem_stats_dealloc_ll(ms, vol)
    type(mem_stats), intent(inout) :: ms
    integer(kind=ll_t) :: vol
    ms%num_dealloc = ms%num_dealloc + vol
    ms%cnt_dealloc = ms%cnt_dealloc + 1
  end subroutine mem_stats_dealloc_ll

  subroutine mem_stats_dealloc_i(ms, vol)
    type(mem_stats), intent(inout) :: ms
    integer :: vol
    ms%num_dealloc = ms%num_dealloc + int(vol,kind=ll_t)
    ms%cnt_dealloc = ms%cnt_dealloc + 1
  end subroutine mem_stats_dealloc_i

  subroutine mem_stats_dealloc_c(ms)
    type(mem_stats), intent(inout) :: ms
    ms%cnt_dealloc = ms%cnt_dealloc + 1
  end subroutine mem_stats_dealloc_c

end module bl_mem_stat_module
