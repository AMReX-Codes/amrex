module amrex_timer_c_module

  implicit none

  interface

     subroutine cpu_second(s) bind(c, name='cpu_second')
       use bl_types
       real(kind=dp_t) :: s
     end subroutine cpu_second

     subroutine wall_second(s) bind(c, name='wall_second')
       use bl_types
       real(kind=dp_t) :: s
     end subroutine wall_second

     subroutine cpu_second_tick(s) bind(c, name='cpu_second_tick')
       use bl_types
       real(kind=dp_t) :: s
     end subroutine cpu_second_tick

     subroutine wall_second_tick(s) bind(c, name='wall_second_tick')
       use bl_types
       real(kind=dp_t) :: s
     end subroutine wall_second_tick

     subroutine sys_abort() bind(c, name='sys_abort')
     end subroutine sys_abort

  end interface

end module amrex_timer_c_module
