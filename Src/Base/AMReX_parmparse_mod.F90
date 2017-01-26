module amrex_parmparse_module

  use iso_c_binding
  use amrex_string_module

  implicit none

  private

  public :: parmparse_build, parmparse_destroy

  type, public :: amrex_parmparse
     type(c_ptr) :: p = c_null_ptr
   contains
     generic :: get      => get_int, get_double, get_logical, get_string
     generic :: query    => query_int, query_double, query_logical, query_string
     procedure, private :: get_int
     procedure, private :: get_double
     procedure, private :: get_logical
     procedure, private :: get_string
     procedure, private :: query_int
     procedure, private :: query_double
     procedure, private :: query_logical
     procedure, private :: query_string
  end type amrex_parmparse

  ! interfaces to cpp functions

  interface
     subroutine amrex_new_parmparse (pp, name) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr) :: pp
       character(c_char), intent(in) :: name(*)
     end subroutine amrex_new_parmparse

     subroutine amrex_delete_parmparse (pp) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
     end subroutine amrex_delete_parmparse

     subroutine amrex_parmparse_get_int (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end subroutine amrex_parmparse_get_int

     subroutine amrex_parmparse_get_double (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       real(c_double) :: v
     end subroutine amrex_parmparse_get_double

     subroutine amrex_parmparse_get_bool (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end subroutine amrex_parmparse_get_bool

     subroutine amrex_parmparse_get_string (pp, name, v, len) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       character(c_char), intent(inout) :: v(*)
       integer :: len
     end subroutine amrex_parmparse_get_string

     subroutine amrex_parmparse_query_int (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end subroutine amrex_parmparse_query_int

     subroutine amrex_parmparse_query_double (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       real(c_double) :: v
     end subroutine amrex_parmparse_query_double

     subroutine amrex_parmparse_query_bool (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end subroutine amrex_parmparse_query_bool

     subroutine amrex_parmparse_query_string (pp, name, v, len) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       character(c_char), intent(inout) :: v(*)
       integer :: len
     end subroutine amrex_parmparse_query_string
  end interface

contains

  subroutine parmparse_build (pp, name)
    type(amrex_parmparse) :: pp
    character(*), intent(in) :: name
    call amrex_new_parmparse(pp%p, string_f_to_c(name))
  end subroutine parmparse_build

  subroutine parmparse_destroy (this)
    type(amrex_parmparse) :: this
    if (c_associated(this%p)) then
       call amrex_delete_parmparse(this%p)
       this%p = c_null_ptr
    end if
  end subroutine parmparse_destroy

  subroutine get_int (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(len=*), intent(in) :: name
    integer :: v
    call amrex_parmparse_get_int (this%p, string_f_to_c(name), v)
  end subroutine get_int

  subroutine get_double (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    double precision :: v
    call amrex_parmparse_get_double (this%p, string_f_to_c(name), v)
  end subroutine get_double

  subroutine get_logical (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    logical :: v
    integer(c_int) :: i
    call amrex_parmparse_get_bool (this%p, string_f_to_c(name), i)
    v = i.eq.1
  end subroutine get_logical

  subroutine get_string (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    character(*), intent(inout) :: v

    ! temporary string for passing back and forth to C -- include NULL
    character(c_char), dimension(len(v)+1) :: v_pass

    call amrex_parmparse_get_string (this%p, string_f_to_c(name), v_pass, len(v)+1)

    ! convert to Fortran string
    v = string_c_to_f(v_pass)
  end subroutine get_string

  subroutine query_int (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(len=*), intent(in) :: name
    integer :: v
    call amrex_parmparse_query_int (this%p, string_f_to_c(name), v)
  end subroutine query_int

  subroutine query_double (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    double precision :: v
    call amrex_parmparse_query_double (this%p, string_f_to_c(name), v)
  end subroutine query_double

  subroutine query_logical (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    logical :: v
    integer(c_int) :: i
    call amrex_parmparse_query_bool (this%p, string_f_to_c(name), i)
    v = i.eq.1
  end subroutine query_logical

  subroutine query_string (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    character(*), intent(inout) :: v

    ! temporary string for passing back and forth to C -- include NULL
    character(c_char), dimension(len(v)+1) :: v_pass

    call amrex_parmparse_query_string (this%p, string_f_to_c(name), v_pass, len(v)+1)

    ! convert to Fortran string
    v = string_c_to_f(v_pass)
  end subroutine query_string

end module amrex_parmparse_module
