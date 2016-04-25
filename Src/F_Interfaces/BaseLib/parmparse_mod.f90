
module parmparse_module

  use iso_c_binding
  use string_module

  implicit none

  private

  public :: parmparse_build

  type, public :: ParmParse
     type(c_ptr) :: p = c_null_ptr
   contains
     generic :: get      => get_int, get_double, get_logical !, get_string
!     generic :: getarr   =>
     generic :: query    => query_int, query_double, query_logical  !, query_string
!     generic :: queryarr =>
     procedure, private :: get_int, get_double, get_logical
     procedure, private :: query_int, query_double, query_logical
     final :: parmparse_destroy
  end type ParmParse

  ! interfaces to cpp functions

  interface
     subroutine fi_new_parmparse (pp, name) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr) :: pp
       character(c_char), intent(in) :: name(*)
     end subroutine fi_new_parmparse

     subroutine fi_delete_parmparse (pp) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
     end subroutine fi_delete_parmparse

     subroutine fi_parmparse_get_int (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end subroutine fi_parmparse_get_int

     subroutine fi_parmparse_get_double (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       real(c_double) :: v
     end subroutine fi_parmparse_get_double

     subroutine fi_parmparse_get_bool (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end subroutine fi_parmparse_get_bool

     subroutine fi_parmparse_query_int (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end subroutine fi_parmparse_query_int

     subroutine fi_parmparse_query_double (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       real(c_double) :: v
     end subroutine fi_parmparse_query_double

     subroutine fi_parmparse_query_bool (pp, name, v) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value :: pp
       character(c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end subroutine fi_parmparse_query_bool
  end interface

contains

  subroutine parmparse_build (pp, name)
    type(ParmParse) :: pp
    character(*), intent(in), optional :: name
    if (present(name)) then
       call fi_new_parmparse(pp%p, string_f_to_c(name))
    else
       call fi_new_parmparse(pp%p, c_char_""//c_null_char)
    end if
  end subroutine parmparse_build

  subroutine parmparse_destroy (this)
    type(ParmParse) :: this
    if (c_associated(this%p)) then
       call fi_delete_parmparse(this%p)
       this%p = c_null_ptr
    end if
  end subroutine parmparse_destroy

  subroutine get_int (this, name, v)
    class(ParmParse), intent(in) :: this
    character(len=*), intent(in) :: name
    integer :: v
    call fi_parmparse_get_int (this%p, string_f_to_c(name), v)
  end subroutine get_int

  subroutine get_double (this, name, v)
    class(ParmParse), intent(in) :: this
    character(*), intent(in) :: name
    double precision :: v
    call fi_parmparse_get_double (this%p, string_f_to_c(name), v)
  end subroutine get_double

  subroutine get_logical (this, name, v)
    class(ParmParse), intent(in) :: this
    character(*), intent(in) :: name
    logical :: v
    integer(c_int) :: i
    call fi_parmparse_get_bool (this%p, string_f_to_c(name), i)
    v = i.eq.1
  end subroutine get_logical

  subroutine query_int (this, name, v)
    class(ParmParse), intent(in) :: this
    character(len=*), intent(in) :: name
    integer :: v
    call fi_parmparse_query_int (this%p, string_f_to_c(name), v)
  end subroutine query_int

  subroutine query_double (this, name, v)
    class(ParmParse), intent(in) :: this
    character(*), intent(in) :: name
    double precision :: v
    call fi_parmparse_query_double (this%p, string_f_to_c(name), v)
  end subroutine query_double

  subroutine query_logical (this, name, v)
    class(ParmParse), intent(in) :: this
    character(*), intent(in) :: name
    logical :: v
    integer(c_int) :: i
    call fi_parmparse_query_bool (this%p, string_f_to_c(name), i)
    v = i.eq.1
  end subroutine query_logical

end module parmparse_module
