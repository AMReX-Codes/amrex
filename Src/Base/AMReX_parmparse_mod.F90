module amrex_parmparse_module

  use iso_c_binding
  use amrex_string_module
  use amrex_error_module
  use amrex_fort_module

  implicit none

  character(len=:), allocatable, public :: amrex_namelist

  private

  public :: amrex_parmparse_build, amrex_parmparse_destroy, amrex_init_namelist, amrex_finalize_namelist

  type, public :: amrex_parmparse
     logical     :: owner = .false.
     type(c_ptr) :: p = c_null_ptr
   contains
     generic :: assignment(=) => amrex_parmparse_assign  ! shallow copy
     generic :: get           => get_int, get_real, get_logical, get_string
     generic :: query         => query_int, query_real, query_logical, query_string
#if defined(__GFORTRAN__) && (__GNUC__ <= 4)
     generic :: getarr        => get_intarr, get_realarr
     generic :: queryarr      => query_intarr, query_realarr
#else
     generic :: getarr        => get_intarr, get_realarr, get_stringarr
     generic :: queryarr      => query_intarr, query_realarr, query_stringarr
#endif
     generic :: add           => add_int, add_real, add_logical, add_string
     generic :: addarr        => add_intarr, add_realarr, add_stringarr
     procedure, private :: amrex_parmparse_assign
     procedure, private :: get_int
     procedure, private :: get_real
     procedure, private :: get_logical
     procedure, private :: get_string
     procedure, private :: get_intarr
     procedure, private :: get_realarr
     procedure, private :: get_stringarr
     procedure, private :: query_int
     procedure, private :: query_real
     procedure, private :: query_logical
     procedure, private :: query_string
     procedure, private :: query_intarr
     procedure, private :: query_realarr
     procedure, private :: query_stringarr
     procedure, private :: add_int
     procedure, private :: add_real
     procedure, private :: add_logical
     procedure, private :: add_string
     procedure, private :: add_intarr
     procedure, private :: add_realarr
     procedure, private :: add_stringarr
#if (!defined(__GFORTRAN__) || (__GNUC__ > 4)) && (!defined(__ibmxl__))
     final :: amrex_parmparse_destroy
#endif
  end type amrex_parmparse

  ! interfaces to cpp functions

  interface
     subroutine amrex_new_parmparse (pp, name) bind(c)
       import
       implicit none
       type(c_ptr) :: pp
       character(kind=c_char), intent(in) :: name(*)
     end subroutine amrex_new_parmparse

     subroutine amrex_delete_parmparse (pp) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
     end subroutine amrex_delete_parmparse

     function amrex_parmparse_get_counts (pp, name) bind(c)
       import
       implicit none
       integer(c_int) :: amrex_parmparse_get_counts
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)       
     end function amrex_parmparse_get_counts

     subroutine amrex_parmparse_get_int (pp, name, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end subroutine amrex_parmparse_get_int

     subroutine amrex_parmparse_get_real (pp, name, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       real(amrex_real) :: v
     end subroutine amrex_parmparse_get_real

     subroutine amrex_parmparse_get_bool (pp, name, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end subroutine amrex_parmparse_get_bool

     subroutine amrex_parmparse_get_string (pp, name, v,  len) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       type(c_ptr), intent(inout) :: v
       integer, intent(out) :: len
     end subroutine amrex_parmparse_get_string

     subroutine amrex_parmparse_delete_cp_char (v, len) bind(c)
       import
       implicit none
       type(c_ptr) :: v(*)
       integer(c_int), value :: len
     end subroutine amrex_parmparse_delete_cp_char

     subroutine amrex_parmparse_get_intarr (pp, name, v, n) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int) :: v(*)
       integer(c_int), value :: n
     end subroutine amrex_parmparse_get_intarr

     subroutine amrex_parmparse_get_realarr (pp, name, v, n) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       real(amrex_real) :: v(*)
       integer(c_int), value :: n
     end subroutine amrex_parmparse_get_realarr

     subroutine amrex_parmparse_get_stringarr (pp, name, v, sv, n) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int), value :: n
       type(c_ptr), intent(inout) :: v(n)
       integer(c_int), intent(inout) :: sv(n)
     end subroutine amrex_parmparse_get_stringarr

     integer function amrex_parmparse_query_int (pp, name, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end function amrex_parmparse_query_int

     integer function amrex_parmparse_query_real (pp, name, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       real(amrex_real) :: v
     end function amrex_parmparse_query_real

     integer function amrex_parmparse_query_bool (pp, name, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int) :: v
     end function amrex_parmparse_query_bool

     integer function amrex_parmparse_query_string (pp, name, v, len) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       type(c_ptr), intent(inout) :: v
       integer(c_int), intent(out) :: len
     end function amrex_parmparse_query_string

     subroutine amrex_parmparse_add_int (pp, name, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int), value :: v
     end subroutine amrex_parmparse_add_int

     subroutine amrex_parmparse_add_real (pp, name, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       real(amrex_real), value :: v
     end subroutine amrex_parmparse_add_real

     subroutine amrex_parmparse_add_bool (pp, name, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int), value :: v
     end subroutine amrex_parmparse_add_bool

     subroutine amrex_parmparse_add_string (pp, name, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       character(kind=c_char), intent(in) :: v(*)
     end subroutine amrex_parmparse_add_string

     subroutine amrex_parmparse_add_intarr (pp, name, v, n) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int), intent(in) :: v(*)
       integer(c_int), value :: n
     end subroutine amrex_parmparse_add_intarr

     subroutine amrex_parmparse_add_realarr (pp, name, v, n) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       real(amrex_real), intent(in) :: v(*)
       integer(c_int), value :: n
     end subroutine amrex_parmparse_add_realarr

     subroutine amrex_parmparse_add_stringarr (pp, name, v, n) bind(c)
       import
       implicit none
       type(c_ptr), value :: pp
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int), value :: n
       character(kind=c_char), intent(in) :: v(*)
     end subroutine amrex_parmparse_add_stringarr
  end interface

contains

  subroutine amrex_parmparse_build (pp, name)
    type(amrex_parmparse) :: pp
    character(*), intent(in), optional :: name
    pp%owner = .true.
    if (present(name)) then
       call amrex_new_parmparse(pp%p, amrex_string_f_to_c(name))
    else
       call amrex_new_parmparse(pp%p, amrex_c_null_char_array)
    end if
  end subroutine amrex_parmparse_build

  subroutine amrex_parmparse_destroy (this)
    type(amrex_parmparse) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_delete_parmparse(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_parmparse_destroy

  subroutine amrex_parmparse_assign (dst, src)
    class(amrex_parmparse), intent(inout) :: dst
    type (amrex_parmparse), intent(in   ) :: src
    call amrex_parmparse_destroy(dst)
    dst%owner = .false.
    dst%p = src%p
  end subroutine amrex_parmparse_assign

  subroutine get_int (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(len=*), intent(in) :: name
    integer :: v
    call amrex_parmparse_get_int (this%p, amrex_string_f_to_c(name), v)
  end subroutine get_int

  subroutine get_real (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    real(amrex_real) :: v
    call amrex_parmparse_get_real (this%p, amrex_string_f_to_c(name), v)
  end subroutine get_real

  subroutine get_logical (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    logical :: v
    integer(c_int) :: i
    call amrex_parmparse_get_bool (this%p, amrex_string_f_to_c(name), i)
    v = i.eq.1
  end subroutine get_logical

  subroutine get_string (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    character(len=:), allocatable, intent(inout) :: v
    ! temporary string for passing back and forth to C -- include NULL
    type(c_ptr) :: cp_pass
    integer :: n_pass
    character(kind=c_char), pointer :: cc(:)
    call amrex_parmparse_get_string (this%p, amrex_string_f_to_c(name), cp_pass, n_pass)
    if (allocated(v)) deallocate(v)
    allocate(character(len=n_pass-1)::v)
    call c_f_pointer(cp_pass, cc, [n_pass])
    v = amrex_string_c_to_f(cc)
    cc => null()
    call amrex_parmparse_delete_cp_char([cp_pass],1)
  end subroutine get_string

  subroutine get_intarr (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, allocatable, intent(inout) :: v(:)
    integer :: n
    n = amrex_parmparse_get_counts(this%p, amrex_string_f_to_c(name))
    if (n .gt. 0) then
       if (allocated(v)) deallocate(v)
       allocate(v(n))
       call amrex_parmparse_get_intarr (this%p, amrex_string_f_to_c(name), v, n)
    else
       call amrex_abort("amrex_parmparse: get_intarr failed to get "//name)
    end if
  end subroutine get_intarr

  subroutine get_realarr (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(len=*), intent(in) :: name
    real(amrex_real), allocatable, intent(inout) :: v(:)
    integer :: n
    n = amrex_parmparse_get_counts(this%p, amrex_string_f_to_c(name))
    if (n .gt. 0) then
       if (allocated(v)) deallocate(v)
       allocate(v(n))
       call amrex_parmparse_get_realarr (this%p, amrex_string_f_to_c(name), v, n)
    else
       call amrex_abort("amrex_parmparse: get_realarr failed to get "//name)
    end if
  end subroutine get_realarr

  subroutine get_stringarr (this, name, v)
    class(amrex_parmparse), intent(in) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(inout) :: v(:)
    integer :: n, i, lenmax
    type(c_ptr), allocatable :: cp_pass(:)
    integer, allocatable :: n_pass(:)
    character(kind=c_char), pointer :: cc(:)
    n = amrex_parmparse_get_counts(this%p, amrex_string_f_to_c(name))
    if (n .gt. 0) then
       allocate(cp_pass(n))
       allocate(n_pass(n))
       call amrex_parmparse_get_stringarr(this%p, amrex_string_f_to_c(name), cp_pass, n_pass, n)
       lenmax = maxval(n_pass)-1
       if (allocated(v)) deallocate(v)
       allocate(character(len=lenmax)::v(n))
       do i = 1, n
          call c_f_pointer(cp_pass(i), cc, [n_pass(i)])
          v(i) = amrex_string_c_to_f(cc)
          cc => null()
       end do
       call amrex_parmparse_delete_cp_char(cp_pass, n)
    else
       call amrex_abort("amrex_parmparse: get_stringarr failed to get "//name)
    end if
  end subroutine get_stringarr

  subroutine query_int (this, name, v, flag)
    class(amrex_parmparse), intent(in) :: this
    character(len=*), intent(in) :: name
    logical, optional, intent(out) :: flag
    integer :: v, iflag
    iflag = amrex_parmparse_query_int (this%p, amrex_string_f_to_c(name), v)
    if (present(flag)) flag = iflag.ne.0
  end subroutine query_int

  subroutine query_real (this, name, v, flag)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    logical, optional, intent(out) :: flag
    real(amrex_real) :: v
    integer :: iflag
    iflag = amrex_parmparse_query_real (this%p, amrex_string_f_to_c(name), v)
    if (present(flag)) flag = iflag.ne.0    
  end subroutine query_real

  subroutine query_logical (this, name, v, flag)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    logical, optional, intent(out) :: flag
    logical :: v
    integer(c_int) :: i, iflag
    iflag = amrex_parmparse_query_bool (this%p, amrex_string_f_to_c(name), i)
    if (iflag.eq.1) then
       v = i.eq.1
    end if
    if (present(flag)) flag = iflag.ne.0    
  end subroutine query_logical

  subroutine query_string (this, name, v, flag)
    class(amrex_parmparse), intent(in) :: this
    character(*), intent(in) :: name
    character(len=:), allocatable, intent(inout) :: v
    logical, optional, intent(out) :: flag
    ! temporary string for passing back and forth to C -- include NULL
    type(c_ptr) :: cp_pass
    integer :: n_pass, iflag
    character(kind=c_char), pointer :: cc(:)
    iflag = amrex_parmparse_query_string (this%p, amrex_string_f_to_c(name), cp_pass, n_pass)
    if (n_pass > 1) then
       if (allocated(v)) deallocate(v)
       allocate(character(len=n_pass-1)::v)
       call c_f_pointer(cp_pass, cc, [n_pass])
       v = amrex_string_c_to_f(cc)
       cc => null()
    end if
    call amrex_parmparse_delete_cp_char([cp_pass],1)
    if (present(flag)) flag = iflag.ne.0
  end subroutine query_string

  subroutine query_intarr (this, name, v, flag)
    class(amrex_parmparse), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, allocatable, intent(inout) :: v(:)
    logical, optional, intent(out) :: flag
    integer :: n
    n = amrex_parmparse_get_counts(this%p, amrex_string_f_to_c(name))
    if (n .gt. 0) then
       call amrex_parmparse_get_intarr (this%p, amrex_string_f_to_c(name), v, n)
    end if
    if (present(flag)) flag = n.gt.0
  end subroutine query_intarr

  subroutine query_realarr (this, name, v, flag)
    class(amrex_parmparse), intent(in) :: this
    character(len=*), intent(in) :: name
    real(amrex_real), allocatable, intent(inout) :: v(:)
    logical, optional, intent(out) :: flag
    integer :: n
    n = amrex_parmparse_get_counts(this%p, amrex_string_f_to_c(name))
    if (n .gt. 0) then
       call amrex_parmparse_get_realarr (this%p, amrex_string_f_to_c(name), v, n)
    end if
    if (present(flag)) flag = n.gt.0
  end subroutine query_realarr

  subroutine query_stringarr (this, name, v, flag)
    class(amrex_parmparse), intent(in) :: this
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(inout) :: v(:)
    logical, optional, intent(out) :: flag
    integer :: n
    n = amrex_parmparse_get_counts(this%p, amrex_string_f_to_c(name))
    if (n .gt. 0) then
       call this%get_stringarr(name, v)
    end if
    if (present(flag)) flag = n.gt.0
  end subroutine query_stringarr

  subroutine add_int (this, name, v)
    class(amrex_parmparse), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: v
    call amrex_parmparse_add_int(this%p, amrex_string_f_to_c(name), v)
  end subroutine add_int

  subroutine add_real (this, name, v)
    class(amrex_parmparse), intent(inout) :: this
    character(*), intent(in) :: name
    real(amrex_real), intent(in) :: v
    call amrex_parmparse_add_real(this%p, amrex_string_f_to_c(name), v)
  end subroutine add_real

  subroutine add_logical (this, name, v)
    class(amrex_parmparse), intent(inout) :: this
    character(*), intent(in) :: name
    logical, intent(in) :: v
    integer(c_int) :: i
    if (v) then
       i = 1
    else
       i = 0
    end if
    call amrex_parmparse_add_bool(this%p, amrex_string_f_to_c(name), i)
  end subroutine add_logical

  subroutine add_string (this, name, v)
    class(amrex_parmparse), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: v
    call amrex_parmparse_add_string(this%p, amrex_string_f_to_c(name), amrex_string_f_to_c(v))
  end subroutine add_string
  
  subroutine add_intarr (this, name, v)
    class(amrex_parmparse), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: v(:)
    call amrex_parmparse_add_intarr(this%p, amrex_string_f_to_c(name), v, size(v))
  end subroutine add_intarr

  subroutine add_realarr (this, name, v)
    class(amrex_parmparse), intent(inout) :: this
    character(*), intent(in) :: name
    real(amrex_real), intent(in) :: v(:)
    call amrex_parmparse_add_realarr(this%p, amrex_string_f_to_c(name), v, size(v))
  end subroutine add_realarr

  subroutine add_stringarr (this, name, v)
    class(amrex_parmparse), intent(inout) :: this
    character(*), intent(in) :: name
    character(len=*), intent(in) :: v(:)
    integer :: n, m, mi, ni, ntot, i
    character(kind=c_char), allocatable :: cs(:)
    n = size(v)
    m = len(v(1))
    ntot = (m+1)*n
    allocate(cs(ntot))
    i = 1
    do ni = 1, n
       do mi = 1, len_trim(v(ni))
         cs(i) = v(ni)(mi:mi)
         i = i+1
       end do
       cs(i) = c_null_char
       i = i+1
    end do
    call amrex_parmparse_add_stringarr(this%p, amrex_string_f_to_c(name), cs, n)
  end subroutine add_stringarr

  subroutine amrex_init_namelist (cstr) bind(c,name='amrex_init_namelist')
    character(kind=c_char), intent(in) :: cstr(*)
    integer :: i, n, oldn
    character(len=:), allocatable :: tmp
    n = 0
    do while (cstr(n+1) .ne. c_null_char)
       n = n+1
    end do
    if (n > 0) then
       if (allocated(amrex_namelist)) then
          oldn = len(amrex_namelist)
          allocate(character(len=oldn)::tmp)
          tmp = amrex_namelist
          deallocate(amrex_namelist)
          allocate(character(len=oldn+n)::amrex_namelist)
          amrex_namelist(1:oldn) = tmp(1:oldn)
          do i = 1, n
             amrex_namelist(i+oldn:i+oldn) = cstr(i)
          end do
       else
          allocate(character(len=n)::amrex_namelist)
          do i = 1, n
             amrex_namelist(i:i) = cstr(i)
          end do
       end if
    end if
  end subroutine amrex_init_namelist

  subroutine amrex_finalize_namelist () bind(c,name='amrex_finalize_namelist')
    if (allocated(amrex_namelist)) deallocate(amrex_namelist)
  end subroutine amrex_finalize_namelist

end module amrex_parmparse_module
