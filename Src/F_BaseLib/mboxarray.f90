module mboxarray_module

  use boxarray_module

  implicit none

  type mboxarray
     integer :: dim = 0
     integer :: nlevel = 0
     integer, pointer :: rr(:,:) => Null()
     type(boxarray), pointer :: bas(:) => Null()
     type(box), pointer :: pd(:) => Null()
  end type mboxarray

  interface destroy
     module procedure mboxarray_destroy
  end interface

  interface build
     module procedure mboxarray_build_ba
     module procedure mboxarray_build_n
  end interface

  interface print
     module procedure mboxarray_print
  end interface

  interface mboxarray_maxsize
     module procedure mboxarray_maxsize_i
     module procedure mboxarray_maxsize_v
  end interface

  type(mem_stats), private, save :: mboxarray_ms

contains

  subroutine mboxarray_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    mboxarray_ms = ms
  end subroutine mboxarray_set_mem_stats

  function mboxarray_mem_stats() result(r)
    type(mem_stats) :: r
    r = mboxarray_ms
  end function mboxarray_mem_stats

  subroutine mboxarray_build_n(mba, nlevel, dim)
    type(mboxarray), intent(out) :: mba
    integer, intent(in) :: nlevel
    integer, intent(in), optional :: dim
    mba%nlevel = nlevel
    allocate(mba%bas(nlevel), mba%pd(nlevel))
    if ( present(dim) ) then
       call mboxarray_alloc_rr(mba, dim)
    end if
    call mem_stats_alloc(mboxarray_ms)
  end subroutine mboxarray_build_n

  subroutine mboxarray_alloc_rr(mba, dim)
    type(mboxarray), intent(inout) :: mba
    integer, intent(in) :: dim
    if ( .not. associated(mba%bas) ) then
       call bl_error("MBOXARRAY_ALLOC_RR: not allocated bas")
    end if
    if ( associated(mba%rr) ) then
       call bl_error("MBOXARRAY_ALLOC_RR: rr already allocated")
    end if
    mba%dim = dim
    allocate(mba%rr(mba%nlevel-1,mba%dim))
  end subroutine mboxarray_alloc_rr

  subroutine mboxarray_build_ba(mba, ba, pd)
    type(mboxarray), intent(out) :: mba
    type(boxarray), intent(in) :: ba
    type(box), intent(in), optional :: pd

    mba%nlevel = 1
    mba%dim = ba%dim
    allocate(mba%rr(0,mba%dim), mba%bas(1), mba%pd(1))
    call boxarray_build_copy(mba%bas(1), ba)
    if ( present(pd) ) then
       mba%pd(1) = pd
    else
       mba%pd(1) = boxarray_bbox(ba)
    end if
    call mem_stats_alloc(mboxarray_ms)
  end subroutine mboxarray_build_ba

  subroutine mboxarray_destroy(mba)
    type(mboxarray), intent(inout) :: mba
    integer i
    if ( associated(mba%bas) ) then
       do i = 1, mba%nlevel
          call boxarray_destroy(mba%bas(i))
       end do
       deallocate(mba%bas, mba%rr, mba%pd)
       call mem_stats_dealloc(mboxarray_ms)
    end if
  end subroutine mboxarray_destroy

  function mboxarray_get_boxarray(mba, n) result(r)
    type(boxarray) :: r
    type(mboxarray), intent(in) :: mba
    integer, intent(in) :: n
    r = mba%bas(n)
  end function mboxarray_get_boxarray

  function mboxarray_get_box(mba, n, i) result(r)
    type(box) :: r
    type(mboxarray), intent(in) :: mba
    integer, intent(in) :: n, i
    r = get_box(mba%bas(n),i)
  end function mboxarray_get_box

  function mboxarray_get_pd(mba, n) result(r)
    type(box) :: r
    type(mboxarray), intent(in) :: mba
    integer, intent(in) :: n
    r = mba%pd(n)
  end function mboxarray_get_pd

  function mboxarray_refrat_n(mba, n) result(r)
    type(mboxarray), intent(in) :: mba
    integer, intent(in) :: n
    integer :: r(mba%dim)
    if ( n < 1 .or. n >= mba%nlevel) &
         call bl_error("MBOXARRAY_REFRAT_N: out of bounds: ", n)
    r = mba%rr(n,:)
  end function mboxarray_refrat_n

  function mboxarray_refrat(mba) result(r)
    type(mboxarray), intent(in) :: mba
    integer :: r(mba%nlevel-1,mba%dim)
    r = mba%rr
  end function mboxarray_refrat

  subroutine mboxarray_maxsize_i(mba, chunk)
    type(mboxarray), intent(inout) :: mba
    integer, intent(in) :: chunk
    integer :: i
    do i = 1, mba%nlevel
       call boxarray_maxsize(mba%bas(i), chunk)
    end do
  end subroutine mboxarray_maxsize_i
  subroutine mboxarray_maxsize_v(mba, chunk)
    type(mboxarray), intent(inout) :: mba
    integer, intent(in) :: chunk(:)
    integer :: i
    do i = 1, mba%nlevel
       call boxarray_maxsize(mba%bas(i), chunk)
    end do
  end subroutine mboxarray_maxsize_v

  function mboxarray_clean(mba) result(r)
    logical :: r
    type(mboxarray), intent(in) :: mba
    integer i
    do i = 1, mba%nlevel
       if ( .not. boxarray_clean(mba%bas(i)%bxs) ) then
          r = .false.
          return
       end if
    end do
    r = .true.
  end function mboxarray_clean

  function mboxarray_properly_nested_new(mba, nproper) result(r)
    logical :: r
    type(mboxarray), intent(in) :: mba
    integer, intent(in), optional :: nproper
    type(boxarray) :: ba
    integer :: i, lnp
    lnp = 0; if ( present(nproper) ) lnp = nproper
    do i = 2, mba%nlevel
       call boxarray_pn_domain_bx_v(ba, mba%bas(i-1), mba%pd(i), lnp, mba%rr(i-1,:))
       call boxarray_diff(ba, mba%bas(i))
       call boxarray_intersection(ba, mba%pd(i))
       if ( .not. empty(ba) ) then
          call boxarray_destroy(ba)
          r = .false.
          return
       end if
       call boxarray_destroy(ba)
    end do
    r = .true.
  end function mboxarray_properly_nested_new

  function mboxarray_properly_nested(mba, ng) result(r)
    logical :: r
    type(mboxarray), intent(in) :: mba
    integer, intent(in), optional :: ng
    type(boxarray) :: ba
    integer :: i, lng
    lng = 0; if ( present(ng) ) lng = ng
    do i = 1, mba%nlevel-1
       call boxarray_build_copy(ba, mba%bas(i+1))
       call boxarray_coarsen(ba, mba%rr(i,:))
       call boxarray_diff(ba, mba%bas(i))
       call boxarray_intersection(ba, mba%pd(i))
       if ( .not. empty(ba) ) then
          call boxarray_destroy(ba)
          r = .false.
          return
       end if
       call boxarray_destroy(ba)
    end do
    r = .true.
  end function mboxarray_properly_nested

  subroutine mboxarray_print(mba, str, unit, skip)
    use bl_IO_module
    type(mboxarray), intent(in) :: mba
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: skip
    integer :: i
    integer :: un
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt = '("MBOXARRAY")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') mba%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NLEVEL  = ",i2)') mba%nlevel
    do i = 1, mba%nlevel
       call unit_skip(un, unit_get_skip(skip)+1)
       write(unit=un, fmt = '("LEVEL ", i2, " PD ")', advance = 'no') i
       call print(mba%pd(i), unit=un)
       call print(mba%bas(i), skip = unit_get_skip(skip) + 2)
    end do
  end subroutine mboxarray_print

end module mboxarray_module

