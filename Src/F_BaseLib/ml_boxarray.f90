module ml_boxarray_module

  use boxarray_module

  implicit none

  type ml_boxarray
     integer :: dim = 0
     integer :: nlevel = 0
     integer, pointer :: rr(:,:) => Null()
     type(boxarray), pointer :: bas(:) => Null()
     type(box), pointer :: pd(:) => Null()
  end type ml_boxarray

  interface destroy
     module procedure ml_boxarray_destroy
  end interface

  interface copy
     module procedure ml_boxarray_copy
  end interface

  interface build
     module procedure ml_boxarray_build_ba
     module procedure ml_boxarray_build_n
  end interface

  interface print
     module procedure ml_boxarray_print
  end interface

  interface ml_boxarray_maxsize
     module procedure ml_boxarray_maxsize_i
     module procedure ml_boxarray_maxsize_v
  end interface

  interface get_box
     module procedure ml_boxarray_get_box
  end interface

  interface get_pd
     module procedure ml_boxarray_get_pd
  end interface

  interface get_nlevel
     module procedure ml_boxarray_get_nlevel
  end interface

  interface get_boxarray
     module procedure ml_boxarray_get_boxarray
  end interface

  type(mem_stats), private, save :: ml_boxarray_ms

contains

  function ml_boxarray_dim(mba) result(r)
    integer :: r
    type(ml_boxarray), intent(in) :: mba
    r = mba%dim
  end function ml_boxarray_dim

  function ml_boxarray_get_nlevel(mba) result(r)
    integer :: r
    type(ml_boxarray), intent(in) :: mba
    r = mba%nlevel
  end function ml_boxarray_get_nlevel

  subroutine ml_boxarray_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    ml_boxarray_ms = ms
  end subroutine ml_boxarray_set_mem_stats

  function ml_boxarray_mem_stats() result(r)
    type(mem_stats) :: r
    r = ml_boxarray_ms
  end function ml_boxarray_mem_stats

  subroutine ml_boxarray_copy(mba, mbai)
    type(ml_boxarray), intent(out) :: mba
    type(ml_boxarray), intent(in)  :: mbai
    integer :: i
    call build(mba, mbai%nlevel, mbai%dim)
    mba%pd = mbai%pd
    mba%rr = mbai%rr
    do i = 1, mba%nlevel
       call copy(mba%bas(i), mbai%bas(i))
    end do
  end subroutine ml_boxarray_copy

  subroutine ml_boxarray_build_n(mba, nlevel, dim)
    type(ml_boxarray), intent(out) :: mba
    integer, intent(in) :: nlevel
    integer, intent(in), optional :: dim
    mba%nlevel = nlevel
    allocate(mba%bas(nlevel), mba%pd(nlevel))
    if ( present(dim) ) then
       call ml_boxarray_alloc_rr(mba, dim)
    end if
    call mem_stats_alloc(ml_boxarray_ms)
  end subroutine ml_boxarray_build_n

  subroutine ml_boxarray_alloc_rr(mba, dim)
    type(ml_boxarray), intent(inout) :: mba
    integer, intent(in) :: dim
    if ( .not. associated(mba%bas) ) then
       call bl_error("ML_BOXARRAY_ALLOC_RR: not allocated bas")
    end if
    if ( associated(mba%rr) ) then
       call bl_error("ML_BOXARRAY_ALLOC_RR: rr already allocated")
    end if
    mba%dim = dim
    allocate(mba%rr(mba%nlevel-1,mba%dim))
  end subroutine ml_boxarray_alloc_rr

  subroutine ml_boxarray_build_ba(mba, ba, pd)
    type(ml_boxarray), intent(out) :: mba
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
    call mem_stats_alloc(ml_boxarray_ms)
  end subroutine ml_boxarray_build_ba

  subroutine ml_boxarray_destroy(mba)
    type(ml_boxarray), intent(inout) :: mba
    integer i
    if ( associated(mba%bas) ) then
       do i = 1, mba%nlevel
          call boxarray_destroy(mba%bas(i))
       end do
       deallocate(mba%bas, mba%pd)
       deallocate(mba%rr)
       call mem_stats_dealloc(ml_boxarray_ms)
    end if
  end subroutine ml_boxarray_destroy

  function ml_boxarray_get_boxarray(mba, n) result(r)
    type(boxarray) :: r
    type(ml_boxarray), intent(in) :: mba
    integer, intent(in) :: n
    r = mba%bas(n)
  end function ml_boxarray_get_boxarray

  function ml_boxarray_get_box(mba, n, i) result(r)
    type(box) :: r
    type(ml_boxarray), intent(in) :: mba
    integer, intent(in) :: n, i
    r = get_box(mba%bas(n),i)
  end function ml_boxarray_get_box

  function ml_boxarray_get_pd(mba, n) result(r)
    type(box) :: r
    type(ml_boxarray), intent(in) :: mba
    integer, intent(in) :: n
    r = mba%pd(n)
  end function ml_boxarray_get_pd

  function ml_boxarray_refrat_n_d(mba, n, dim) result(r)
    type(ml_boxarray), intent(in) :: mba
    integer, intent(in) :: n
    integer, intent(in) :: dim
    integer :: r
    if ( n < 1 .or. n >= mba%nlevel) &
         call bl_error("ML_BOXARRAY_REFRAT_N: out of bounds: ", n)
    if ( dim < 1 .or. dim > mba%dim) &
         call bl_error("ML_BOXARRAY_REFRAT_N_D: dim out of bounds: ", dim)
    r = mba%rr(n,dim)
  end function ml_boxarray_refrat_n_d

  function ml_boxarray_refrat_n(mba, n) result(r)
    type(ml_boxarray), intent(in) :: mba
    integer, intent(in) :: n
    integer :: r(mba%dim)
    if ( n < 1 .or. n >= mba%nlevel) &
         call bl_error("ML_BOXARRAY_REFRAT_N: out of bounds: ", n)
    r = mba%rr(n,:)
  end function ml_boxarray_refrat_n

  function ml_boxarray_refrat(mba) result(r)
    type(ml_boxarray), intent(in) :: mba
    integer :: r(mba%nlevel-1,mba%dim)
    r = mba%rr
  end function ml_boxarray_refrat

  subroutine ml_boxarray_maxsize_i(mba, chunk)
    type(ml_boxarray), intent(inout) :: mba
    integer, intent(in) :: chunk
    integer :: i
    do i = 1, mba%nlevel
       call boxarray_maxsize(mba%bas(i), chunk)
    end do
  end subroutine ml_boxarray_maxsize_i
  subroutine ml_boxarray_maxsize_v(mba, chunk)
    type(ml_boxarray), intent(inout) :: mba
    integer, intent(in) :: chunk(:)
    integer :: i
    do i = 1, mba%nlevel
       call boxarray_maxsize(mba%bas(i), chunk)
    end do
  end subroutine ml_boxarray_maxsize_v

  function ml_boxarray_clean(mba) result(r)
    logical :: r
    type(ml_boxarray), intent(in) :: mba
    integer i
    do i = 1, mba%nlevel
       if ( .not. boxarray_clean(mba%bas(i)%bxs) ) then
          r = .false.
          return
       end if
    end do
    r = .true.
  end function ml_boxarray_clean

  function ml_boxarray_properly_nested_new(mba, nproper) result(r)
    logical :: r
    type(ml_boxarray), intent(in) :: mba
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
  end function ml_boxarray_properly_nested_new

  function ml_boxarray_properly_nested(mba, ng) result(r)
    logical :: r
    type(ml_boxarray), intent(in) :: mba
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
  end function ml_boxarray_properly_nested

  subroutine ml_boxarray_print(mba, str, unit, skip)
    use bl_IO_module
    type(ml_boxarray), intent(in) :: mba
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: skip
    integer :: i
    integer :: un
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt = '("ML_BOXARRAY")', advance = 'no')
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
  end subroutine ml_boxarray_print

end module ml_boxarray_module

