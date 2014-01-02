module ml_layout_module

  use layout_module
  use multifab_module
  use ml_boxarray_module

  implicit none

  type ml_layout
     integer                  :: dim    = 0
     integer                  :: nlevel = 0
     type(ml_boxarray)        :: mba
     type(layout)   , pointer ::    la(:) => Null()
     type(lmultifab), pointer ::  mask(:) => Null() ! cell-centered mask
     logical        , pointer :: pmask(:) => Null() ! periodic mask
     !
     ! Have to destroy all layouts, not just the coarsest one, if we're
     ! build using ml_layout_build_n().
     !
     logical                  :: destroy_all_layouts = .false.
  end type ml_layout

  interface build
     module procedure ml_layout_build
     module procedure ml_layout_build_n
     module procedure ml_layout_build_la_array
  end interface

  interface destroy
     module procedure ml_layout_destroy
  end interface

  interface operator(.eq.)
     module procedure ml_layout_equal
  end interface
  interface operator(.ne.)
     module procedure ml_layout_not_equal
  end interface

  interface print
     module procedure ml_layout_print
  end interface

  interface nlevels
     module procedure ml_layout_nlevels
  end interface

  interface nboxes
     module procedure ml_layout_nboxes
  end interface

  interface get_box
     module procedure ml_layout_get_box
  end interface

  interface built_q
     module procedure ml_layout_built_q
  end interface

contains

  function ml_layout_built_q(mla) result(r)
    logical :: r
    type(ml_layout), intent(in) :: mla
    r = associated(mla%la)
  end function ml_layout_built_q

  function ml_layout_nlevels(mla) result(r)
    integer :: r
    type(ml_layout), intent(in) :: mla
    r = mla%nlevel
  end function ml_layout_nlevels

  function ml_layout_nboxes(mla, lev) result(r)
    integer :: r
    type(ml_layout), intent(in) :: mla
    integer, intent(in) :: lev
    r = nboxes(mla%mba, lev)
  end function ml_layout_nboxes

  function ml_layout_equal(mla1, mla2) result(r)
    logical :: r
    type(ml_layout), intent(in) :: mla1, mla2
    r = associated(mla1%la, mla2%la)
  end function ml_layout_equal
  
  function ml_layout_not_equal(mla1, mla2) result(r)
    logical :: r
    type(ml_layout), intent(in) :: mla1, mla2
    r = .not. associated(mla1%la, mla2%la)
  end function ml_layout_not_equal

  function ml_layout_get_layout(mla, n) result(r)
    type(layout) :: r
    type(ml_layout), intent(in) :: mla
    integer, intent(in) :: n
    r = mla%la(n)
  end function ml_layout_get_layout

  function ml_layout_get_pd(mla, n) result(r)
    type(box) :: r
    type(ml_layout), intent(in) :: mla
    integer, intent(in) :: n
    r = ml_boxarray_get_pd(mla%mba, n)
  end function ml_layout_get_pd

  function ml_layout_get_box(mla, lev, n) result(r)
    type(box) :: r
    type(ml_layout), intent(in) :: mla
    integer, intent(in) :: n, lev
    r = get_box(mla%la(lev), n)
  end function ml_layout_get_box

  subroutine ml_layout_build_n(mla, nlevel, dm)
    type(ml_layout), intent(out) :: mla
    integer, intent(in) :: nlevel, dm

    mla%nlevel              = nlevel
    mla%dim                 = dm
    mla%destroy_all_layouts = .true.
    allocate(mla%pmask(mla%dim))
    allocate(mla%la(mla%nlevel))
    allocate(mla%mask(mla%nlevel-1))
    call build(mla%mba, nlevel, dm)
  end subroutine ml_layout_build_n

  subroutine ml_layout_build_la_array(mla, mba, la_array, pmask, nlevel)

    type(ml_layout  ), intent(  out) :: mla
    type(ml_boxarray), intent(in   ) :: mba
    type(   layout  ), intent(inout) :: la_array(:)
    integer,           intent(in   ) :: nlevel
    logical                        :: pmask(:)

    type(boxarray) :: bac
    integer        :: n

    mla%destroy_all_layouts = .true.

    mla%nlevel = nlevel
    mla%dim    = get_dim(mba)

    ! Copy only nlevel levels of the mba
    call build(mla%mba,nlevel,mla%dim)

    mla%mba%pd(1:nlevel) = mba%pd(1:nlevel)
    do n = 1, mla%nlevel-1
      mla%mba%rr(n,:) = mba%rr(n,:)
    end do

    do n = 1, mla%nlevel
      call copy(mla%mba%bas(n),mba%bas(n))
    end do

    ! Build the pmask
    allocate(mla%pmask(mla%dim))
    mla%pmask  = pmask

    ! Point to the existing la_array(:)
    allocate(mla%la(mla%nlevel))
    do n = 1, nlevel
       mla%la(n) = la_array(n)
    end do

    allocate(mla%mask(mla%nlevel-1))

    do n = mla%nlevel-1,  1, -1
       call lmultifab_build(mla%mask(n), mla%la(n), nc = 1, ng = 0)
       call setval(mla%mask(n), val = .TRUE.)
       call copy(bac, mba%bas(n+1))
       call boxarray_coarsen(bac, mba%rr(n,:))
       call setval(mla%mask(n), .false., bac)
       call destroy(bac)
    end do

  end subroutine ml_layout_build_la_array

  subroutine ml_layout_build(mla, mba, pmask)

    type(ml_layout)  , intent(inout) :: mla
    type(ml_boxarray), intent(in   ) :: mba
    logical, optional                :: pmask(:)

    type(boxarray) :: bac
    integer :: n
    logical :: lpmask(mba%dim)

    lpmask = .false.; if (present(pmask)) lpmask = pmask
    allocate(mla%pmask(mba%dim))
    mla%pmask  = lpmask

    mla%nlevel = mba%nlevel
    mla%dim    = mba%dim

    call copy(mla%mba, mba)
    allocate(mla%la(mla%nlevel))
    allocate(mla%mask(mla%nlevel-1))
    call build(mla%la(1), mba%bas(1), mba%pd(1), pmask=lpmask)
    do n = 2, mla%nlevel
       call layout_build_pn(mla%la(n), mla%la(n-1), mba%bas(n), mba%rr(n-1,:))
    end do
    do n = mla%nlevel-1,  1, -1
       call lmultifab_build(mla%mask(n), mla%la(n), nc = 1, ng = 0)
       call setval(mla%mask(n), val = .TRUE.)
       call copy(bac, mba%bas(n+1))
       call boxarray_coarsen(bac, mba%rr(n,:))
       call setval(mla%mask(n), .false., bac)
       call destroy(bac)
    end do

  end subroutine ml_layout_build

  subroutine ml_layout_restricted_build(mla, mba, nlevs, pmask)

    ! this subroutine is the same thing as ml_layout_build except that
    ! the mla will only have nlevs instead of mba%nlevel

    type(ml_layout)  , intent(inout) :: mla
    type(ml_boxarray), intent(in   ) :: mba
    integer          , intent(in   ) :: nlevs
    logical, optional                :: pmask(:)

    type(boxarray) :: bac
    integer :: n
    logical :: lpmask(mba%dim)

    lpmask = .false.; if (present(pmask)) lpmask = pmask
    allocate(mla%pmask(mba%dim))
    mla%pmask  = lpmask

    mla%nlevel = nlevs
    mla%dim    = mba%dim

!   Have to copy only nlevs of the mba
!   Replace 
!   call copy(mla%mba, mba)
!   by these lines
    call build(mla%mba,nlevs,mla%dim)
    mla%mba%pd(1:nlevs) = mba%pd(1:nlevs)
    do n = 1, mla%nlevel-1
      mla%mba%rr(n,:) = mba%rr(n,:)
    end do
    do n = 1, mla%nlevel
      call copy(mla%mba%bas(n),mba%bas(n))
    end do

    allocate(mla%la(mla%nlevel))
    allocate(mla%mask(mla%nlevel-1))
    call build(mla%la(1), mba%bas(1), mba%pd(1), pmask=lpmask)
    do n = 2, mla%nlevel
       call layout_build_pn(mla%la(n), mla%la(n-1), mba%bas(n), mba%rr(n-1,:))
    end do
    do n = mla%nlevel-1,  1, -1
       call lmultifab_build(mla%mask(n), mla%la(n), nc = 1, ng = 0)
       call setval(mla%mask(n), val = .TRUE.)
       call copy(bac, mba%bas(n+1))
       call boxarray_coarsen(bac, mba%rr(n,:))
       call setval(mla%mask(n), .false., bac)
       call destroy(bac)
    end do

  end subroutine ml_layout_restricted_build

  subroutine ml_layout_destroy(mla)
    type(ml_layout), intent(inout) :: mla
    integer :: n
    do n = 1, mla%nlevel-1
       if (built_q(mla%mask(n))) call destroy(mla%mask(n))
    end do
    call destroy(mla%mba)

    if ( mla%destroy_all_layouts ) then
       do n = 1, mla%nlevel
          call destroy(mla%la(n))
       end do
    else
       !
       ! We need only delete the coarsest level layout
       ! since it 'owns' the refined levels.
       !
       call destroy(mla%la(1))
    end if

    deallocate(mla%la, mla%mask)
    mla%dim = 0
    mla%nlevel = 0
    deallocate(mla%pmask)
  end subroutine ml_layout_destroy

  subroutine ml_layout_print(mla, str, unit, skip)
    use bl_IO_module
    type(ml_layout), intent(in) :: mla
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: skip
    integer :: i, j
    integer :: un
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt = '("MLLAYOUT[(*")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(" ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') mla%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NLEVEL  = ",i2)') mla%nlevel
    call unit_skip(un, skip)
    write(unit=un, fmt='(" *) {")')
    do i = 1, mla%nlevel
       call unit_skip(un, unit_get_skip(skip)+1)
       write(unit=un, fmt = '("(* LEVEL ", i2)') i
       call unit_skip(un, unit_get_skip(skip)+1)
       write(unit=un, fmt = '(" PD = ")', advance = 'no')
       call print(mla%mba%pd(i), unit=un, advance = 'NO')
       write(unit=un, fmt = '(" *) {")')
       do j = 1, nboxes(mla%mba%bas(i))
           call unit_skip(un, unit_get_skip(skip)+2)
           write(unit=un, fmt = '("{")', advance = 'no')
           call print(get_box(mla%mba%bas(i),j), unit = unit, advance = 'NO')
           write(unit=un, fmt = '(", ", I0, "}")', advance = 'no') get_proc(mla%la(i), j)
           if ( j == nboxes(mla%mba%bas(i)) ) then
              call unit_skip(un, unit_get_skip(skip)+1)
              write(unit=un, fmt = '("}")')
           else
              write(unit=un, fmt = '(",")')
           end if
       end do
       if ( i == mla%nlevel ) then
          call unit_skip(un, skip)
          write(unit=un, fmt = '("}]")')
       else
          write(unit=un, fmt = '(",")')
       end if
    end do
  end subroutine ml_layout_print

end module ml_layout_module
