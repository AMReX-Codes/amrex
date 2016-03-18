module ml_boxarray_module

  use boxarray_module

  implicit none

  type ml_boxarray
     logical                 :: owner   = .true.
     integer                 :: dim     =  0
     integer                 :: nlevel  =  0
     integer, pointer        :: rr(:,:) => Null()
     type(boxarray), pointer :: bas(:)  => Null()
     type(box), pointer      :: pd(:)   => Null()
  end type ml_boxarray

  interface get_dim
     module procedure ml_boxarray_dim
  end interface

  interface destroy
     module procedure ml_boxarray_destroy
  end interface

  interface copy
     module procedure ml_boxarray_copy
  end interface

  interface build
     module procedure ml_boxarray_build_ba
     module procedure ml_boxarray_build_n
     module procedure ml_boxarray_build_alias
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

  interface nboxes
     module procedure ml_boxarray_get_nboxes
  end interface

  interface get_boxarray
     module procedure ml_boxarray_get_boxarray
  end interface

  interface amr_ref_ratio_init
     module procedure amr_ref_ratio_init_s
     module procedure amr_ref_ratio_init_m
  end interface amr_ref_ratio_init

  type(mem_stats), private, save :: ml_boxarray_ms

  integer, allocatable, private, save :: amr_ref_ratio(:,:)
  integer, private, parameter :: amr_ref_ratio_default = 2

  private :: boxarray_intersection_bx

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

  function ml_boxarray_get_nboxes(mba, lev) result(r)
    integer :: r
    type(ml_boxarray), intent(in) :: mba
    integer, intent(in) :: lev
    r = nboxes(mba%bas(lev))
  end function ml_boxarray_get_nboxes

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
       call boxarray_build_copy(mba%bas(i), mbai%bas(i))
    end do
  end subroutine ml_boxarray_copy

  subroutine ml_boxarray_build_n(mba, nlevel, dim)
    type(ml_boxarray), intent(out) :: mba
    integer, intent(in) :: nlevel
    integer, intent(in), optional :: dim
    mba%owner  = .true.
    mba%nlevel = nlevel
    allocate(mba%bas(nlevel), mba%pd(nlevel))
    if ( present(dim) ) then
       call ml_boxarray_alloc_rr(mba, dim)
    end if
    call mem_stats_alloc(ml_boxarray_ms)
  end subroutine ml_boxarray_build_n

  subroutine ml_boxarray_alloc_rr(mba, dim)
    use bl_error_module
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

    mba%owner  = .true.
    mba%nlevel = 1
    mba%dim = get_dim(ba)
    allocate(mba%rr(0,mba%dim), mba%bas(1), mba%pd(1))
    call boxarray_build_copy(mba%bas(1), ba)
    if ( present(pd) ) then
       mba%pd(1) = pd
    else
       mba%pd(1) = boxarray_bbox(ba)
    end if
    call mem_stats_alloc(ml_boxarray_ms)
  end subroutine ml_boxarray_build_ba

  subroutine ml_boxarray_build_alias(mba, mba_orig, min_lev, max_lev)

    ! This subroutine builds a shallow copy of the origianl mba with levels starting
    ! at the original min_lev and ending at the original max_lev.  Note that the new mba
    ! starts at level 1, which corresponds to level min_lev in the original mba

    type(ml_boxarray), intent(inout) :: mba
    type(ml_boxarray), intent(in   ) :: mba_orig
    integer          , intent(in   ) :: min_lev, max_lev

    integer :: i, nlevs

    nlevs = max_lev - min_lev + 1

    mba%owner = .false. !!!

    mba%dim    = mba_orig%dim
    mba%nlevel = nlevs

    allocate(mba%rr(nlevs-1,mba%dim))
    do i = 1, nlevs-1
       mba%rr(i,:) = mba_orig%rr(i-1+min_lev,:)
    end do

    allocate(mba%bas(nlevs))
    do i = 1, nlevs
       mba%bas(i) = mba_orig%bas(i-1+min_lev)
    end do

    allocate(mba%pd(nlevs))
    do i = 1, nlevs
       mba%pd(i) = mba_orig%pd(i-1+min_lev)
    end do

  end subroutine ml_boxarray_build_alias

  subroutine ml_boxarray_destroy(mba)
    type(ml_boxarray), intent(inout) :: mba
    integer i
    if ( associated(mba%bas) ) then
       if (mba%owner) then
          do i = 1, mba%nlevel
             call boxarray_destroy(mba%bas(i))
          end do
          call mem_stats_dealloc(ml_boxarray_ms)
       end if
       deallocate(mba%bas, mba%pd)
       deallocate(mba%rr)
       mba%owner  = .true.
       mba%dim    = 0
       mba%nlevel = 0
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
    use bl_error_module
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
    use bl_error_module
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
    use layout_module
    logical :: r
    type(ml_boxarray), intent(in) :: mba
    type(box), pointer :: pbxs(:)
    integer i
    do i = 1, mba%nlevel
       pbxs => dataptr(mba%bas(i))
       if ( .not. boxarray_clean(pbxs) ) then
          r = .false.
          return
       end if
    end do
    r = .true.
  end function ml_boxarray_clean

  function ml_boxarray_properly_nested(mba, nproper, pmask, min_fine_level, max_fine_level) result(r)
    use layout_module
    use bl_prof_module
    use bl_error_module
    logical                       :: r
    type(ml_boxarray), intent(in) :: mba
    integer, intent(in), optional :: nproper
    logical, intent(in), optional :: pmask(:)
    integer, intent(in), optional :: min_fine_level
    integer, intent(in), optional :: max_fine_level

    integer        :: i, j, k, lnp, cnt, shft(3**mba%dim,mba%dim)
    integer        :: lev_min, lev_max
    logical        ::  lpmask(mba%dim), nodal(mba%dim)
    type(boxarray) :: ba, ba_fine, ba_crse_fine
    type(box     ) :: bx, bxs(3**mba%dim)
    type(layout  ) :: fla
    type(list_box) :: pbl, pieces, leftover, extra

    type(list_box_node),   pointer :: bln
    type(box_intersector), pointer :: bi(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_boxarray_properly_nested")

    if ( present(min_fine_level) ) then
      if ( .not. present(max_fine_level) ) &
        call bl_error('must specify both min_fine_level and max_fine_level or neither')
      if ( max_fine_level .gt. mba%nlevel ) &
        call bl_error('max_fine_level into ml_boxarray_properly_nested is too big')
    end if

    lnp    = 0;       if ( present(nproper) ) lnp    = nproper
    lpmask = .false.; if ( present(pmask) )   lpmask = pmask
    lev_min =          2; if (present(min_fine_level)) lev_min = min_fine_level
    lev_max = mba%nlevel; if (present(max_fine_level)) lev_max = max_fine_level

    nodal = .false.
    !
    ! First make sure grids are all contained in problem domain.
    !
    do i = lev_min-1, lev_max
       if (.not. box_contains(mba%pd(i), boxarray_bbox(mba%bas(i)))) then
          call bl_error("boxarray not contained in problem domain")
       end if
    end do
    !
    ! Check that level 1 grids cover all of the problem domain.
    !
    if (lev_min .eq. 2) then
       call boxarray_build_bx(ba,mba%pd(1))
       call boxarray_diff(ba, mba%bas(1))
       if ( .not. empty(ba) ) then
          call print(mba%pd(1),'Problem Domain')
          call print(mba%bas(1),'Level 1 boxarray')
          call bl_error('Level 1 grids must cover entire domain')
       end if
       call boxarray_destroy(ba)
    end if
    !
    ! We now can assume that the level 1 (lowest level) grids cover the entire domain.
    ! Given this, the level 2 grids are automatically properly nested.
    ! So we start the loop at level 3.
    !
    do i = lev_min, lev_max
       !
       ! This part of the test ignores periodic boundaries.
       !
       call boxarray_build_copy(ba_crse_fine,mba%bas(i-1))
       call boxarray_refine(ba_crse_fine,mba%rr(i-1,:))
       call boxarray_build_copy(ba_fine,mba%bas(i))
       call boxarray_grow(ba_fine, lnp)
       call boxarray_intersection_bx(ba_fine, mba%pd(i))
       r = contains(ba_crse_fine, ba_fine)
       call boxarray_destroy(ba_fine)
       if ( .not. r ) then
          call boxarray_destroy(ba_crse_fine)
          call destroy(bpt)
          return
       end if

       if ( any(lpmask) .and. lnp > 0) then
          !
          ! Collect additional boxes that contribute to periodically 
          ! filling fine ghost cells.
          !
          call layout_build_ba(fla, mba%bas(i), mba%pd(i), pmask = lpmask, mapping = LA_LOCAL)
          ! LA_LOCAL ==> bypass processor distribution calculation.

          do j = 1, nboxes(mba%bas(i))
             bx = get_box(mba%bas(i),j)
             call box_periodic_shift(mba%pd(i), bx, nodal, lpmask, lnp, shft, cnt, bxs)
             do k = 1, cnt
                call push_back(pbl, bxs(k))
             end do
             bln => begin(pbl)
             do while (associated(bln))
                bx =  value(bln)
                bi => layout_get_box_intersector(fla, bx)
                do k = 1, size(bi)
                   call push_back(pieces, bi(k)%bx)
                end do
                deallocate(bi)
                leftover = boxlist_boxlist_diff(bx, pieces)
                call splice(extra, leftover)
                call destroy(pieces)
                bln => next(bln)
             end do
             call destroy(pbl)
          end do

          call layout_destroy(fla)
          !
          ! Check to see if extra boxes are also covered by coarser region.
          !
          if ( size(extra) > 0 ) then
             call boxarray_build_l(ba, extra)
             call destroy(extra)
             r = contains(ba_crse_fine, ba)
             call boxarray_destroy(ba)
             if ( .not. r ) then
                call boxarray_destroy(ba_crse_fine)
                call destroy(bpt)
                return
             end if
          end if
       end if

       call boxarray_destroy(ba_crse_fine)

    end do

    call destroy(bpt)

    r = .true.
    
  end function ml_boxarray_properly_nested

  ! NOTE: this only works as a "contains" operator, it doesn't take a buffer width
  function ml_boxarray_properly_nested_old(mba) result(r)
    logical :: r
    type(ml_boxarray), intent(in) :: mba
    type(boxarray) :: ba
    integer :: i
    do i = 1, mba%nlevel-1
       call boxarray_build_copy(ba, mba%bas(i+1))
       call boxarray_coarsen(ba, mba%rr(i,:))
       call boxarray_diff(ba, mba%bas(i))
       call boxarray_intersection_bx(ba, mba%pd(i))
       if ( .not. empty(ba) ) then
          call boxarray_destroy(ba)
          r = .false.
          return
       end if
       call boxarray_destroy(ba)
    end do
    r = .true.
  end function ml_boxarray_properly_nested_old

  !
  ! This is a very naive implementation that works perfectly here
  ! because bx is the physical domain box that intersects with every boxes in ba.
  !
  subroutine boxarray_intersection_bx(ba, bx)
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bx
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
      ba%bxs(i) = intersection(ba%bxs(i), bx)
    end do
    !$OMP END PARALLEL DO
    call boxarray_simplify(ba)
  end subroutine boxarray_intersection_bx

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
    write(unit=un, fmt = '("MLBOXARRAY[(*")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(" ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') mba%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NLEVEL  = ",i2)') mba%nlevel
    call unit_skip(un, skip)
    write(unit=un, fmt='(" *) {")')
    do i = 1, mba%nlevel
       call unit_skip(un, unit_get_skip(skip)+1)
       write(unit=un, fmt = '("(* LEVEL ", i2)') i
       call unit_skip(un, unit_get_skip(skip)+1)
       write(unit=un, fmt = '(" PD = ")', advance = 'no')
       call print(mba%pd(i), unit=un, advance = 'NO')
       write(unit=un, fmt = '(" *)")')
       call print(mba%bas(i), skip = unit_get_skip(skip) + 2)
       if ( i == mba%nlevel ) then
          call unit_skip(un, skip)
          write(unit=un, fmt = '("}]")')
       else
          write(unit=un, fmt = '(",")')
       end if
    end do
  end subroutine ml_boxarray_print

  subroutine amr_ref_ratio_init_s(max_levs,dim,rr)
    integer, intent(in) :: max_levs, dim, rr
    allocate(amr_ref_ratio(max_levs-1,dim))
    amr_ref_ratio = rr
  end subroutine amr_ref_ratio_init_s

  subroutine amr_ref_ratio_init_m(rr)
    integer, intent(in) :: rr(:,:)
    allocate(amr_ref_ratio(size(rr,1),size(rr,2)))
    amr_ref_ratio = rr
  end subroutine amr_ref_ratio_init_m

  subroutine ml_boxarray_set_ref_ratio(mba)
    use bl_error_module
    type(ml_boxarray), intent(inout) :: mba
    if (.not. associated(mba%rr)) then
       call bl_error("ML_BOXARRAY_SET_REF_RATIO: mba%rr not allocated")
    end if
    if (allocated(amr_ref_ratio)) then
       mba%rr = amr_ref_ratio(1:mba%nlevel-1,1:mba%dim)
    else
       mba%rr = amr_ref_ratio_default
    end if
  end subroutine ml_boxarray_set_ref_ratio

end module ml_boxarray_module

