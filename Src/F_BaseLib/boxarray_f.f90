!!
!! A _BoxArray_ is an array of boxes.
!!
module boxarray_module

  use bl_types
  use box_module
  use list_box_module
  use bl_mem_stat_module

  implicit none

  type boxarray
     private
     integer :: dim = 0
     integer :: nboxes = 0
     type(box), pointer :: bxs(:) => Null()
  end type boxarray

  interface dataptr
     module procedure boxarray_dataptr
  end interface

  interface get_dim
     module procedure boxarray_dim
  end interface

  interface empty
     module procedure boxarray_empty
  end interface

  interface built_q
     module procedure boxarray_built_q
  end interface

  interface copy
     module procedure boxarray_build_copy
     module procedure boxarray_build_copy_l
  end interface

  interface build
     module procedure boxarray_build_v
     module procedure boxarray_build_l
     module procedure boxarray_build_bx
  end interface

  interface destroy
     module procedure boxarray_destroy
  end interface

  interface nboxes
     module procedure boxarray_nboxes
     module procedure boxlist_nboxes
  end interface

  interface volume
     module procedure boxarray_volume
  end interface

  interface dvolume
     module procedure boxarray_dvolume
  end interface

  interface set_box
     module procedure boxarray_set_box
  end interface

  interface get_box
     module procedure boxarray_get_box
  end interface

  interface boxarray_maxsize
     module procedure boxarray_maxsize_i
     module procedure boxarray_maxsize_v
  end interface

  interface boxarray_coarsen
     module procedure boxarray_coarsen_v
     module procedure boxarray_coarsen_i
  end interface

  interface boxarray_refine
     module procedure boxarray_refine_v
     module procedure boxarray_refine_i
  end interface

  interface boxarray_shift
     module procedure boxarray_shift_v
     module procedure boxarray_shift_i
  end interface

  interface boxarray_intersection
     module procedure boxarray_intersection_bx
  end interface

  interface boxarray_grow
     module procedure boxarray_grow_n
     module procedure boxarray_grow_n_f
     module procedure boxarray_grow_n_d_f
     module procedure boxarray_grow_v
     module procedure boxarray_grow_v_f
  end interface

  interface bbox
     module procedure boxarray_bbox
  end interface

  interface print
     module procedure boxarray_print
  end interface

  interface contains
     module procedure boxarray_box_contains
     module procedure boxarray_boxarray_contains
  end interface

  interface equal
     module procedure boxarray_equal
  end interface
  interface operator( .eq. )
     module procedure boxarray_equal
  end interface

  interface not_equal
     module procedure boxarray_not_equal
  end interface
  interface operator( .ne. )
     module procedure boxarray_not_equal
  end interface

  private :: boxarray_maxsize_l
  private :: boxlist_nboxes

  type(mem_stats), private, save :: boxarray_ms

contains

  function boxarray_dataptr(ba) result(r)
    type(boxarray), intent(in) :: ba
    type(box), pointer :: r(:)
    r => ba%bxs
  end function boxarray_dataptr
  
  pure function boxarray_equal(ba1, ba2) result(r)
    type(boxarray), intent(in) :: ba1, ba2
    logical :: r
    r = associated(ba1%bxs, ba2%bxs)
  end function boxarray_equal

  pure function boxarray_not_equal(ba1, ba2) result(r)
    type(boxarray), intent(in) :: ba1, ba2
    logical :: r
    r = .not. associated(ba1%bxs, ba2%bxs)
  end function boxarray_not_equal

  pure function boxarray_same_q(ba1, ba2) result(r)
    type(boxarray), intent(in) :: ba1, ba2
    logical :: r
    integer :: i
    if ( ba1 == ba2 ) then
       r = .true.
    else if ( ba1%dim /= ba2%dim .or. ba1%nboxes /= ba2%nboxes ) then
       r = .false.
    else
       do i = 1, ba1%nboxes
         if ( ba1%bxs(i) /= ba2%bxs(i) ) then
            r = .false.
            return
         end if
       end do
       r = .true.
    end if
  end function boxarray_same_q

  subroutine boxarray_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    boxarray_ms = ms
  end subroutine boxarray_set_mem_stats

  function boxarray_mem_stats() result(r)
    type(mem_stats) :: r
    r = boxarray_ms
  end function boxarray_mem_stats

  pure function boxarray_empty(ba) result(r)
    logical :: r
    type(boxarray), intent(in) :: ba
    r = ba%nboxes == 0
  end function boxarray_empty

  pure function boxarray_built_q(ba) result(r)
    logical :: r
    type(boxarray), intent(in) :: ba
    r = ba%dim /= 0
  end function boxarray_built_q

  pure function boxarray_dim(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer :: r
    r = ba%dim
  end function boxarray_dim

  subroutine boxarray_set_box(ba, i, bx)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: i
    type(box), intent(in) :: bx
    ba%bxs(i) = bx
  end subroutine boxarray_set_box

  pure function boxarray_get_box(ba, i) result(r)
    type(boxarray), intent(in) :: ba
    integer, intent(in) :: i
    type(box) :: r
    r = ba%bxs(i)
  end function boxarray_get_box

  subroutine boxarray_build_copy(ba, ba1)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(boxarray), intent(in) :: ba1
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_COPY: already built")
    if ( .not. built_q(ba1) ) return
    ba%nboxes = size(ba1%bxs)
    allocate(ba%bxs(size(ba1%bxs)))
    ba%bxs = ba1%bxs
    ba%dim = ba1%dim
    call boxarray_verify_dim(ba)
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_copy

  subroutine boxarray_build_copy_l(ba, bl)
    type(boxarray), intent(inout) :: ba
    type(list_box), intent(in) :: bl
    if ( built_q(ba) ) call destroy(ba)
    call boxarray_build_l(ba, bl)
  end subroutine boxarray_build_copy_l

  subroutine boxarray_build_v(ba, bxs, sort)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(box), intent(in), dimension(:) :: bxs
    logical, intent(in), optional :: sort
    logical :: lsort
    
    lsort = .false. ; if (present(sort)) lsort = sort
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_V: already built")
    ba%nboxes = size(bxs)
    allocate(ba%bxs(size(bxs)))
    ba%bxs = bxs
    if ( ba%nboxes > 0 ) then
       ba%dim = ba%bxs(1)%dim
    end if
    call boxarray_verify_dim(ba)
    if (lsort) call boxarray_sort(ba) !! make sure all grids are sorted
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_v

  subroutine boxarray_build_bx(ba, bx)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bx
    
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_BX: already built")
    ba%nboxes = 1
    allocate(ba%bxs(1))
    ba%bxs(1) = bx
    ba%dim = bx%dim
    call boxarray_verify_dim(ba)
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_bx

  subroutine boxarray_build_l(ba, bl, sort)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(list_box), intent(in) :: bl
    logical, intent(in), optional :: sort
    type(list_box_node), pointer :: bln
    logical :: lsort
    integer :: i
    !
    ! Default is to sort.
    !
    lsort = .true. ; if ( present(sort) ) lsort = sort
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_L: already built")
    ba%nboxes = size(bl)
    allocate(ba%bxs(ba%nboxes))
    bln => begin(bl)
    i = 1
    do while (associated(bln))
       ba%bxs(i) = value(bln)
       i = i + 1
       bln=>next(bln)
    end do
    if ( ba%nboxes > 0 ) then
       ba%dim = ba%bxs(1)%dim
    end if
    call boxarray_verify_dim(ba)
    if ( lsort ) call boxarray_sort(ba)
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_l

  subroutine boxarray_destroy(ba)
    type(boxarray), intent(inout) :: ba
    if ( associated(ba%bxs) ) then
       call mem_stats_dealloc(boxarray_ms, ba%nboxes)
       deallocate(ba%bxs) 
       ba%bxs => Null()
    end if
    ba%dim = 0
    ba%nboxes = 0
  end subroutine boxarray_destroy

  subroutine boxarray_sort(ba)
    use sort_box_module
    type(boxarray), intent(inout) :: ba
    call box_sort(ba%bxs)
  end subroutine boxarray_sort

  subroutine boxarray_verify_dim(ba, stat)
    use bl_error_module
    type(boxarray), intent(in) :: ba
    integer, intent(out), optional :: stat
    integer :: i, dm
    if ( present(stat) ) stat = 0
    if ( ba%nboxes < 1 ) return
    dm = ba%dim
    if ( dm == 0 ) then
       dm = ba%bxs(1)%dim
    end if
    if ( dm == 0 ) then
       call bl_error("BOXARRAY_VERIFY_DIM: dim is zero!")
    end if
    do i = 1, ba%nboxes
       if ( ba%dim /= ba%bxs(i)%dim ) then
          if ( present(stat) ) then
             stat = 1
             return
          else
             call bl_error("BOXARRAY_VERIFY_DIM: " // &
                  "ba%dim not equal to some boxes dim: ", ba%dim)
          end if
       end if
    end do
  end subroutine boxarray_verify_dim

  subroutine boxarray_grow_v(ba, rv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:)
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), rv)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_grow_v
  subroutine boxarray_grow_v_f(ba, rv, face)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:), face
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), rv, face)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_grow_v_f
  subroutine boxarray_grow_n(ba, n, allow_empty)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: n
    logical, intent(in), optional :: allow_empty
    integer :: i
    logical :: lallow
    lallow = .false. ; if (present(allow_empty)) lallow = allow_empty
    if (lallow) then
      !$OMP PARALLEL DO
      do i = 1, ba%nboxes
         if (empty(ba%bxs(i))) cycle
         ba%bxs(i) = grow(ba%bxs(i), n)
      end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO
      do i = 1, ba%nboxes
         ba%bxs(i) = grow(ba%bxs(i), n)
      end do
      !$OMP END PARALLEL DO
    endif
  end subroutine boxarray_grow_n
  subroutine boxarray_grow_n_f(ba, n, face)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: n, face
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), n, face)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_grow_n_f
  subroutine boxarray_grow_n_d_f(ba, n, dim, face)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: n, face, dim
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), n, dim, face)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_grow_n_d_f

  subroutine boxarray_nodalize(ba, nodal)
    type(boxarray), intent(inout) :: ba
    logical, intent(in), optional :: nodal(:)
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
       ba%bxs(i) = box_nodalize(ba%bxs(i), nodal)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_nodalize

  pure function boxarray_projectable(ba, rr) result(r)
    logical :: r
    type(boxarray), intent(in) :: ba
    integer, intent(in) :: rr(:)
    integer :: i
    r = .true.
    do i = 1, nboxes(ba)
       if ( .not. box_projectable(ba%bxs(i), rr) ) then
          r = .false.
          exit
       end if
    end do
  end function boxarray_projectable

  subroutine boxarray_coarsen_v_m(ba, cv, mask)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: cv(:)
    logical, intent(in) :: mask(:)
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), cv, mask)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_coarsen_v_m
  subroutine boxarray_coarsen_v(ba, cv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: cv(:)
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), cv)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_coarsen_v
  subroutine boxarray_coarsen_i_m(ba, ci, mask)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ci
    logical, intent(in) :: mask(:)
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), ci, mask)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_coarsen_i_m
  subroutine boxarray_coarsen_i(ba, ci)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ci
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), ci)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_coarsen_i

  subroutine boxarray_shift_v(ba, rv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:)
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
      ba%bxs(i) = shift(ba%bxs(i), rv)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_shift_v
  subroutine boxarray_shift_i(ba, ri)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ri
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
      ba%bxs(i) = shift(ba%bxs(i), ri)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_shift_i
    
  subroutine boxarray_refine_v(ba, rv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:)
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
      ba%bxs(i) = refine(ba%bxs(i), rv)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_refine_v
  subroutine boxarray_refine_i(ba, ri)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ri
    integer :: i
    !$OMP PARALLEL DO
    do i = 1, ba%nboxes
      ba%bxs(i) = refine(ba%bxs(i), ri)
    end do
    !$OMP END PARALLEL DO
  end subroutine boxarray_refine_i
  !
  ! This is a very naive implementation.
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

  subroutine boxarray_box_boundary_n(bao, bx, n)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: n
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_n(bao, baa, n)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_n

  subroutine boxarray_boundary_n_d_f(bao, ba, n, dim, face)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer, intent(in) :: n, face, dim
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, n, dim, face)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_n_d_f
  subroutine boxarray_boundary_n(bao, ba, n)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer,        intent(in)  :: n
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, n)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_n

  pure function boxarray_nboxes(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer :: r
    r = ba%nboxes
  end function boxarray_nboxes

  pure function boxlist_nboxes(bl) result(r)
    type(list_box), intent(in) :: bl
    integer :: r
    r = size(bl)
  end function boxlist_nboxes

  function boxarray_volume(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer(kind=ll_t) :: r
    integer :: i
    r = 0_ll_t
    !$OMP PARALLEL DO REDUCTION(+:r)
    do i = 1, ba%nboxes
       r = r + box_volume(ba%bxs(i))
    end do
    !$OMP END PARALLEL DO
  end function boxarray_volume

  function boxarray_dvolume(ba) result(r)
    type(boxarray), intent(in) :: ba
    real(dp_t) :: r
    integer :: i
    r = 0_dp_t
    !$OMP PARALLEL DO REDUCTION(+:r)
    do i = 1, ba%nboxes
       r = r + box_dvolume(ba%bxs(i))
    end do
    !$OMP END PARALLEL DO
  end function boxarray_dvolume

  pure function boxarray_bbox(ba) result(r)
    type(boxarray), intent(in) :: ba
    type(box) :: r
    integer :: i
    r = nobox(ba%dim)
    do i = 1, ba%nboxes
       r = bbox(r, ba%bxs(i))
    end do
  end function boxarray_bbox

  subroutine boxarray_box_diff(ba, b1, b2)
    type(boxarray), intent(out) :: ba
    type(box), intent(in) :: b1, b2
    type(list_box) :: bl
    bl = boxlist_box_diff(b1, b2)
    call boxarray_build_l(ba, bl)
    call destroy(bl)
  end subroutine boxarray_box_diff

  subroutine boxarray_diff(bao, ba)
    type(boxarray), intent(inout) :: bao
    type(boxarray), intent(in) :: ba
    type(list_box) :: bl, bl1, bl2
    integer :: i
    call build(bl1, ba%bxs)
    do i = 1, bao%nboxes
       bl2 = boxlist_boxlist_diff(bao%bxs(i), bl1)
       call splice(bl, bl2)
    end do
    call boxarray_destroy(bao)
    call boxarray_build_l(bao, bl)
    call destroy(bl)
    call destroy(bl1)
  end subroutine boxarray_diff

  subroutine boxarray_maxsize_i(bxa, chunk) 
    type(boxarray), intent(inout) :: bxa
    integer, intent(in) :: chunk
    integer :: vchunk(bxa%dim) 
    vchunk = chunk
    call boxarray_maxsize_v(bxa, vchunk)
  end subroutine boxarray_maxsize_i

  subroutine boxarray_maxsize_v(bxa, chunk) 
    type(boxarray), intent(inout) :: bxa
    integer, intent(in), dimension(:) :: chunk
    type(list_box) :: bl
    bl = boxarray_maxsize_l(bxa, chunk)
    call boxarray_destroy(bxa)
    call boxarray_build_l(bxa, bl)
    call destroy(bl)
  end subroutine boxarray_maxsize_v

  function boxarray_maxsize_l(bxa, chunk) result(r)
    type(list_box) :: r
    type(boxarray), intent(in) ::  bxa
    integer, intent(in), dimension(:) :: chunk
    integer :: i,k
    type(list_box_node), pointer :: li
    integer :: len(bxa%dim)
    integer :: nl, bs, rt, nblk, sz, ex, ks, ps
    type(box) :: bxr, bxl

    do i = 1, bxa%nboxes
       call push_back(r, bxa%bxs(i))
    end do
    li => begin(r)
    do while ( associated(li) )
       len = extent(value(li))
       do i = 1, bxa%dim
          if ( len(i) > chunk(i) ) then
             rt = 1
             bs = chunk(i)
             nl = len(i)
             do while ( mod(bs,2) == 0 .AND. mod(nl,2) == 0)
                rt = rt * 2
                bs = bs/2
                nl = nl/2
             end do
             nblk = nl/bs
             if ( mod(nl,bs) /= 0 ) nblk = nblk + 1
             sz   = nl/nblk
             ex   = mod(nl,nblk)
             do k = 0, nblk-2
                if ( k < ex ) then
                   ks = (sz+1)*rt
                else
                   ks = sz*rt
                end if
                ps = upb(value(li), i) - ks + 1
                call box_chop(value(li), bxr, bxl, i, ps)
                call set(li, bxr)
                call push_back(r, bxl)
             end do
          end if
       end do
       li => next(li)
    end do

  end function boxarray_maxsize_l

  subroutine boxarray_simplify(bxa)
    type(boxarray), intent(inout) :: bxa
    type(list_box) :: bxl
    call build(bxl, bxa%bxs)
    call boxlist_simplify(bxl)
    call boxarray_destroy(bxa)
    call boxarray_build_l(bxa, bxl)
    call destroy(bxl)
  end subroutine boxarray_simplify

  subroutine boxarray_print(ba, str, unit, legacy, skip)
    use bl_IO_module
    type(boxarray), intent(in) :: ba
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: legacy
    integer, intent(in), optional :: skip
    integer :: i
    integer :: un
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt = '("BOXARRAY[(*")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(" ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i5)') ba%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NBOXES  = ",i5)') ba%nboxes
    call unit_skip(un, skip)
    write(unit=un, fmt='(" *) {")')
    do i = 1, ba%nboxes
       call print(ba%bxs(i), unit=unit, advance = 'NO', &
            legacy = legacy, skip = unit_get_skip(skip)+ 1)
       if ( i == ba%nboxes ) then
          write(unit=un, fmt='("}]")')
       else
          write(unit=un, fmt='(",")')
       end if
    end do
  end subroutine boxarray_print

  function boxarray_clean(boxes) result(r)
    logical :: r
    type(box), intent(in), dimension(:) :: boxes
    integer :: i, j

    do i = 1, size(boxes)-1
       do j = i+1, size(boxes)
          if ( intersects(boxes(i),boxes(j)) ) then
             r = .FALSE.
             return
          end if
       end do
    end do
    r = .TRUE.

  end function boxarray_clean

  subroutine boxarray_add_clean(ba, bx)

    use bl_prof_module

    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bx
    type(list_box) :: check, tmp, tmpbl, bl
    type(list_box_node), pointer :: cp, lp
    type(bl_prof_timer), save :: bpt

    if ( empty(ba) ) then
       call boxarray_build_bx(ba, bx)
       return
    end if
    call build(bpt, "ba_add_clean")
    call list_build_v_box(bl, ba%bxs)
    call push_back(check, bx)
    lp => begin(bl)
    do while ( associated(lp) )
       cp => begin(check)
       do while ( associated(cp) )
          if ( intersects(value(cp), value(lp)) ) then
             tmpbl = boxlist_box_diff(value(cp), value(lp))
             call splice(tmp, tmpbl)
             cp => erase(check, cp)
          else
             cp => next(cp)
          end if
       end do
       call splice(check, tmp)
       lp => next(lp)
    end do
    call splice(bl, check)
    call boxlist_simplify(bl)
    call boxarray_build_copy_l(ba, bl)
    call list_destroy_box(bl)
    call destroy(bpt)

  end subroutine boxarray_add_clean

  subroutine boxarray_add_clean_boxes(ba, bxs, simplify)

    use bl_prof_module

    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bxs(:)
    logical, intent(in), optional :: simplify
    logical :: lsimplify
    type(list_box) :: check, tmp, tmpbl, bl
    type(list_box_node), pointer :: cp, lp
    integer :: i
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ba_add_clean_boxes")

    lsimplify = .true.; if ( present(simplify) ) lsimplify = simplify

    if ( size(bxs) .eq. 0 ) return
    
    if ( empty(ba) ) call boxarray_build_bx(ba, bxs(1))

    call list_build_v_box(bl, ba%bxs)

    do i = 1, size(bxs)
       call list_build_v_box(check, bxs(i:i))
       lp => begin(bl)
       do while ( associated(lp) )
          cp => begin(check)
          do while ( associated(cp) )
             if ( intersects(value(cp), value(lp)) ) then
                tmpbl = boxlist_box_diff(value(cp), value(lp))
                call splice(tmp, tmpbl)
                cp => erase(check, cp)
             else
                cp => next(cp)
             end if
          end do
          call splice(check, tmp)
          lp => next(lp)
       end do
       call splice(bl, check)
    end do

    if ( lsimplify ) call boxlist_simplify(bl)

    call boxarray_build_copy_l(ba, bl)
    call list_destroy_box(bl)
    call destroy(bpt)

  end subroutine boxarray_add_clean_boxes

  subroutine boxarray_to_domain(ba)
    type(boxarray), intent(inout) :: ba
    type(boxarray) :: ba1
    call boxarray_add_clean_boxes(ba1, ba%bxs)
    call boxarray_destroy(ba)
    ba = ba1
  end subroutine boxarray_to_domain

  function boxarray_box_contains(ba, bx) result(r)
    use bl_error_module
    logical                    :: r
    type(boxarray), intent(in) :: ba
    type(box),      intent(in) :: bx

    type(list_box) :: bl1, bl
    type(box)      :: bx1
    integer        :: i

    if ( nboxes(ba) .eq. 0 ) &
       call bl_error('Empty boxarray in boxarray_box_contains')

    call build(bl1)
    do i = 1, nboxes(ba)
       bx1 = intersection(bx, get_box(ba,i))
       if ( empty(bx1) ) cycle
       call push_back(bl1, bx1)
    end do
    bl = boxlist_boxlist_diff(bx, bl1)
    r = empty(bl)
    call destroy(bl)
    call destroy(bl1)

  end function boxarray_box_contains

  function boxarray_boxarray_contains(ba1, ba2, allow_empty) result(r)
    use bl_error_module
    logical :: r
    type(boxarray), intent(in) :: ba1, ba2
    logical, intent(in), optional :: allow_empty

    integer :: i
    logical :: lallow

    !Note that allow_empty refers to permitting empty boxes, not empty boxarrays
    lallow = .false.; if (present(allow_empty)) lallow=allow_empty

    if ( nboxes(ba1) .eq. 0 ) &
       call bl_error('Empty boxarray ba1 in boxarray_boxarray_contains')

    if ( nboxes(ba2) .eq. 0 ) &
       call bl_error('Empty boxarray ba2 in boxarray_boxarray_contains')

    if ( lallow) then
       do i = 1, nboxes(ba2)
          if (empty(get_box(ba2,i))) cycle !ignore empty boxes
          r = boxarray_box_contains(ba1, get_box(ba2,i)) 
          if ( .not. r ) return
       end do
    else
       do i = 1, nboxes(ba2)
          r = boxarray_box_contains(ba1, get_box(ba2,i)) 
          if ( .not. r ) return
       end do
    endif

    r = .true.

  end function boxarray_boxarray_contains

end module boxarray_module
