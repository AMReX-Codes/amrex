!! A _BoxArray_ is an array of boxes.
module boxarray_module

  use bl_types
  use box_module
  use list_box_module
  use bl_mem_stat_module

  implicit none

  type boxarray
     integer :: dim = 0
     integer :: nboxes = 0
     type(box), pointer :: bxs(:) => Null()
  end type boxarray

  interface dim
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
  end interface

  interface build
     module procedure boxarray_build_v
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
     module procedure boxarray_intersection_ba
     module procedure boxarray_intersection_bx
  end interface

  interface boxarray_grow
!     module procedure boxarray_grow_m
     module procedure boxarray_grow_n
     module procedure boxarray_grow_n_f
     module procedure boxarray_grow_n_d_f
     module procedure boxarray_grow_v
     module procedure boxarray_grow_v_f
  end interface

  interface boxarray_boundary
!    module procedure boxarray_boundary_m
!    module procedure boxarray_boundary_n
!    module procedure boxarray_boundary_n_f
!    module procedure boxarray_boundary_n_d_f
!    module procedure boxarray_boundary_v
!    module procedure boxarray_boundary_v_f
  end interface

  interface boxarray_box_boundary
!    module procedure boxarray_box_boundary_m
!    module procedure boxarray_box_boundary_n
!    module procedure boxarray_box_boundary_n_f
!    module procedure boxarray_box_boundary_n_d_f
!    module procedure boxarray_box_boundary_v
!    module procedure boxarray_box_boundary_v_f
  end interface

  interface bbox
     module procedure boxarray_bbox
  end interface

  interface print
     module procedure boxarray_print
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

  private :: boxlist_simplify
  private :: boxarray_maxsize_l
  private :: boxlist_box_diff
  private :: boxlist_boxlist_diff
  private :: boxlist_build_a
  private :: boxlist_nboxes
  private :: boxlist_verify_dim
  private :: boxarray_build_l
  private :: boxarray_build_copy_l

  type(mem_stats), private, save :: boxarray_ms

contains
  
  function boxarray_equal(ba1, ba2) result(r)
    type(boxarray), intent(in) :: ba1, ba2
    logical :: r
    r = associated(ba1%bxs, ba2%bxs)
  end function boxarray_equal

  function boxarray_not_equal(ba1, ba2) result(r)
    type(boxarray), intent(in) :: ba1, ba2
    logical :: r
    r = .not. associated(ba1%bxs, ba2%bxs)
  end function boxarray_not_equal

  function boxarray_same_q(ba1, ba2) result(r)
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
            exit
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

  function boxarray_empty(ba) result(r)
    logical :: r
    type(boxarray), intent(in) :: ba
    r = ba%nboxes == 0
  end function boxarray_empty

  function boxarray_built_q(ba) result(r)
    logical :: r
    type(boxarray), intent(in) :: ba
    r = ba%dim /= 0
  end function boxarray_built_q

  function boxarray_dim(ba) result(r)
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

  function boxarray_get_box(ba, i) result(r)
    type(boxarray), intent(in) :: ba
    integer, intent(in) :: i
    type(box) :: r
    r = ba%bxs(i)
  end function boxarray_get_box

  subroutine boxarray_build_copy(ba, ba1)
    type(boxarray), intent(inout) :: ba
    type(boxarray), intent(in) :: ba1
    
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_COPY: already built")
    ! If ba1 is not built then return
    if ( .not. built_q(ba1) ) return
    ba%nboxes = size(ba1%bxs)
    allocate(ba%bxs(size(ba1%bxs)))
    ba%bxs = ba1%bxs
    ba%dim = ba1%dim
    call boxarray_verify_dim(ba)
    ! call boxarray_sort(ba) ! shouldn't be needed, ba1 same sort as ba
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_copy

  subroutine boxarray_build_copy_l(ba, bl)
    type(boxarray), intent(inout) :: ba
    type(list_box), intent(in) :: bl
    
    if ( built_q(ba) ) call destroy(ba)
    call boxarray_build_l(ba, bl)
    ! call boxarray_sort(ba) ! boxarray_build_l already does it
    ! call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_copy_l

  subroutine boxarray_build_v(ba, bxs, sort)
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
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bx
    
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_BX: already built")
    ba%nboxes = 1
    allocate(ba%bxs(1))
    ba%bxs(1) = bx
    ba%dim = bx%dim
    call boxarray_verify_dim(ba)
    ! call boxarray_sort(ba)  !! already sorted, one box
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_bx

  subroutine boxarray_build_l(ba, bl)
    type(boxarray), intent(inout) :: ba
    type(list_box), intent(in) :: bl
    type(list_box_node), pointer :: bln
    integer :: i

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
    call boxarray_sort(ba)
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_l

  subroutine boxarray_destroy(ba)
    type(boxarray), intent(inout) :: ba
!    if ( associated(ba%bxs) ) then
       call mem_stats_dealloc(boxarray_ms, ba%nboxes)
       deallocate(ba%bxs)
       ba%dim = 0
       ba%nboxes = 0
!    end if
  end subroutine boxarray_destroy

  subroutine boxlist_build_a(bl, ba)
    type(boxarray), intent(in) :: ba
    type(list_box), intent(out) :: bl
    integer :: i
    do i = 1, ba%nboxes
       call push_back(bl, ba%bxs(i))
    end do
  end subroutine boxlist_build_a

  subroutine boxarray_sort(ba)
    use sort_box_module
    type(boxarray), intent(inout) :: ba
    call box_sort(ba%bxs)
  end subroutine boxarray_sort

  subroutine boxarray_verify_dim(ba, stat)
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

  subroutine boxlist_verify_dim(bl, stat)
    type(list_box), intent(in) :: bl
    integer, intent(out), optional :: stat
    type(list_box_node), pointer :: bln
    type(box) :: bx
    integer :: dm
    if ( present(stat) ) stat = 0
    if ( size(bl) < 1 ) return
    bln => begin(bl)
    bx = value(bln)
    dm = bx%dim
    do while (associated(bln))
       bx = value(bln)
       if ( bx%dim /= dm ) then
          if ( present(stat) ) then
             stat = 1
             return
          else
             call bl_error("BOXLIST_VERIFY_DIM:" // &
                  "some box's dim not equal to the first box's dim: ", dm)
          end if
       end if
    end do
  end subroutine boxlist_verify_dim

!   subroutine boxarray_grow_m(ba, mat)
!     type(boxarray), intent(inout) :: ba
!     integer, intent(in) :: mat(:,:)
!     integer :: i
!     do i = 1, ba%nboxes
!        ba%bxs(i) = grow(ba%bxs(i), mat)
!     end do
!   end subroutine boxarray_grow_m
  subroutine boxarray_grow_v(ba, rv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:)
    integer i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), rv)
    end do
  end subroutine boxarray_grow_v
  subroutine boxarray_grow_v_f(ba, rv, face)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:), face
    integer i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), rv, face)
    end do
  end subroutine boxarray_grow_v_f
  subroutine boxarray_grow_n(ba, n)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: n
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), n)
    end do
  end subroutine boxarray_grow_n
  subroutine boxarray_grow_n_f(ba, n, face)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: n, face
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), n, face)
    end do
  end subroutine boxarray_grow_n_f
  subroutine boxarray_grow_n_d_f(ba, n, dim, face)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: n, face, dim
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), n, dim, face)
    end do
  end subroutine boxarray_grow_n_d_f

  subroutine boxarray_nodalize(ba, nodal)
    type(boxarray), intent(inout) :: ba
    logical, intent(in), optional :: nodal(:)
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = box_nodalize(ba%bxs(i), nodal)
    end do
  end subroutine boxarray_nodalize

  subroutine boxarray_coarsen_v_m(ba, cv, mask)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: cv(:)
    logical, intent(in) :: mask(:)
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), cv, mask)
    end do
  end subroutine boxarray_coarsen_v_m
  subroutine boxarray_coarsen_v(ba, cv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: cv(:)
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), cv)
    end do
  end subroutine boxarray_coarsen_v
  subroutine boxarray_coarsen_i_m(ba, ci, mask)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ci
    logical, intent(in) :: mask(:)
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), ci, mask)
    end do
  end subroutine boxarray_coarsen_i_m
  subroutine boxarray_coarsen_i(ba, ci)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ci
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), ci)
    end do
  end subroutine boxarray_coarsen_i

  subroutine boxarray_shift_v(ba, rv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:)
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = shift(ba%bxs(i), rv)
    end do
  end subroutine boxarray_shift_v
  subroutine boxarray_shift_i(ba, ri)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ri
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = shift(ba%bxs(i), ri)
    end do
  end subroutine boxarray_shift_i
    
  subroutine boxarray_refine_v(ba, rv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:)
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = refine(ba%bxs(i), rv)
    end do
  end subroutine boxarray_refine_v
  subroutine boxarray_refine_i(ba, ri)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ri
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = refine(ba%bxs(i), ri)
    end do
  end subroutine boxarray_refine_i

  subroutine boxarray_intersection_ba(bai, ba)
    type(boxarray), intent(inout) :: bai
    type(boxarray), intent(in)    :: ba
    integer :: i, j
    do i = 1, ba%nboxes
       do j = 1, bai%nboxes
          bai%bxs(j) = intersection(bai%bxs(j), ba%bxs(i))
       end do
    end do
    call boxarray_simplify(bai)
  end subroutine boxarray_intersection_ba

  subroutine boxarray_intersection_bx(ba, bx)
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bx
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = intersection(ba%bxs(i), bx)
    end do
    call boxarray_simplify(ba)
  end subroutine boxarray_intersection_bx
    
!   subroutine boxarray_box_boundary_m(bao, bx, mat)
!     type(boxarray), intent(out) :: bao
!     type(box), intent(in)  :: bx
!     integer, intent(in) :: mat(:,:)
!     type(boxarray) :: baa
!     call boxarray_build_bx(baa, bx)
!     call boxarray_boundary(bao, baa, mat)
!     call boxarray_destroy(baa)
!   end subroutine boxarray_box_boundary_m
  subroutine boxarray_box_boundary_v_f(bao, bx, nv, face)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: nv(:), face
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_v_f(bao, baa, nv, face)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_v_f
  subroutine boxarray_box_boundary_v(bao, bx, nv)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: nv(:)
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_v(bao, baa, nv)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_v
  subroutine boxarray_box_boundary_n_f(bao, bx, n, face)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: n, face
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_n_f(bao, baa, n, face)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_n_f
  subroutine boxarray_box_boundary_n_d_f(bao, bx, n, dim, face)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: n, face, dim
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_n_d_f(bao, baa, n, dim, face)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_n_d_f
  subroutine boxarray_box_boundary_n(bao, bx, n)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: n
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_n(bao, baa, n)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_n

!   subroutine boxarray_boundary_m(bao, ba, mat)
!     type(boxarray), intent(out) :: bao
!     type(boxarray), intent(in)  :: ba
!     integer, intent(in) :: mat(:,:)
!     call boxarray_build_copy(bao, ba)
!     call boxarray_grow(bao, mat)
!     call boxarray_diff(bao, ba)
!   end subroutine boxarray_boundary_m
  subroutine boxarray_boundary_v_f(bao, ba, nv, face)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer, intent(in) :: nv(:), face
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, nv, face)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_v_f
  subroutine boxarray_boundary_v(bao, ba, nv)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer, intent(in) :: nv(:)
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, nv)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_v
  subroutine boxarray_boundary_n_d_f(bao, ba, n, dim, face)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer, intent(in) :: n, face, dim
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, n, dim, face)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_n_d_f
  subroutine boxarray_boundary_n_f(bao, ba, n, face)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer, intent(in) :: n, face
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, n, face)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_n_f
  subroutine boxarray_boundary_n(bao, ba, n)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer, intent(in) :: n
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, n)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_n

  function boxarray_nboxes(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer :: r
    r = ba%nboxes
  end function boxarray_nboxes

  function boxlist_nboxes(bl) result(r)
    type(list_box), intent(in) :: bl
    integer :: r
    r = size(bl)
  end function boxlist_nboxes

  function boxarray_volume(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer(kind=ll_t) :: r
    integer :: i
    r = 0_ll_t
    do i = 1, ba%nboxes
       r = r + box_volume(ba%bxs(i))
    end do
  end function boxarray_volume

  function boxarray_dvolume(ba) result(r)
    type(boxarray), intent(in) :: ba
    real(dp_t) :: r
    integer :: i
    r = 0
    do i = 1, ba%nboxes
       r = r + box_dvolume(ba%bxs(i))
    end do
  end function boxarray_dvolume

  function boxarray_bbox(ba) result(r)
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

  subroutine boxarray_boxarray_diff(ba, bx, ba1)
    type(boxarray), intent(out) :: ba
    type(boxarray), intent(in) :: ba1
    type(box) :: bx
    type(list_box) :: bl1, bl
    call build(bl1, ba1%bxs)
    bl = boxlist_boxlist_diff(bx, bl1)
    call boxarray_build_l(ba, bl)
    call destroy(bl)
    call destroy(bl1)
  end subroutine boxarray_boxarray_diff

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
    integer nl, bs, rt, nblk, sz, ex, ks, ps
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

  function boxlist_box_diff(bx1, b2) result(r)
    type(box), intent(in) :: bx1, b2
    type(list_box) :: r
    type(box) :: b1, bn
    integer, dimension(bx1%dim) :: b2lo, b2hi, b1lo, b1hi
    integer i, dm

    dm = bx1%dim
    b1 = bx1
    if ( .not. contains(b2,b1) ) then
       if ( .not. intersects(b1,b2) ) then
          call push_back(r, b1)
       else
          b2lo = lwb(b2); b2hi = upb(b2)
          do i = 1, dm
             b1lo = lwb(b1); b1hi = upb(b1)
             if ( b1lo(i) < b2lo(i) .AND. b2lo(i) <= b1hi(i) ) then
                bn = b1
                call set_lwb(bn, i, b1lo(i))
                call set_upb(bn, i, b2lo(i)-1)
                call push_back(r, bn)
                call set_lwb(b1, i, b2lo(i))
             end if
             if ( b1lo(i) <= b2hi(i) .AND. b2hi(i) < b1hi(i) ) then
                bn = b1
                call set_lwb(bn, i, b2hi(i)+1)
                call set_upb(bn, i, b1hi(i))
                call push_back(r, bn)
                call set_upb(b1, i, b2hi(i))
             end if
          end do
       end if
    end if

  end function boxlist_box_diff

  ! r = bx - bxl
  function boxlist_boxlist_diff(bx, bxl) result(r)
    type(box), intent(in) :: bx
    type(list_box), intent(in) :: bxl
    type(list_box) :: r, bl
    type(list_box_node), pointer :: blp, bp

    call push_back(r, bx)
    blp => begin(bxl)
    do while ( associated(blp) .AND. .NOT. empty(r) )
       bp => begin(r)
       do while ( associated(bp) )
          if ( intersects(value(bp), value(blp)) ) then
             bl = boxlist_box_diff(value(bp), value(blp))
             call splice(r, bl)
             bp => erase(r, bp)
          else
             bp => next(bp)
          end if
       end do
       blp => next(blp)
    end do

  end function boxlist_boxlist_diff

  subroutine boxarray_simplify(bxa)
    type(boxarray), intent(inout) :: bxa
    type(list_box) :: bxl
    integer :: i
    do i = 1, bxa%nboxes
       call push_back(bxl, bxa%bxs(i))
    end do
    call boxlist_simplify(bxl)
    call boxarray_destroy(bxa)
    call boxarray_build_l(bxa, bxl)
    call destroy(bxl)
  end subroutine boxarray_simplify

  !! Takes a boxlist and normalizes it by attaching boxes
  !! that abut across periodic edges
  subroutine boxlist_normalize(bxl, pd, pmask)
    type(list_box), intent(inout) :: bxl
    type(box), intent(in) :: pd
    logical, intent(in) :: pmask(:)
    integer :: dm
    integer :: ext(pd%dim)
    if ( size(bxl) == 0 ) return
    dm = dim(front(bxl))
    ext = extent(pd)
    call norm
  contains
    subroutine norm
      type(list_box_node), pointer :: ba, bb
      integer :: i
      type(box) :: lbx
      integer :: lo(dm), hi(dm)
      ba => begin(bxl)
      do while ( associated(ba) )
         bb => next(ba)
         do while ( associated(bb) )
            do i = 1, dm
               lbx = shift(value(bb),i,-ext(i))
               if ( box_abut(value(ba), lbx, i) ) then
                  lo = lwb(value(ba))
                  hi = upb(value(ba))
                  lo(i) = min(lwb(value(ba),i),lwb(lbx,i))
                  hi(i) = max(upb(value(ba),i),upb(lbx,i))
                  print *, 'got one'
               end if
               lbx = shift(value(bb),i,+ext(i))
            end do
         end do
      end do
    end subroutine norm

  end subroutine boxlist_normalize

  subroutine boxlist_simplify(bxl)
    type(list_box), intent(inout) :: bxl
    integer :: dm
    
    if ( size(bxl) == 0 ) return
    dm = dim(front(bxl))
    do while ( simp_old() > 0 )
    end do

  contains

    function simp() result(cnt)
      integer :: cnt
      integer :: joincnt
      integer, dimension(dm) :: lo, hi, alo, ahi, blo, bhi
      type(list_box_node), pointer :: ba, bb
      type(box) :: bx
      logical match, canjoin, t
      integer :: i

      if ( size(bxl) == 0 ) return
      cnt = 0
      ba => begin(bxl)

      do while ( associated(ba) )
         match = .FALSE.
         bb => next(ba)
         do while ( associated(bb) )
            do i = 1, dm
               t = box_merge(bx, value(ba), value(bb), i)
               if ( t ) exit
            end do
            if ( t ) then
               call set(bb, bx)
               ba => erase(bxl, ba)
               cnt = cnt + 1
               match = .TRUE.
               exit
            else
               ! No match found, try next element.
               bb => next(bb)
            end if
         end do
         ! If a match was found, a was already advanced in the list.
         if (.not. match) then
            ba => next(ba)
         end if
      end do

    end function simp

    function simp_old() result(cnt)
      integer :: cnt
      integer :: joincnt
      integer, dimension(dm) :: lo, hi, alo, ahi, blo, bhi
      type(list_box_node), pointer :: ba, bb
      type(box) :: bx
      logical match, canjoin
      integer :: i

      if ( size(bxl) == 0 ) return
      cnt = 0
      ba => begin(bxl)

      do while ( associated(ba) )
         alo = lwb(value(ba)); ahi = upb(value(ba))
         match = .FALSE.
         bb => next(ba)
         do while ( associated(bb) )
            blo = lwb(value(bb)); bhi = upb(value(bb))
            canjoin = .TRUE.
            joincnt = 0
            do i = 1, dm
               if (alo(i) == blo(i) .AND. ahi(i)==bhi(i)) then
                  lo(i) = alo(i)
                  hi(i) = ahi(i)
               else if (alo(i)<=blo(i) .AND. blo(i)<=ahi(i)+1) then
                  lo(i) = alo(i)
                  hi(i) = max(ahi(i),bhi(i))
                  joincnt = joincnt + 1
               else if (blo(i)<=alo(i) .AND. alo(i)<=bhi(i)+1) then
                  lo(i) = blo(i)
                  hi(i) = max(ahi(i),bhi(i))
                  joincnt = joincnt + 1
               else
                  canjoin = .FALSE.
                  exit
               end if
            end do
            if (canjoin .AND. (joincnt <= 1)) then
               ! Modify b and remove a from the list.
               call build(bx, lo, hi)
               call set(bb, bx)
               ba => erase(bxl, ba)
               cnt = cnt + 1
               match = .TRUE.
               exit
            else
               ! No match found, try next element.
               bb => next(bb)
            end if
         end do
         ! If a match was found, a was already advanced in the list.
         if (.not. match) then
            ba => next(ba)
         end if
      end do

    end function simp_old
  end subroutine boxlist_simplify

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
    write(unit=un, fmt = '("BOXARRAY")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(": ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i5)') ba%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NBOXES  = ",i5)') ba%nboxes
    do i = 1, ba%nboxes
       call print(ba%bxs(i), unit=unit, advance = 'yes', &
            legacy = legacy, skip = unit_get_skip(skip)+ 1)
    end do
  end subroutine boxarray_print

  subroutine boxlist_print(bl, str, unit, legacy)
    use bl_IO_module
    type(list_box), intent(in) :: bl
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: legacy
    type(list_box_node), pointer :: bn
    integer :: un
    un = unit_stdout(unit)
    if ( present(str) ) then
       write(unit=un, fmt='("(*",A,"*)")') str
    end if
    write(unit=un, fmt = '("boxlist[ ", i10, ", {")') size(bl)
    bn => begin(bl)
    do while ( associated(bn) )
       call print(value(bn), unit=unit)
       bn => next(bn)
    end do
  end subroutine boxlist_print

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
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bx
    type(list_box) :: check, tmp, tmpbl, bl
    type(list_box_node), pointer :: cp, lp

    if ( empty(ba) ) then
       call boxarray_build_bx(ba, bx)
       return
    end if
    call build(bl, ba%bxs)
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
    call destroy(bl)

  end subroutine boxarray_add_clean

  subroutine boxarray_add_clean_boxes(ba, bxs)
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bxs(:)
    type(list_box) :: check, tmp, tmpbl, bl
    type(list_box_node), pointer :: cp, lp
    integer :: i
    if ( empty(ba) ) then
       call boxarray_build_bx(ba, bxs(1))
    end if
    call build(bl, ba%bxs)
    do i = 1, size(bxs)
       call build(check, bxs(i:i))
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
    call boxlist_simplify(bl)
    call boxarray_build_copy_l(ba, bl)
    call destroy(bl)

  end subroutine boxarray_add_clean_boxes

  subroutine boxarray_to_domain_1(ba)
    type(boxarray), intent(inout) :: ba
    type(boxarray) :: ba1
    integer :: i
    do i = 1, ba%nboxes
       call boxarray_add_clean(ba1, ba%bxs(i))
    end do
    call boxarray_destroy(ba)
    ba = ba1
  end subroutine boxarray_to_domain_1

  subroutine boxarray_to_domain(ba)
    type(boxarray), intent(inout) :: ba
    type(boxarray) :: ba1
    call boxarray_add_clean_boxes(ba1, ba%bxs)
    call boxarray_destroy(ba)
    ba = ba1
  end subroutine boxarray_to_domain

  subroutine boxarray_pn_domain_ba(bap, ba, pd, nproper, rr)
    type(boxarray), intent(out) :: bap
    type(boxarray), intent(in)  :: ba
    type(boxarray), intent(in)  :: pd
    integer, intent(in)         :: nproper
    integer, intent(in), optional :: rr
    integer :: rrv(pd%dim)
    if ( present(rr) ) then
       rrv = rr
       call boxarray_pn_domain_ba_v(bap, ba, pd, nproper, rrv)
    else
       call boxarray_pn_domain_ba_v(bap, ba, pd, nproper)
    end if
  end subroutine boxarray_pn_domain_ba

  subroutine boxarray_pn_domain_bx(bap, ba, pd, nproper, rr)
    type(boxarray), intent(out) :: bap
    type(boxarray), intent(in)  :: ba
    type(box), intent(in)  :: pd
    integer, intent(in)         :: nproper
    type(boxarray) :: pda
    integer, intent(in), optional :: rr
    call build(pda, pd)
    call boxarray_pn_domain_ba(bap, ba, pda, nproper, rr)
    call destroy(pda)
  end subroutine boxarray_pn_domain_bx

  !! This is the one that actually does things.
  subroutine boxarray_pn_domain_ba_v(bap, ba, pd, nproper, rr)
    type(boxarray), intent(out) :: bap
    type(boxarray), intent(in)  :: ba
    type(boxarray), intent(in)  :: pd
    integer, intent(in)         :: nproper
    type(boxarray) :: cdl, dl
    integer, intent(in), optional :: rr(:)
    call copy(dl, ba)
    if ( present(rr) ) then
       call boxarray_refine(dl, rr)
    end if
    call copy(cdl, pd); call boxarray_diff(cdl,  dl)
    call boxarray_grow(cdl, nproper)
    call boxarray_intersection(cdl, pd)
    call copy(bap, pd); call boxarray_diff(bap, cdl)
    call destroy(cdl)
    call destroy(dl)
  end subroutine boxarray_pn_domain_ba_v

  subroutine boxarray_pn_domain_bx_v(bap, ba, pd, nproper, rr)
    type(boxarray), intent(out) :: bap
    type(boxarray), intent(in)  :: ba
    type(box), intent(in)  :: pd
    integer, intent(in)         :: nproper
    type(boxarray) :: pda
    integer, intent(in), optional :: rr(:)
    call build(pda, pd)
    call boxarray_pn_domain_ba_v(bap, ba, pda, nproper, rr)
    call destroy(pda)
  end subroutine boxarray_pn_domain_bx_v

  function boxlist_decompose_bx_bx(bx1, bx2) result(r)
    type(list_box) :: r
    type(box), intent(in) :: bx1, bx2
    type(box) :: tb(27)
    integer :: n, i
    call box_decompose(tb, n, bx1, bx2)
    do i = 1, n
       call push_back(r, tb(i))
    end do
  end function boxlist_decompose_bx_bx

  function boxlist_decompose_bl_bx(bl, bx) result(r)
    type(list_box) :: r
    type(list_box), intent(in) :: bl
    type(list_box) :: tbl
    type(box), intent(in) :: bx
    type(list_box_node), pointer :: bn
    bn => begin(bl)
    do while ( associated(bn) )
       tbl = boxlist_decompose_bx_bx(value(bn), bx)
       call splice(r, tbl)
       bn => next(bn)
    end do
  end function boxlist_decompose_bl_bx

  function boxlist_decompose_bl_bl(bl1, bl2) result(r)
    type(list_box) :: r
    type(list_box), intent(in) :: bl1, bl2
    type(list_box) :: tbl
    type(list_box_node), pointer :: bn
    if ( empty(bl2) ) then
       call list_copy_box(r, bl1)
       return
    end if
    bn => begin(bl2)
    r = boxlist_decompose_bl_bx(bl1, value(bn))
    bn => next(bn)
    do while ( associated(bn) )
       tbl = boxlist_decompose_bl_bx(r, value(bn))
       call destroy(r)
       r = tbl
       bn => next(bn)
    end do
  end function boxlist_decompose_bl_bl

  subroutine boxarray_decompose(bad, ba)
    type(boxarray), intent(out) :: bad
    type(boxarray), intent(in)  :: ba
    type(list_box):: bl, tbl

    call boxlist_build_a(bl, ba)
    tbl = boxlist_decompose_bl_bl(bl, bl)
    call boxarray_build_copy_l(bad, tbl)
    call destroy(tbl)
    call destroy(bl)
  end subroutine boxarray_decompose

end module boxarray_module
