!! Support for box calculus
!!
!! _Boxs_ are define rectangular domains in space.  They can also be
!! thought of as providing a definition of a Cell-Centered index region.
module box_module

  use bl_types
  use bl_space
  use bl_error_module

  implicit none

  !! Boxes consist of
  !!   - A dimensionality, and
  !!   - A low and high end of cells
  !!
  !! Initially, each box is completely empty, with a low end at positive
  !! infinity and a low end at negative infinity.  Almost all box operations
  !! are ill-defined on default boxes.
  type box
     integer :: dim  = 0
     integer :: lo(MAX_SPACEDIM) = Huge(1)
     integer :: hi(MAX_SPACEDIM) = -Huge(1)
  end type box

  integer :: box_print_fixed_width = 5

  !! Builds a box
  interface build
     module procedure box_build_2
     module procedure box_build_1
  end interface

  !! Returns a default box
  interface nobox
     module procedure box_nobox
  end interface

  !! Returns a box that covers all index space [-Inf,+Inf]
  interface allbox
     module procedure box_allbox
  end interface

  !! Returns the dimensionality of the box.
!  interface dim
!     module procedure box_dim
!  end interface

  !! Returns the volume of the box as a ll_t
  interface volume
     module procedure box_volume
  end interface

  !! Returns the volume of the box as REAL(kind=dp_t)
  interface dvolume
     module procedure box_dvolume
  end interface

  !! Returns the upper bound of a box, or of a box in
  !! a given direction
  interface upb
     module procedure box_upb
     module procedure box_upb_d
  end interface
  !! Sets the upper bound of a box, or sets the upper
  !! bound of a box in a given direction
  interface set_upb
     module procedure box_set_upb
     module procedure box_set_upb_d
  end interface set_upb

  !! Returns the lower bound of a box, or of a box in
  !! a given direction
  interface lwb
     module procedure box_lwb
     module procedure box_lwb_d
  end interface
  !! Sets the lower bound of a box, or sets the lower
  !! bound of a box in a given direction
  interface set_lwb
     module procedure box_set_lwb
     module procedure box_set_lwb_d
  end interface set_lwb

  !! returns the refinement of a box
  interface refine
     module procedure box_refine_v_m
     module procedure box_refine_v
     module procedure box_refine_i_m
     module procedure box_refine_i
  end interface

  !! returns a coarsened version of a box
  interface coarsen
     module procedure box_coarsen_v_m
     module procedure box_coarsen_v
     module procedure box_coarsen_i_m
     module procedure box_coarsen_i
  end interface

  !! grows a box
  interface grow
!    module procedure box_grow_m
     module procedure box_grow_n
     module procedure box_grow_n_f
     module procedure box_grow_n_d_f
     module procedure box_grow_v
     module procedure box_grow_v_f
  end interface

  !! shifts a box
  interface shift
     module procedure box_shift_d
     module procedure box_shift_i
     module procedure box_shift_v
  end interface shift

  !! Returns true if two boxes are equal
  interface equal
     module procedure box_equal
  end interface
  interface operator(.EQ.)
     module procedure box_equal
  end interface

  !! Returns true if two boxes are not equal
  interface not_equal
     module procedure box_not_equal
  end interface
  interface operator(.NE.)
     module procedure box_not_equal
  end interface

  !! Returns true if the lwb(bx1) < lwb(bx2)
  interface operator(.LT.)
     module procedure box_less
  end interface

  !! Returns the intersection of two boxes
  interface intersection
     module procedure box_intersection
  end interface
  interface operator(.INTERSECT.)
     module procedure box_intersection
  end interface

  !! Returns true if two boxes intersect
  interface intersects
     module procedure box_intersects
  end interface

  !! returns the bounding box of two boxes; that is
  !! the smallest box the contains them both.
  interface bbox
     module procedure box_bbox
  end interface

  !! Returns the size, or extent, of a box or of a box 
  !! in a given direction.
  interface extent
     module procedure box_extent
     module procedure box_extent_d
  end interface

  !! Returns the boundary of a box as a matrix of boxes.
!  interface boundary
!     module procedure box_boundary_n
!    module procedure box_boundary_m
!  end interface

  !! Prints a box.
  interface print
     module procedure box_print
  end interface

  interface contains
     module procedure box_contains
     module procedure box_contains_iv
  end interface

  interface empty
     module procedure box_empty
  end interface

  interface reduce
     module procedure box_reduce
  end interface

contains

  elemental function int_coarsen(v, i) result(r)
    integer, intent(in) :: v
    integer, intent(in) :: i
    integer :: r
    if ( v < 0 ) then
       r = -abs(v+1)/i - 1
    else
       r = v/i
    end if
  end function int_coarsen

  !! Makes a box given a low and high end.
  !! If given a single point (index), it returs a unit box located
  !! at the that point.
  function make_box(lo, hi) result(r)
    integer, intent(in), dimension(:) :: lo
    integer, intent(in), dimension(:), optional :: hi
    type(box) :: r
    if ( present(hi) ) then
       call box_build_2(r, lo, hi)
    else
       call box_build_2(r, lo, lo)
    end if
  end function make_box

  !! Returns a box defined on the entire index domain [-Inf,+Inf] of
  !! dimension _dim_.
  function box_allbox(dim) result(r)
    integer, intent(in) :: dim
    type(box) :: r
    r%dim = dim
    r%lo(1:dim) = -Huge(1)+1
    r%hi(1:dim) =  Huge(1)-1
  end function box_allbox

  !! Returns a completely empty box [+Inf,-Inf]; this box is also
  !! the default box of dimension _dim_.
  function box_nobox(dim) result(r)
    integer, intent(in) :: dim
    type(box) :: r
    r%dim = dim
    r%lo(1:dim) =  Huge(1)
    r%hi(1:dim) = -Huge(1)
  end function box_nobox

  !! Returns a unit_box at the origin of index space of dimension _dim_.
  function unit_box(dim) result(r)
    type(box) :: r
    integer, intent(in) :: dim
    r%dim = dim
    r%lo(1:dim) = 0
    r%hi(1:dim) = 0
  end function unit_box

  !! Returns .TRUE. if the box _bx_ contains no cells, or if it is of dimension
  !! zero.
  function box_empty(bx) result(r)
    logical :: r
    type(box), intent(in) :: bx
    r = ( bx%dim == 0) .or. ( any(bx%lo(1:bx%dim) > bx%hi(1:bx%dim)) )
  end function box_empty

  !! Builds a box given a low and high end of index space.  It is an error
  !! for the sizes of _lo_ and _hi_ to be unequal or for the size of _lo_
  !! to exceed MAX_SPACEDIM.
  subroutine box_build_2(bx, lo, hi)
    type(box), intent(out) :: bx
    integer, intent(in) :: lo(:), hi(:)
    if ( size(lo) /= size(hi) ) then
       call bl_error("BOX_BUILD_2: lo /= hi")
    end if
    if ( size(lo) > MAX_SPACEDIM ) then
       call bl_error("BOX_BUILD_2: size(lo,hi) > MAX_SPACEDIM")
    end if
    bx%dim = size(lo)
    bx%lo(1:bx%dim) = lo
    bx%hi(1:bx%dim) = hi
  end subroutine box_build_2

  !! Builds a unit box at position _lo_.  It is an error for size of _lo_
  !! to exceed MAX_SPACEDIM.
  subroutine box_build_1(bx, lo)
    type(box), intent(out) :: bx
    integer, intent(in) :: lo(:)
    call box_build_2(bx, lo, lo)
  end subroutine box_build_1

  !! Returns the dimensionality of the box _bx_.
  function box_dim(bx) result(r)
    type(box), intent(in) :: bx
    integer :: r
    r = bx%dim
  end function box_dim

  !! Returns a coarsened version of the box _bx_, coarsened by the vector
  !! _cv_, but only on those components with mask true
  function box_coarsen_v_m(bx, cv, mask) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: cv(:)
    logical, intent(in) :: mask(:)
    type(box) :: r
    r%dim = bx%dim
    where ( mask(1:bx%dim) ) 
       r%lo = int_coarsen(bx%lo, cv)
       r%hi = int_coarsen(bx%hi, cv)
    end where
  end function box_coarsen_v_m
  !! Returns a coarsened version of the box _bx_, coarsened by the vector
  !! _cv_.
  function box_coarsen_v(bx, cv) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: cv(:)
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) = int_coarsen(bx%lo(1:bx%dim), cv)
    r%hi(1:r%dim) = int_coarsen(bx%hi(1:bx%dim), cv)
  end function box_coarsen_v
  !! Returns a coarsened version of the box _bx_, coarsened by the integer
  !! _ci_, but only on components with mask true
  function box_coarsen_i_m(bx, ci, mask) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: ci
    logical, intent(in) :: mask(:)
    type(box) :: r
    r%dim = bx%dim
    where ( mask(1:bx%dim) )
       r%lo = int_coarsen(bx%lo, ci)
       r%hi = int_coarsen(bx%hi, ci)
    end where
  end function box_coarsen_i_m
  !! Returns a coarsened version of the box _bx_, coarsened by the integer
  !! _ci_.
  function box_coarsen_i(bx, ci) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: ci
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) = int_coarsen(bx%lo(1:bx%dim), ci)
    r%hi(1:r%dim) = int_coarsen(bx%hi(1:bx%dim), ci)
  end function box_coarsen_i

  function box_refine_v_m(bx, rv, mask) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: rv(:)
    logical, intent(in) :: mask(:)
    type(box) :: r
    r%dim = bx%dim
    where ( mask(1:bx%dim) )
       r%lo =  bx%lo*rv
       r%hi = (bx%hi+1)*rv-1
    end where
  end function box_refine_v_m
  function box_refine_v(bx, rv) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: rv(:)
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) =  bx%lo(1:bx%dim)*rv
    r%hi(1:r%dim) = (bx%hi(1:bx%dim)+1)*rv-1
  end function box_refine_v
  function box_refine_i_m(bx, ri, mask) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: ri
    logical, intent(in) :: mask(:)
    type(box) :: r
    r%dim = bx%dim
    where ( mask(1:bx%dim) )
       r%lo =  bx%lo*ri
       r%hi = (bx%hi+1)*ri-1
    end where
  end function box_refine_i_m
  function box_refine_i(bx, ri) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: ri
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) =  bx%lo(1:bx%dim)*ri
    r%hi(1:r%dim) = (bx%hi(1:bx%dim)+1)*ri-1
  end function box_refine_i

  !! Box_nodalize, perhaps inaptly named, grows a box in the
  !! on the hiside by one if the correspending directions of
  !! nodal are true, or else not at all if nodal is not present.
  !! this last because by default fabs, multifabs, and etc, are cell-
  !! centered.
  function box_nodalize(bx, nodal) result(r)
    type(box) :: r
    type(box), intent(in) :: bx
    logical, intent(in), optional :: nodal(:)
    integer :: i
    r = bx
    if ( .not. present(nodal) ) return
    do i = 1, bx%dim
       if ( nodal(i) ) r = grow(r, 1, i, +1)
    end do
  end function box_nodalize

!   function box_grow_m(bx, mat) result(r)
!     type(box), intent(in) :: bx
!     integer, intent(in) :: mat(:,:)
!     type(box) :: r
!     r%dim = bx%dim
!     r%lo(1:r%dim) =  bx%lo(1:bx%dim) - mat(1:bx%dim,1)
!     r%hi(1:r%dim) =  bx%hi(1:bx%dim) + mat(1:bx%dim,2)
!   end function box_grow_m
  function box_grow_v_f(bx, rv, face) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: rv(:)
    integer, intent(in) :: face
    type(box) :: r
    r = bx
    if ( face == -1 ) then
       r = bx
       r%lo(1:r%dim) =  bx%lo(1:bx%dim) - rv
    else if ( face == 1 ) then
       r = bx
       r%hi(1:r%dim) =  bx%hi(1:bx%dim) + rv
    else 
       call bl_error("BOX_GROW_V_F: unexpected face: ", face)
    end if
  end function box_grow_v_f
  function box_grow_n_d_f(bx, ri, dim, face) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: ri, dim
    integer, intent(in) :: face
    type(box) :: r
    r = bx
    if ( face == -1 ) then
       r = bx
       r%lo(dim) = bx%lo(dim) - ri
    else if ( face == 1 ) then
       r = bx
       r%hi(dim) = bx%hi(dim) + ri
    else 
       call bl_error("BOX_GROW_N_D_F: unexpected face: ", face)
    end if
  end function box_grow_n_d_f
  function box_grow_n_f(bx, ri, face) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: ri
    integer, intent(in) :: face
    type(box) :: r
    r = bx
    if ( face == -1 ) then
       r%lo(1:r%dim) = bx%lo(1:bx%dim) - ri
    else if ( face == 1 ) then
       r%hi(1:r%dim) = bx%hi(1:bx%dim) + ri
    else 
       call bl_error("BOX_GROW_N_F: unexpected face: ", face)
    end if
  end function box_grow_n_f
  function box_grow_v_m(bx, rv, mask) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: rv(:)
    logical, intent(in) :: mask(:)
    type(box) :: r
    r%dim = bx%dim
    where ( mask(1:bx%dim) ) 
       r%lo = bx%lo - rv
       r%hi = bx%hi + rv
    end where
  end function box_grow_v_m
  function box_grow_v(bx, rv) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: rv(:)
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) =  bx%lo(1:bx%dim) - rv
    r%hi(1:r%dim) =  bx%hi(1:bx%dim) + rv
  end function box_grow_v
  function box_grow_n_m(bx, ri, mask) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: ri
    logical, intent(in) :: mask(:)
    type(box) :: r
    r%dim = bx%dim
    where ( mask(1:bx%dim) )
       r%lo = bx%lo - ri
       r%hi = bx%hi + ri
    end where
  end function box_grow_n_m
  function box_grow_n(bx, ri) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: ri
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) = bx%lo(1:bx%dim) - ri
    r%hi(1:r%dim) = bx%hi(1:bx%dim) + ri
  end function box_grow_n

  
  !! lexigraphic ordering in the lower-bound
  function box_less(bx1, bx2) result(r)
    logical :: r
    type(box), intent(in) :: bx1, bx2
    integer :: i
    do i = 1, bx1%dim
       if ( bx1%lo(i) == bx2%lo(i) ) cycle
       if ( bx1%lo(i) < bx2%lo(i) ) r = .TRUE.
       if ( bx1%lo(i) > bx2%lo(i) ) r = .FALSE.
       return
    end do
    r = .FALSE.
  end function box_less

  function box_equal(bx1, bx2) result(r)
    type(box), intent(in) :: bx1, bx2
    logical :: r
    r = all(bx1%lo == bx2%lo .and. bx1%hi == bx2%hi)
  end function box_equal
  function box_not_equal(bx1, bx2) result(r)
    type(box), intent(in) :: bx1, bx2
    logical :: r
    r = .NOT. box_equal(bx1, bx2)
  end function box_not_equal

  function box_intersection(bx1, bx2) result(r)
    type(box), intent(in) :: bx1, bx2
    type(box) :: r
    r%dim = min(bx1%dim, bx2%dim)
    r%lo(1:r%dim) = max(bx1%lo(1:r%dim), bx2%lo(1:r%dim))
    r%hi(1:r%dim) = min(bx1%hi(1:r%dim), bx2%hi(1:r%dim))
  end function box_intersection

  function box_intersects(bx1, bx2) result(r)
    type(box), intent(in) :: bx1, bx2
    logical :: r
    integer :: dm
    dm = min(bx1%dim, bx2%dim)
    r = all( (min(bx1%hi(1:dm),bx2%hi(1:dm)) &
         - max(bx1%lo(1:dm),bx2%lo(1:dm))) >= 0 )
  end function box_intersects

  function box_contains(bx1, bx2) result(r)
    type(box), intent(in) :: bx1, bx2
    logical :: r
    r = all(bx1%lo(1:bx1%dim) <= bx2%lo(1:bx2%dim)) &
         .and. all(bx1%hi(1:bx1%dim) >= bx2%hi(1:bx2%dim))
  end function box_contains

  function box_contains_strict(bx1, bx2) result(r)
    type(box), intent(in) :: bx1, bx2
    logical :: r
    r = all(bx1%lo(1:bx1%dim) < bx2%lo(1:bx2%dim)) &
         .and. all(bx1%hi(1:bx1%dim) > bx2%hi(1:bx2%dim))
  end function box_contains_strict

  function box_contains_iv(bx1, iv) result(r)
    type(box), intent(in) :: bx1
    integer, intent(in) :: iv(:)
    logical :: r
    r = all(bx1%lo(1:bx1%dim) <= iv) .and. all(bx1%hi(1:bx1%dim) >= iv)
  end function box_contains_iv

  function box_shift_v_m(bx, iv, mask) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: iv(:)
    logical, intent(in) :: mask(:)
    type(box) :: r
    r%dim = bx%dim
    where ( mask(1:bx%dim) )
       r%lo = bx%lo + iv
       r%hi = bx%hi + iv
    end where
  end function box_shift_v_m
  function box_shift_v(bx, iv) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: iv(:)
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) = bx%lo(1:r%dim) + iv(1:r%dim)
    r%hi(1:r%dim) = bx%hi(1:r%dim) + iv(1:r%dim)
  end function box_shift_v
  function box_shift_i_m(bx, i, mask) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: i
    logical, intent(in) :: mask(:)
    type(box) :: r
    r%dim = bx%dim
    where ( mask(1:bx%dim) )
       r%lo = bx%lo + i
       r%hi = bx%hi + i
    end where
  end function box_shift_i_m
  function box_shift_i(bx, i) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: i
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) = bx%lo(1:r%dim) + i
    r%hi(1:r%dim) = bx%hi(1:r%dim) + i
  end function box_shift_i
  function box_shift_d(bx, i, dim) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: i
    integer, intent(in) :: dim
    type(box) :: r
    r = bx
    r%lo(dim) = bx%lo(dim) + i
    r%hi(dim) = bx%hi(dim) + i
  end function box_shift_d

  function box_extent(bx) result(r)
    type (box), intent(in) :: bx
    integer, dimension(bx%dim) :: r
    r = bx%hi(1:bx%dim) - bx%lo(1:bx%dim) + 1
  end function box_extent

  function box_extent_d(bx, dim) result(r)
    type (box), intent(in) :: bx
    integer, intent(in) :: dim
    integer :: r
    r = bx%hi(dim)-bx%lo(dim)+1
  end function box_extent_d

  function box_lwb(bx) result(r)
    type (box), intent(in) :: bx
    integer, dimension(bx%dim) :: r
    r = Huge(1)
    r(1:bx%dim) = bx%lo(1:bx%dim)
  end function box_lwb
  function box_lwb_d(bx, dim) result(r)
    type (box), intent(in) :: bx
    integer, intent(in) :: dim
    integer :: r
    r = bx%lo(dim)
  end function box_lwb_d
  subroutine box_set_lwb(bx, iv)
    integer, intent(in), dimension(:) :: iv
    type (box), intent(inout) :: bx
    bx%lo = iv
  end subroutine box_set_lwb
  subroutine box_set_lwb_d(bx, dim, v)
    integer, intent(in) :: v
    integer, intent(in) :: dim
    type (box), intent(inout) :: bx
    bx%lo(dim) = v
  end subroutine box_set_lwb_d

  function box_upb(bx) result(r)
    type (box), intent(in) :: bx
    integer, dimension(bx%dim) :: r
    r = -Huge(1)
    r(1:bx%dim) = bx%hi(1:bx%dim)
  end function box_upb
  function box_upb_d(bx, dim) result(r)
    type (box), intent(in) :: bx
    integer, intent(in) :: dim
    integer :: r
    r = bx%hi(dim)
  end function box_upb_d
  subroutine box_set_upb(bx, iv)
    integer, intent(in), dimension(:) :: iv
    type (box), intent(inout) :: bx
    bx%hi = iv
  end subroutine box_set_upb
  subroutine box_set_upb_d(bx, dim, v)
    integer, intent(in) :: v
    integer, intent(in) :: dim
    type (box), intent(inout) :: bx
    bx%hi(dim) = v
  end subroutine box_set_upb_d

  function box_mod(bx1, bx2, pmask) result(r)
    type(box) :: r
    type(box), intent(in) :: bx1
    type(box), intent(in) :: bx2
    logical, intent(in) :: pmask(:)
    integer :: dm
    integer, dimension(bx1%dim) :: l1, l2, h1, h2, a1, a2
    dm = bx1%dim
    l1 = bx1%lo(1:dm); h1 = bx1%hi(1:dm); a1 = h1-l1+1
    l2 = bx2%lo(1:dm); h2 = bx2%hi(1:dm); a2 = h2-l2+1
    where ( pmask )
       l1 = modulo(l1-l2, a2) + l2
       h1 = l1 + a1 - 1
    end where
    call build(r, l1, h1)
  end function box_mod

  subroutine box_decompose_mod(boxes, n, bx, bxi, pd, pmask)
    type(box), intent(in) :: bx
    type(box), intent(in) :: bxi
    type(box), intent(in) :: pd
    type(box), intent(out) :: boxes(:,:)
    logical, intent(in) :: pmask(:)
    integer, intent(out) :: n
    integer :: i
    call box_decompose(boxes(:,1), n, bx, bxi)
    do i = 1, n
       boxes(i,2) = box_mod(boxes(i,1), pd, pmask)
    end do
  end subroutine box_decompose_mod
  
  !! Box_decompose decompoes a box, bx1, into boxes that are
  !! either completely inside, or completely outside the regions of
  !! index spaces delineated by bx2.  Up to 3**dim boxes are returned;
  !! and failure results if not enough space is passed.

  subroutine box_decompose(boxes, n, bx1, bx2)
    type(box), intent(in) :: bx1
    type(box), intent(in) :: bx2
    type(box), intent(out) :: boxes(:)
    integer, intent(out) :: n
    type(box) :: r(27)
    type(box) :: bxl, bxr
    integer :: dm, i

    dm = bx1%dim
    bxl = bx1
    call box_chop(bxl, r(1), bxr,  1, lwb(bx2,1))
    call box_chop(bxr, r(2), r(3), 1, upb(bx2,1)+1)
    if ( dm >= 2 ) then
       r(1:9:3) = r(1:3)
       do i = 1, 9, 3
          bxl = r(i)
          call box_chop(bxl, r(i  ), bxr,    2, lwb(bx2,2))
          call box_chop(bxr, r(i+1), r(i+2), 2, upb(bx2,2)+1)
       end do
    end if
    if ( dm >= 3 ) then
       r(1:27:3) = r(1:9)
       do i = 1, 27, 3
          bxl = r(i)
          call box_chop(bxl, r(i  ), bxr,    3, lwb(bx2,3))
          call box_chop(bxr, r(i+1), r(i+2), 3, upb(bx2,3)+1)
       end do
    end if
    n = 0
    do i = 1, 3**dm
       if ( empty(r(i)) ) cycle
       n = n + 1
       if ( n > size(boxes) ) then
          call bl_error("BOX_DECOMPOSE: no room")
       end if
       boxes(n) = r(i)
    end do

  end subroutine box_decompose

  !! Box_chop chops a box into two boxes in the coordinate direction
  !! dir at the cell ichop.  The resulting boxes are  bxl and bxr where
  !! bxl is to the left of the chop line and bxr is to the right.  If
  !! the chop line is outside the extent of the box; the coresponding 
  !! is an empty, default box.
  subroutine box_chop(bx, bxl, bxr, dim, ichop)
    type(box), intent(in) :: bx
    type(box), intent(out) :: bxl, bxr
    integer, intent(in) :: dim, ichop
    if ( ichop <= lwb(bx,dim) ) then
       bxr = bx
    else if ( ichop > upb(bx,dim) ) then
       bxl = bx
    else
       bxl = bx
       bxr = bx
       bxl%hi(dim) = ichop-1
       bxr%lo(dim) = ichop
    end if
  end subroutine box_chop

  subroutine box_peel(bx, dim, face, thk, bx1)
    type(box), intent(inout) :: bx
    integer, intent(in) :: dim, thk, face
    type(box), intent(out) :: bx1
    bx1 = bx
    if ( face == -1 ) then
       bx%lo(dim)  = bx%lo(dim)  +  thk
       bx1%hi(dim) = bx1%lo(dim) + (thk-1)
    else if ( face == 1 ) then
       bx%hi(dim)  = bx%hi(dim)  -  thk
       bx1%lo(dim) = bx1%hi(dim) - (thk-1)
    else
       call bl_error("BOX_PEEL: called with face /= {-1,+1}: ", face)
    end if
  end subroutine box_peel

  elemental function box_volume(bx) result(r)
    type(box), intent(in) :: bx
    integer(kind=ll_t) :: r
    integer :: i
    r = 1_ll_t
    do i = 1, bx%dim
       r = r*(bx%hi(i)-bx%lo(i)+1)
    end do
  end function box_volume

  elemental function box_dvolume(bx) result(r)
    type(box), intent(in) :: bx
    real(dp_t) :: r
    integer :: i
    r = 1
    do i = 1, bx%dim
       r = r*real(bx%hi(i)-bx%lo(i)+1, kind=dp_t)
    end do
  end function box_dvolume

  function box_bbox(b1, b2) result(r)
    type(box), intent(in) :: b1, b2
    type(box) :: r
    r%dim = b1%dim
    r%lo(1:r%dim) = min(b1%lo(1:b1%dim),b2%lo(1:b2%dim))
    r%hi(1:r%dim) = max(b1%hi(1:b1%dim),b2%hi(1:b2%dim))
  end function box_bbox

  function box_boundary_n(bx, igrow) result(r)
    type(box), intent(in) :: bx
    integer, intent(in)  :: igrow
    type(box), dimension(MAX_SPACEDIM, 2) :: r
    type(box) :: tbx
    integer :: i, j, jj
    if ( igrow < 0 ) then
       call bl_error('box_boundary_n: grow < 0: ', igrow)
    end if
    tbx = grow(bx, igrow)
    do i = 1, bx%dim
       do j = -1, 1, 2
          jj = (3+j)/2
          call box_peel(tbx, i, j, igrow, r(i,jj))
       end do
    end do
    if ( tbx /= bx ) then
       call bl_error('box_boundary_n: tbx /= bx')
    end if
  end function box_boundary_n

!   function box_boundary_m(bx, mat) result(r)
!     type(box), intent(in) :: bx
!     integer, intent(in)  :: mat(:,:)
!     type(box), dimension(MAX_SPACEDIM, 2) :: r
!     type(box) :: tbx
!     integer :: i, j, jj
!     tbx = grow(bx, mat)
!     do i = 1, bx%dim
!        do j = -1, 1, 2
!           jj = (3+j)/2
!           call box_peel(tbx, i, j, mat(i,jj), r(i,jj))
!        end do
!     end do
!     if ( tbx /= bx ) then
!        call bl_error('box_boundary_all: tbx /= bx')
!     end if
!   end function box_boundary_m

  function box_reduce(bx, dim) result(r)
    type(box), intent(in) :: bx
    type(box) :: r
    integer, intent(in) :: dim
    r%dim = bx%dim-1
    r%lo(:dim-1) = bx%lo(:dim-1)
    r%lo(dim:r%dim)   = bx%lo(dim+1:bx%dim)
    r%hi(:dim-1) = bx%hi(:dim-1)
    r%hi(dim:r%dim)   = bx%hi(dim+1:bx%dim)
  end function box_reduce

  function box_touch_d(a, b, dim, face) result(r)
    type(box), intent(in) :: a, b
    integer, intent(in) :: dim, face
    logical :: r
    if ( face == -1 ) then
       r = a%lo(dim) == b%hi(dim) + 1
    else if ( face == 1 ) then
       r = a%hi(dim) + 1 == b%lo(dim)
    end if
    if ( a%dim > 1 ) then
       r = r .and. intersects(reduce(a,dim), reduce(b,dim))
    end if
  end function box_touch_d

  function box_touch(a, b, dim) result(r)
    type(box), intent(in) :: a, b
    integer, intent(in) :: dim
    logical :: r
    r = a%hi(dim) + 1 == b%lo(dim) .or. a%lo(dim) == b%hi(dim) + 1
    if ( a%dim > 1 ) then
       r = r .and. intersects(reduce(a,dim), reduce(b,dim))
    end if
  end function box_touch

  function box_abut_d(a, b, dim, face) result(r)
    type(box), intent(in) :: a, b
    integer, intent(in) :: dim, face
    logical :: r
    if ( face == -1 ) then
       r = a%lo(dim) == b%hi(dim) + 1
    else if ( face == 1 ) then
       r = a%hi(dim) + 1 == b%lo(dim)
    end if
    if ( a%dim > 1 ) then
       r = r .and. reduce(a,dim) == reduce(b,dim)
    end if
  end function box_abut_d

  function box_merge(bb, b1, b2, dim) result(r)
    logical :: r
    type(box), intent(out) :: bb
    type(box), intent(in)  :: b1, b2
    integer, intent(in) :: dim
    integer :: lo(b1%dim), hi(b1%dim)
    if ( b1 == b2 ) then
       bb = b1
       r = .true.
       return
    end if
    r = box_abut(b1, b2, dim)
    if ( r ) then
       lo = lwb(b1); hi = upb(b1)
       if ( lwb(b1, dim) < lwb(b2, dim) ) then
          hi(dim) = upb(b2, dim)
       else
          lo(dim) = lwb(b2, dim)
       end if
       bb = make_box(lo, hi)
    end if
  end function box_merge

  function box_abut(a, b, dim) result(r)
    type(box), intent(in) :: a, b
    integer, intent(in) :: dim
    logical :: r
    r = a%hi(dim) + 1 == b%lo(dim) .or. a%lo(dim) == b%hi(dim) + 1
    if ( a%dim > 1 ) then
       r = r .and. reduce(a,dim) == reduce(b,dim)
    end if
  end function box_abut

  subroutine box_print(bx, str, unit, advance, legacy, nodal, skip)
    use bl_IO_module
    use bl_string_module
    type(box), intent(in) :: bx
    character(len=*), intent(in), optional :: str
    character(len=*), intent(in), optional :: advance
    integer, intent(in), optional :: skip
    logical, intent(in), optional :: legacy
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: nodal(:)
    integer :: un
    character(len=3) :: adv
    integer :: i
    logical :: llegacy
!   logical :: lfixw
!   integer :: maxdig
    character(len=2) :: findex
    ! how many digits, including possible sign
!    maxdig = max(maxval(abs(lwb(bx))),maxval(upb(bx)))
!    maxdig = ceiling(log10(real(maxdig+1)))
!    lfixw = maxdig <= box_print_fixed_width
!    if ( lfixw ) then
!       write(findex, fmt='("I",i1)') box_print_fixed_width
!    else
       findex = 'I0'
!    end if
    un = unit_stdout(unit)
    adv = unit_advance(advance)
    llegacy = .FALSE.; if ( present(legacy) ) llegacy = legacy
    if ( llegacy ) then
       call write_a_legacy_box(nodal)
    else
       call write_a_box()
    end if
    if ( eq_i(adv, 'YES') ) write(unit=un, fmt='()')
  contains
    subroutine write_a_box
      call unit_skip(un, skip)
      write(unit=un, fmt='("BOX")', advance='no')
      if ( present(str) ) then
         write(unit=un, fmt='("{", A, "}")', advance = 'no') str
      end if
      if ( bx%dim > 0 ) then
         write(unit=un, fmt='("[{", 3('//findex//',:,", "))', advance = 'no') bx%lo(1:bx%dim)
         write(unit=un, fmt='("}, {", 3('//findex//',:,", "))', advance = 'no') bx%hi(1:bx%dim)
         write(unit=un, fmt='("}]")', advance = 'no' )
      end if
    end subroutine write_a_box
    subroutine write_a_legacy_box(nodal)
      logical, intent(in), optional :: nodal(:)
      integer :: tp(bx%dim), thi(bx%dim)
      tp = 0
      if ( present(nodal) ) then
         do i = 1, bx%dim
            if ( nodal(i) ) tp(i) = 1
         end do
      end if
      if ( bx%dim == 0 ) then
         write(unit=un, fmt='("()")', advance = 'no')
      else
         thi = bx%hi(1:bx%dim)
         if ( present(nodal) ) then
            where ( nodal ) thi = thi + 1
         end if
         write(unit=un, fmt='("(")', advance='no')

         call write_a_legacy_intvect(bx%lo); write(unit=un, fmt='(" ")', advance='no')
         call write_a_legacy_intvect(thi); write(unit=un, fmt='(" ")', advance='no')
         call write_a_legacy_intvect(tp)
         write(unit=un, fmt='(")")', advance='no')
      end if
    end subroutine write_a_legacy_box

    subroutine write_a_legacy_intvect(vc)
      integer, intent(in) :: vc(:)
      write(unit=un, fmt='("(")', advance='no')
      do i = 1, bx%dim-1
         write(unit=un, fmt='('//findex//', ",")', advance = 'no' ) vc(i)
      end do
      write(unit=un, fmt='('//findex//', ")")', advance = 'no' ) vc(bx%dim)
    end subroutine write_a_legacy_intvect
  end subroutine box_print

  subroutine box_read(bx, unit, stat)
    use bl_IO_module
    use bl_stream_module
    type(box), intent(out) :: bx
    type(bl_stream) :: strm
    integer, intent(in), optional :: unit
    logical, intent(out), optional :: stat
    character :: c
    if ( present(stat) ) stat = .true.
    call build(strm, unit_stdin(unit))
    c = bl_stream_peek_chr(strm)
    if ( c == 'b' ) then
       call read_box()
    else if ( c == '(' ) then
       call read_a_legacy_box()
    else 
       call bl_error("BOX_READ: expected b or (, got ", c)
    end if
    call destroy(strm)
  contains
    subroutine read_box
       call bl_stream_expect(strm, "box")
       call bl_stream_expect(strm, (/"[", "{"/))
       call bl_stream_expect(strm, (/"}", "]"/))
    end subroutine read_box
    subroutine read_a_legacy_box
      character(len=1) :: c
      integer :: i
      integer :: dim
      integer :: it(MAX_SPACEDIM)

      call bl_stream_expect(strm, '(')
      call bl_stream_expect(strm, '(')
      bx%lo(1) = bl_stream_scan_int(strm)
      do i = 2, MAX_SPACEDIM+1
         c = bl_stream_scan_chr(strm)
         if ( c == ')') then
            call bl_stream_putback_chr(strm, c)
            dim = i-1
            exit
         else if ( c /= ',' ) then
            call bl_error("READ_BOX: expected ',' or ')' got: ", '"'//c//'"')
         end if
         bx%lo(i) = bl_stream_scan_int(strm)
      end do
      call bl_stream_expect(strm, ')')

      call bl_stream_expect(strm, '(')
      bx%hi(1) = bl_stream_scan_int(strm)
      do i = 2, dim
         call bl_stream_expect(strm, ',')
         bx%hi(i) = bl_stream_scan_int(strm)
      end do
      call bl_stream_expect(strm, ')')

      it = 0
      call bl_stream_expect(strm, '(')
      it(1) = bl_stream_scan_int(strm)
      do i = 2, dim
         call bl_stream_expect(strm, ',')
         it(i) = bl_stream_scan_int(strm)
      end do
      call bl_stream_expect(strm, ')')
      call bl_stream_expect(strm, ')')
      if ( any(it /= 0) ) then
         if ( present(stat) ) then
            call bl_warn("READ_BOX: throws away information: a nodal box")
            stat = .false.
            return
         else
            call bl_error("READ_BOX: throws away information: a nodal box")
         end if
      end if
      bx%dim = dim

    end subroutine read_a_legacy_box

  end subroutine box_read

  subroutine box_periodic_shift(dmn,b,nodal,pmask,ng,shft,bxs,cnt)

    type(box), intent(in)  :: dmn,b
    logical,   intent(in)  :: nodal(:),pmask(:)
    integer,   intent(in)  :: ng
    integer,   intent(out) :: shft(:,:),cnt
    type(box), intent(out) :: bxs(:)

    type(box) :: domain, bx, src
    integer   :: nbeg(3),nend(3),ldomain(3),r(3),ri,rj,rk

    cnt = 0

    if (all(pmask .eqv. .false.)) return

    bx     = grow(box_nodalize(b,nodal),ng)
    domain = box_nodalize(dmn,nodal)

    if (contains(domain,bx)) return

    nbeg = 0
    nend = 0
    nbeg(1:bx%dim) = -1
    nend(1:bx%dim) = +1

    ldomain = (/ extent(domain,1),extent(domain,2),extent(domain,3) /)

    do ri = nbeg(1), nend(1)
       if (ri /= 0 .and. (.not. is_periodic(1))) cycle
       if (ri /= 0 .and. is_periodic(1)) bx = shift(bx,ri*ldomain(1),1)

       do rj = nbeg(2), nend(2)
          if (rj /= 0 .and. (.not. is_periodic(2))) cycle
          if (rj /= 0 .and. is_periodic(2)) bx = shift(bx,rj*ldomain(2),2)

          do rk = nbeg(3), nend(3)
             if (rk /= 0 .and. (.not. is_periodic(3))) cycle
             if (rk /= 0 .and. is_periodic(3)) bx = shift(bx,rk*ldomain(3),3)

             if (ri == 0 .and. rj == 0 .and. rk == 0) cycle

             src = intersection(domain,bx)

             if (.not. empty(src)) then
                cnt = cnt + 1
                bxs(cnt) = src
                r = (/ri,rj,rk/)
                r = r*ldomain
                shft(cnt,1:bx%dim) = r(1:bx%dim)
             end if

             if (rk /= 0 .and. is_periodic(3)) bx = shift(bx,-rk*ldomain(3),3)
          end do

          if (rj /= 0 .and. is_periodic(2)) bx = shift(bx,-rj*ldomain(2),2)
       end do

       if (ri /= 0 .and. is_periodic(1)) bx = shift(bx,-ri*ldomain(1),1)
    end do

    contains

      function is_periodic(i) result(r)
        integer, intent(in) :: i
        logical             :: r
        r = .false.
        if (i >= 1 .and. i <= bx%dim) r = pmask(i) .eqv. .true.
      end function is_periodic

  end subroutine box_periodic_shift

end module box_module
