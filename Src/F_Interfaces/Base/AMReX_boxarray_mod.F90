
module amrex_boxarray_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_long
  use amrex_box_module

  implicit none

  private

  public :: amrex_boxarray_build, amrex_boxarray_destroy, amrex_print

  type, public :: amrex_boxarray
     logical     :: owner = .false.
     type(c_ptr) :: p = c_null_ptr
   contains
     generic   :: assignment(=) => amrex_boxarray_assign, amrex_boxarray_install  ! shallow copy
     procedure :: clone         => amrex_boxarray_clone    ! deep copy
     procedure :: move          => amrex_boxarray_move     ! transfer ownership
     generic   :: maxSize       => amrex_boxarray_maxsize_int, &  ! make the boxes smaller
          &                        amrex_boxarray_maxsize_int3, amrex_boxarray_maxsize_iv
     procedure :: get_box       => amrex_boxarray_get_box
     procedure :: nodal_type    => amrex_boxarray_nodal_type  ! get index type
     procedure :: num_pts       => amrex_boxarray_num_pts
     procedure :: intersects    => amrex_boxarray_intersects_box
     procedure, private :: amrex_boxarray_assign
     procedure, private :: amrex_boxarray_install
     procedure, private :: amrex_boxarray_maxsize_int
     procedure, private :: amrex_boxarray_maxsize_int3
     procedure, private :: amrex_boxarray_maxsize_iv
     procedure, private :: amrex_boxarray_intersects_box
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_boxarray_destroy
#endif
  end type amrex_boxarray

  interface amrex_boxarray_build
     module procedure amrex_boxarray_build_bx
     module procedure amrex_boxarray_build_bxs
  end interface amrex_boxarray_build

  interface amrex_print
     module procedure amrex_boxarray_print
  end interface amrex_print

  ! interfaces to cpp functions

  interface
     subroutine amrex_fi_new_boxarray (ba,lo,hi) bind(c)
       import
       implicit none
       type(c_ptr) :: ba
       integer(c_int), intent(in) :: lo(3), hi(3)
     end subroutine amrex_fi_new_boxarray

     subroutine amrex_fi_new_boxarray_from_bxfarr (ba,bxs,nsides,ndims,nbxs) bind(c)
       import
       implicit none
       type(c_ptr) :: ba
       integer(c_int),intent(in) :: bxs(*)
       integer(c_int), value :: nsides, ndims, nbxs
     end subroutine amrex_fi_new_boxarray_from_bxfarr

     subroutine amrex_fi_delete_boxarray (ba) bind(c)
       import
       implicit none
       type(c_ptr), value :: ba
     end subroutine amrex_fi_delete_boxarray

     subroutine amrex_fi_clone_boxarray (bao, bai) bind(c)
       import
       implicit none
       type(c_ptr) :: bao
       type(c_ptr), value :: bai
     end subroutine amrex_fi_clone_boxarray

     subroutine amrex_fi_boxarray_maxsize (ba,s) bind(c)
       import
       implicit none
       type(c_ptr), value :: ba
       integer(c_int), intent(in) :: s(3)
     end subroutine amrex_fi_boxarray_maxsize

     subroutine amrex_fi_boxarray_get_box (ba,i,lo,hi) bind(c)
       import
       implicit none
       type(c_ptr), value :: ba
       integer(c_int), value :: i
       integer, intent(inout) :: lo(3), hi(3)
     end subroutine amrex_fi_boxarray_get_box

     subroutine amrex_fi_print_boxarray (ba) bind(c)
       import
       implicit none
       type(c_ptr), value :: ba
     end subroutine amrex_fi_print_boxarray

     pure subroutine amrex_fi_boxarray_nodal_type (ba, inodal) bind(c)
       import
       implicit none
       type(c_ptr), value, intent(in) :: ba
       integer, intent(inout) :: inodal(3)
     end subroutine amrex_fi_boxarray_nodal_type

     pure function amrex_fi_boxarray_numpts (ba) bind(c)
       import
       implicit none
       type(c_ptr), value, intent(in) :: ba
       integer(amrex_long) :: amrex_fi_boxarray_numpts
     end function amrex_fi_boxarray_numpts

     pure integer function amrex_fi_boxarray_intersects_box (ba, lo, hi) bind(c)
       import
       implicit none
       type(c_ptr), value, intent(in) :: ba
       integer, intent(in) :: lo(*), hi(*)
     end function amrex_fi_boxarray_intersects_box
  end interface

contains

  subroutine amrex_boxarray_build_bx (ba, bx)
    type(amrex_boxarray) :: ba
    type(amrex_box), intent(in ) :: bx
    ba%owner = .true.
    call amrex_fi_new_boxarray(ba%p, bx%lo, bx%hi)
  end subroutine amrex_boxarray_build_bx

  subroutine amrex_boxarray_build_bxs (ba, bxs)
    type(amrex_boxarray) :: ba
    integer,intent(in) :: bxs(:,:,:) ! (lo:hi,dim,#ofboxs)
    ba%owner = .true.
    call amrex_fi_new_boxarray_from_bxfarr(ba%p, bxs, size(bxs,1), size(bxs,2), size(bxs,3))
  end subroutine amrex_boxarray_build_bxs

  impure elemental subroutine amrex_boxarray_destroy (this)
    type(amrex_boxarray), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_boxarray(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_boxarray_destroy

  subroutine amrex_boxarray_assign (dst, src)
    class(amrex_boxarray), intent(inout) :: dst
    type (amrex_boxarray), intent(in   ) :: src
    call amrex_boxarray_destroy(dst)
    dst%owner = .false.
    dst%p = src%p
  end subroutine amrex_boxarray_assign

  subroutine amrex_boxarray_install (this, p)
    class(amrex_boxarray), intent(inout) :: this
    type(c_ptr), intent(in) :: p
    this%owner = .false.
    this%p     = p
  end subroutine amrex_boxarray_install

  subroutine amrex_boxarray_clone (dst, src)
    class(amrex_boxarray), intent(inout) :: dst
    type (amrex_boxarray), intent(in   ) :: src
    dst%owner = .true.
    call amrex_fi_clone_boxarray(dst%p, src%p)
  end subroutine amrex_boxarray_clone

  subroutine amrex_boxarray_move (dst, src)
    class(amrex_boxarray) :: dst, src
    call amrex_boxarray_destroy(dst)
    dst%owner = src%owner
    dst%p = src%p
    src%owner = .false.
    src%p = c_null_ptr
  end subroutine amrex_boxarray_move

  subroutine amrex_boxarray_maxsize_int (this, s)
    class(amrex_boxarray) :: this
    integer, intent(in)   :: s
    call amrex_fi_boxarray_maxsize(this%p, [s,s,s])
  end subroutine amrex_boxarray_maxsize_int

  subroutine amrex_boxarray_maxsize_int3 (this, sx, sy, sz)
    class(amrex_boxarray) :: this
    integer, intent(in)   :: sx, sy, sz
    call amrex_fi_boxarray_maxsize(this%p, [sx,sy,sz])
  end subroutine amrex_boxarray_maxsize_int3

  subroutine amrex_boxarray_maxsize_iv (this, s)
    class(amrex_boxarray) :: this
    integer, intent(in)   :: s(3)
    call amrex_fi_boxarray_maxsize(this%p, s)
  end subroutine amrex_boxarray_maxsize_iv

  function amrex_boxarray_get_box (this, i) result(bx)
    class(amrex_boxarray) :: this
    integer, intent(in)   :: i
    type(amrex_box) :: bx
    call amrex_fi_boxarray_get_box(this%p, i, bx%lo, bx%hi)
  end function amrex_boxarray_get_box

  subroutine amrex_boxarray_print (ba)
    type(amrex_boxarray), intent(in) :: ba
    call amrex_fi_print_boxarray(ba%p)
  end subroutine amrex_boxarray_print

  pure function amrex_boxarray_nodal_type (this) result(nodal)
    class(amrex_boxarray), intent(in) :: this
    logical, dimension(3) :: nodal
    integer, dimension(3) :: inodal
    inodal = 0
    call amrex_fi_boxarray_nodal_type(this%p,inodal)
    nodal(1) = inodal(1) .ne. 0
    nodal(2) = inodal(2) .ne. 0
    nodal(3) = inodal(3) .ne. 0
  end function amrex_boxarray_nodal_type

  pure function amrex_boxarray_num_pts (this) result(n)
    class(amrex_boxarray), intent(in) :: this
    integer(amrex_long) :: n
    n = amrex_fi_boxarray_numpts(this%p)
  end function amrex_boxarray_num_pts

  pure function amrex_boxarray_intersects_box (this, bx) result(r)
    class(amrex_boxarray), intent(in) :: this
    type(amrex_box), intent(in) :: bx
    logical :: r
    integer :: ir
    ir = amrex_fi_boxarray_intersects_box(this%p, bx%lo, bx%hi)
    r = ir .ne. 0
  end function amrex_boxarray_intersects_box

end module amrex_boxarray_module
