module ml_layout_module

  use layout_module
  use multifab_module
  use ml_boxarray_module

  implicit none

  type ml_layout
     integer :: dim = 0
     integer :: nlevel = 0
     type(ml_boxarray) :: mba
     type(layout), pointer :: la(:) => Null()
     type(lmultifab), pointer :: mask(:) => Null()
  end type ml_layout

  interface destroy
     module procedure ml_layout_destroy
  end interface

  interface operator(.eq.)
     module procedure ml_layout_equal
  end interface
  interface operator(.ne.)
     module procedure ml_layout_not_equal
  end interface

contains

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

  subroutine ml_layout_destroy(mla)
    type(ml_layout), intent(inout) :: mla
    integer :: n
    do n = 1, mla%nlevel-1
       call destroy(mla%mask(n))
    end do
    call destroy(mla%mba)
    call destroy(mla%la(1))
    deallocate(mla%la, mla%mask)
    mla%dim = 0
    mla%nlevel = 0
  end subroutine ml_layout_destroy

  function ml_get_layout(mla, n) result(r)
    type(layout) :: r
    type(ml_layout), intent(in) :: mla
    integer, intent(in) :: n
    r = mla%la(n)
  end function ml_get_layout

  function ml_get_pd(mla, n) result(r)
    type(box) :: r
    type(ml_layout), intent(in) :: mla
    integer, intent(in) :: n
    r = ml_boxarray_get_pd(mla%mba, n)
  end function ml_get_pd

  subroutine ml_layout_build(mla, mba)
    type(ml_layout), intent(inout) :: mla
    type(ml_boxarray), intent(in) :: mba
    type(boxarray) :: bac
    integer :: n
    mla%nlevel = mba%nlevel
    mla%dim = mba%dim
    call copy(mla%mba, mba)
    allocate(mla%la(mla%nlevel))
    allocate(mla%mask(mla%nlevel-1))
    call build(mla%la(1), mba%bas(1))
    do n = 2, mba%nlevel
       call layout_build_pn(mla%la(n), mla%la(n-1), mba%bas(n), mba%rr(n-1,:))
    end do
    do n = mba%nlevel-1,  1, -1
       call lmultifab_build(mla%mask(n), mla%la(n), nc = 1, ng = 0)
       call setval(mla%mask(n), val = .TRUE.)
       call copy(bac, mba%bas(n+1))
       call boxarray_coarsen(bac, mba%rr(n,:))
       call setval(mla%mask(n), .false., bac)
       call destroy(bac)
    end do
  end subroutine ml_layout_build

end module ml_layout_module
