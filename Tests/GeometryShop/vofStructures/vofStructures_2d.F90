
module stencil_test_module

  use amrex_fort_module, only : dim=>bl_spacedim
  implicit none

  type, bind(c) :: ccgn
     integer :: Nnbr(2,dim)
     integer :: nbr(4,2,dim)
     integer :: ebID
  end type ccgn

  type, bind(c) :: node
     integer :: nCells
     integer :: iv(dim)
     type(ccgn) :: cells(NCELLMAX)
  end type node

  private

  public :: do_eb_work

contains

  subroutine do_eb_work(nodes,num,&
       mask,mask_l1,mask_l2,mask_h1,mask_h2) bind(C,name="do_eb_work")

    integer, intent(in   ) :: num
    integer, intent(in   ) :: mask_l1,mask_l2,mask_h1,mask_h2
    integer, intent(in   ) :: mask(mask_l1:mask_h1,mask_l2:mask_h2)
    type(node), intent(in) :: nodes(num) 

    integer :: i

    do i = 1, num
       print *, nodes(i)%nCells, nodes(i)%iv, nodes(i)%cells(1)%ebID
    end do

  end subroutine do_eb_work

end module stencil_test_module
