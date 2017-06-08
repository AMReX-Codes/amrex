
module stencil_test_module

  use amrex_fort_module, only : dim=>bl_spacedim
  implicit none

  type, bind(c) :: ccgn
     integer :: Nnbr(0:1,0:dim-1)
     integer :: nbr(0:NCELLMAX-1,0:1,0:dim-1)
     integer :: ebID
  end type ccgn

  type, bind(c) :: node
     integer :: nCells
     integer :: iv(0:dim-1)
     type(ccgn) :: cells(0:NCELLMAX-1)
  end type node

  private

  public :: do_eb_work

contains

  subroutine do_eb_work(lo, hi, nodes, num,&
       mask, mask_l1, mask_l2, mask_h1, mask_h2) bind(C,name="do_eb_work")

    integer,    intent(in   ) :: lo(dim), hi(dim)
    integer,    intent(in   ) :: num
    integer,    intent(in   ) :: mask_l1,mask_l2,mask_h1,mask_h2
    integer,    intent(in   ) :: mask(mask_l1:mask_h1,mask_l2:mask_h2)
    type(node), intent(in   ) :: nodes(0:num-1) 

    integer :: i, j, mC, mL, nC, nL

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          mC = mask(i  ,j)
          mL = mask(i-1,j)
          if (mC .eq. -2) then         ! C is regular
             if (mL .eq. -2) then      ! L is regular
                nL = 0
             else if (mL .eq. -1) then ! L is covered
                call bl_pd_abort('ERROR: cannot have covered next to regular cell')
             else                      ! L is irregular
                nC = 0
                nL = nodes(mL) % nCells
                print *,i,j,'i r'
             endif
          else if (mC .ne. -1) then    ! C is irregular
             if (mC .lt. 0  .or.  mC .ge. num) then
                call bl_pd_abort('ERROR: cut cell index OOB')
             endif
             nC = nodes(mC) % nCells
             if (nC .gt. 1) then
                call bl_pd_abort('ERROR: Should not be any mv cut cells here')
             endif
             nL = nodes(mC) % cells(1) % Nnbr(0,1)
             print *,i,j,'r i'
          endif
       enddo
    enddo

  end subroutine do_eb_work

end module stencil_test_module
