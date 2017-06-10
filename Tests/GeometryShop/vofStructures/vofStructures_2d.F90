#include <Node.H>

module stencil_test_module

  use amrex_fort_module, only : dim=>bl_spacedim
  implicit none

  type, bind(c) :: cutcell
     integer :: Nnbr(0:1,0:dim-1)
     integer :: nbr(0:NCELLMAX-1,0:1,0:dim-1)
     integer :: faceID(0:NCELLMAX-1,0:1,0:dim-1)
     integer :: ebCellID
  end type cutcell

  type, bind(c) :: cnode
     integer :: nCells
     integer :: iv(0:dim-1)
     type(cutcell) :: cells(0:NCELLMAX-1)
  end type cnode

  private

  public :: do_eb_work

contains

  subroutine do_eb_work(lo, hi, cnodes, num,&
       mask, mask_l1, mask_l2, mask_h1, mask_h2) bind(C,name="do_eb_work")

    integer,     intent(in) :: lo(dim), hi(dim)
    integer,     intent(in) :: num
    integer,     intent(in) :: mask_l1,mask_l2,mask_h1,mask_h2
    integer,     intent(in) :: mask(mask_l1:mask_h1,mask_l2:mask_h2)
    type(cnode), intent(in) :: cnodes(0:num-1) 

    integer :: i, j, mC, mL, nC, nL, iL, iC, iLL

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          mC = mask(i  ,j)
          mL = mask(i-1,j)
          if (mC .eq. REGULAR_CELL) then         ! C is regular
             iC = REGULAR_CELL
             if (mL .eq. REGULAR_CELL) then      ! ... and L is regular
                nL = 0
                iL = REGULAR_CELL
             else if (mL .eq. COVERED_CELL) then ! ... and L is covered
                call bl_pd_abort('ERROR: cannot have covered next to regular cell')
             else                                ! ... and L is irregular
                nC = 0
                nL = cnodes(mL)%nCells
                if (nL .ne. 1) then
                   call bl_pd_abort('ERROR: Should not be any mv cut cells here')
                endif
                iL = cnodes(mL)%cells(0)%ebCellID
             endif
          else if (mC .eq. COVERED_CELL) then    ! C is covered
             iC = COVERED_CELL
             iL = mask(i-1,j)
          else                                   ! C is irregular
             if (mC .lt. 0  .or.  mC .ge. num) then
                call bl_pd_abort('ERROR: cut cell index OOB')
             endif
             nC = cnodes(mC)%nCells
             if (nC .gt. 1) then
                call bl_pd_abort('ERROR: Should not be any mv cut cells here')
             endif
             iC = cnodes(mC)%cells(0)%ebCellID
             if (mask(i-1,j) .eq. REGULAR_CELL) then
                iL = REGULAR_CELL
             else if (mask(i-1,j) .eq. COVERED_CELL) then
                iL = COVERED_CELL
             else
                nL = cnodes(mask(i-1,j))%nCells
                if (nL .ge. 0) then                 ! ... and L is irregular
                   if (mask(i-1,j) .lt. 0) then
                      call bl_pd_abort('ERROR: C thinks L is irregular, but it isnt')
                   endif
                   iLL = cnodes(mask(i,j))%cells(0)%nbr(0,0,0)
                   if (iLL .ne. 0) then
                      call bl_pd_abort('ERROR: C thinks his neighbor is a mv cell')
                   endif
                   iL = cnodes(mask(i-1,j))%cells(iLL)%ebCellID
                else if (nL .eq. REGULAR_CELL) then
                   iL = REGULAR_CELL
                else
                   if (mask(i-1,j) .ne. COVERED_CELL) then
                      call bl_pd_abort('ERROR: C thinks L is covered but it isnt')
                   endif
                endif
             endif
          endif
          print *,i,j,'left=',iL,' cent=',iC
       enddo
    enddo

  end subroutine do_eb_work

end module stencil_test_module
