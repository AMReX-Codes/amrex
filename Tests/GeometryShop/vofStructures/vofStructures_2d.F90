#include <Node.H>

module stencil_test_module

  use amrex_fort_module, only : dim=>bl_spacedim
  use amrex_ebstruct_module, only : cnode, fnode

  implicit none

  private

  public :: do_eb_work

contains

  subroutine do_eb_work(lo, hi, cnodes, cnum, fnodesI, fnodesJ, fnum,&
       mask, mask_l1, mask_l2, mask_h1, mask_h2) bind(C,name="do_eb_work")

    integer,     intent(in) :: lo(dim), hi(dim)
    integer,     intent(in) :: cnum, fnum(dim)
    integer,     intent(in) :: mask_l1,mask_l2,mask_h1,mask_h2
    integer,     intent(in) :: mask(mask_l1:mask_h1,mask_l2:mask_h2)
    type(cnode), intent(in) :: cnodes(0:cnum-1)
    type(fnode), intent(in) :: fnodesI(0:fnum(1)-1)
    type(fnode), intent(in) :: fnodesJ(0:fnum(2)-1)

    integer :: i, j, mC, mL, nC, nL, iL, iC, iFL

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
                if (cnodes(mL)%cells(0)%Nnbr(0,1) .ne. 1) then
                   call bl_pd_abort('ERROR: Irreg L does not see reg C neighbor')
                else
                   if (cnodes(mL)%cells(0)%nbr(0,1,0) .ne. REGULAR_CELL) then
                      call bl_pd_abort('ERROR: Irreg L does not see reg C neighbor as regular cell')
                   endif
                endif
                iFL = cnodes(mL)%cells(0)%faceID(0,1,0);
                if (iFL .ne. REGULAR_FACE) then
                   call bl_pd_abort('ERROR: Should not be a cut face between regular and cut cells')
                endif
             endif
          else if (mC .eq. COVERED_CELL) then    ! C is covered
             iC = COVERED_CELL
             iL = mask(i-1,j)
          else                                   ! C is irregular
             if (mC .lt. 0  .or.  mC .ge. cnum) then
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
                if (mL .eq. REGULAR_CELL) then          ! ... and L is regular
                   call bl_pd_abort('should have seen this already')
                else if (mL .eq. COVERED_CELL) then     ! ... and L is covered
                   nL = cnodes(mC)%cells(0)%Nnbr(0,0)
                   if (nL .ne. 0) then
                      call bl_pd_abort('ERROR: C thinks a face connect it to a covered cell')
                   endif
                else                                    ! ... and L is irregular
                   nL = cnodes(mC)%cells(0)%Nnbr(0,0)
                   iL = cnodes(mL)%cells(0)%ebCellID
                   if (nL .ne. 1) then
                      call bl_pd_abort('ERROR: C thinks it has mutliple L neighbors')
                   endif
                   if (cnodes(mC)%cells(0)%nbr(0,0,0) .ne. 0) then
                      call bl_pd_abort('ERROR: C thinks L has index != 0')
                   endif
                   iFL = cnodes(mC)%cells(0)%faceID(0,0,0)
                   if (iFL .lt. 0) then
                      call bl_pd_abort('ERROR: C thinks face between cut cells is not cut')
                   endif
                   if (iFL .ne. cnodes(mL)%cells(0)%faceID(0,1,0)) then
                      call bl_pd_abort('ERROR: inconsistent faceID')
                   endif
                   print *,i,j,'left=',iL,' cent=',iC,'faceID',iFL
                endif
             endif
          endif
       enddo
    enddo

  end subroutine do_eb_work

end module stencil_test_module
