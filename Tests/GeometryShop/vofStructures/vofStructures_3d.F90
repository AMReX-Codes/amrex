#include <Node.H>

module stencil_test_module

  use amrex_fort_module, only : dim=>bl_spacedim
  use amrex_ebstruct_module, only : cnode, fnode

  implicit none

  private

  public :: do_eb_work

contains

  subroutine do_eb_work(lo, hi, cnodes, cnum, fnodesI, fnodesJ, fnodesK, fnum,&
       mask, mask_l1, mask_l2, mask_l3, mask_h1, mask_h2, mask_h3) bind(C,name="do_eb_work")

    integer,     intent(in) :: lo(dim), hi(dim)
    integer,     intent(in) :: cnum, fnum(dim)
    integer,     intent(in) :: mask_l1,mask_l2,mask_l3,mask_h1,mask_h2,mask_h3
    integer,     intent(in) :: mask(mask_l1:mask_h1,mask_l2:mask_h2,mask_l3:mask_h3)
    type(cnode), intent(in) :: cnodes(0:cnum-1) 
    type(fnode), intent(in) :: fnodesI(0:fnum(1)-1)
    type(fnode), intent(in) :: fnodesJ(0:fnum(2)-1)
    type(fnode), intent(in) :: fnodesK(0:fnum(3)-1)

    integer :: i, j, k, mC, mL, nC, nL, iL, iC, iLL

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mC = mask(i  ,j,k)
             mL = mask(i-1,j,k)
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
                iL = mask(i-1,j,k)
             else                                   ! C is irregular
                if (mC .lt. 0  .or.  mC .ge. cnum) then
                   call bl_pd_abort('ERROR: cut cell index OOB')
                endif
                nC = cnodes(mC)%nCells
                if (nC .gt. 1) then
                   call bl_pd_abort('ERROR: Should not be any mv cut cells here')
                endif
                iC = cnodes(mC)%cells(0)%ebCellID
                if (mask(i-1,j,k) .eq. REGULAR_CELL) then
                   iL = REGULAR_CELL
                else if (mask(i-1,j,k) .eq. COVERED_CELL) then
                   iL = COVERED_CELL
                else
                   nL = cnodes(mask(i-1,j,k))%nCells
                   if (nL .ge. 0) then                 ! ... and L is irregular
                      if (mask(i-1,j,k) .lt. 0) then
                         call bl_pd_abort('ERROR: C thinks L is irregular, but it isnt')
                      endif
                      iLL = cnodes(mask(i,j,k))%cells(0)%nbr(0,0,0)
                      if (iLL .ne. 0) then
                         call bl_pd_abort('ERROR: C thinks his neighbor is a mv cell')
                      endif
                      iL = cnodes(mask(i-1,j,k))%cells(iLL)%ebCellID
                   else if (nL .eq. REGULAR_CELL) then
                      iL = REGULAR_CELL
                   else
                      if (mask(i-1,j,k) .ne. COVERED_CELL) then
                         call bl_pd_abort('ERROR: C thinks L is covered but it isnt')
                      endif
                   endif
                endif
             endif
             print *,i,j,k,'left=',iL,' cent=',iC
          enddo
       enddo
    enddo

  end subroutine do_eb_work

end module stencil_test_module
