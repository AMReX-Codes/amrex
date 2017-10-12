
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_BC_TYPES.H"

#define SDIM 1

! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for cell-centered data.  It knows how to extrapolate
! ::: and reflect data and is used to supplement the problem-specific
! ::: fill functions which call it.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: q           <=  array to fill
! ::: lo,hi        => index extent of q array
! ::: domlo,domhi  => index extent of problem domain
! ::: dx           => cell spacing
! ::: xlo          => physical location of lower left hand
! :::	              corner of q array
! ::: bc	   => array of boundary flags bc(SPACEDIM,lo:hi)
! ::: 
! ::: NOTE: all corner as well as edge data is filled if not EXT_DIR
! ::: -----------------------------------------------------------

subroutine filcc(q,q_l1,q_h1,domlo,domhi,dx,xlo,bc)

  implicit none

  integer    q_l1, q_h1
  integer    domlo(SDIM), domhi(SDIM)
  integer    bc(SDIM,2)
  REAL_T     xlo(SDIM), dx(SDIM)
  REAL_T     q(q_l1:q_h1)

  integer    nlft, nrgt
  integer    ilo, ihi
  integer    i
  integer    is, ie
  
  nlft = max(0,domlo(1)-q_l1)
  nrgt = max(0,q_h1-domhi(1))

  is = max(q_l1,domlo(1))
  ie = min(q_h1,domhi(1))

  !     ::::: first fill sides
  if (nlft .gt. 0) then
     ilo = domlo(1)

     if (bc(1,1) .eq. EXT_DIR) then
        !     set all ghost cell values to a prescribed dirichlet
        !     value; in this example, we have chosen 1
        !	    do i = 1, nlft
        !	       q(ilo-i) = 1.d0
        !	    end do
     else if (bc(1,1) .eq. FOEXTRAP) then
        do i = 1, nlft
           q(ilo-i) = q(ilo)
        end do
     else if (bc(1,1) .eq. HOEXTRAP) then
        do i = 2, nlft
           q(ilo-i) = q(ilo) 
        end do
        if (ilo+2 .le. ie) then 
           q(ilo-1) = (fifteen*q(ilo) - ten*q(ilo+1) + &
                three*q(ilo+2)) * eighth
        else 
           q(ilo-1) = half*(three*q(ilo) - q(ilo+1))
        end if
     else if (bc(1,1) .eq. REFLECT_EVEN) then
        do i = 1, nlft
           q(ilo-i) = q(ilo+i-1)
        end do
     else if (bc(1,1) .eq. REFLECT_ODD) then
        do i = 1, nlft
           q(ilo-i) = -q(ilo+i-1)
        end do
     end if
  end if

  if (nrgt .gt. 0) then
     ihi = domhi(1)

     if (bc(1,2) .eq. EXT_DIR) then
        !	    do i = 1, nrgt
        !	       q(ihi+i) = 1.d0
        !	    end do
     else if (bc(1,2) .eq. FOEXTRAP) then
        do i = 1, nrgt
           q(ihi+i) = q(ihi)
        end do
     else if (bc(1,2) .eq. HOEXTRAP) then
        do i = 2, nrgt
           q(ihi+i) = q(ihi)
        end do
        if (ihi-2 .ge. is) then
           q(ihi+1) = (fifteen*q(ihi) - ten*q(ihi-1) + &
                three*q(ihi-2)) * eighth
        else
           q(ihi+1) = half*(three*q(ihi) - q(ihi-1))
        end if
     else if (bc(1,2) .eq. REFLECT_EVEN) then
        do i = 1, nrgt
           q(ihi+i) = q(ihi-i+1)
        end do
     else if (bc(1,2) .eq. REFLECT_ODD) then
        do i = 1, nrgt
           q(ihi+i) = -q(ihi-i+1)
        end do
     end if
  end if

  return
end subroutine filcc





module filcc_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  AMREX_DEVICE subroutine filccn(blo, bhi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)

    implicit none

    integer,  intent(in   ) :: blo(3), bhi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: domlo(1), domhi(1)
    integer,  intent(in   ) :: bc(1,2,ncomp)
    real(rt), intent(in   ) :: xlo(1), dx(1)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),ncomp)
    integer,  intent(in   ) :: ncomp

    integer :: ilo, ihi
    integer :: i, n
    integer :: is, ie

    is = max(q_lo(1), domlo(1))
    ie = min(q_hi(1), domhi(1))

    ilo = domlo(1)
    ihi = domhi(1)

    do n = 1, ncomp

       if (bc(1,1,n) .eq. EXT_DIR) then

          ! Do nothing.

       else if (bc(1,1,n) .eq. FOEXTRAP) then

          do i = blo(1), bhi(1)

             if (i < ilo) then
                q(i,n) = q(ilo,n)
             end if

          end do

       else if (bc(1,1,n) .eq. HOEXTRAP) then

          do i = blo(1), bhi(1)

             if (i < ilo - 1) then
                q(i,n) = q(ilo,n)
             else if (i == ilo - 1) then
                if (ilo+2 .le. ie) then
                   q(i,n) = eighth * (fifteen*q(ilo,n) - ten*q(ilo+1,n) + three*q(ilo+2,n))
                else
                   q(i,n) = half * (three*q(ilo,n) - q(ilo+1,n))
                end if
             end if

          end do

       else if (bc(1,1,n) .eq. REFLECT_EVEN) then

          do i = blo(1), bhi(1)

             if (i < ilo) then
                q(i,n) = q(ilo+(ilo-i)-1,n)
             end if

          end do

       else if (bc(1,1,n) .eq. REFLECT_ODD) then

          do i = blo(1), bhi(1)

             if (i < ilo) then
                q(i,n) = -q(ilo+(ilo-i)-1,n)
             end if

          end do

       end if

       if (bc(1,2,n) .eq. EXT_DIR) then

          ! Do nothing.

       else if (bc(1,2,n) .eq. FOEXTRAP) then

          do i = blo(1), bhi(1)

             if (i > ihi) then
                q(i,n) = q(ihi,n)
             end if

          end do

       else if (bc(1,2,n) .eq. HOEXTRAP) then

          do i = blo(1), bhi(1)

             if (i > ihi + 1) then
                q(i,n) = q(ihi,n)
             else if (i == ihi + 1) then
                if (ihi-2 .ge. is) then
                   q(i,n) = eighth * (fifteen*q(ihi,n) - ten*q(ihi-1,n) + three*q(ihi-2,n))
                else
                   q(i,n) = half * (three*q(ihi,n) - q(ihi-1,n))
                end if
             end if

          end do

       else if (bc(1,2,n) .eq. REFLECT_EVEN) then

          do i = blo(1), bhi(1)

             if (i > ihi) then
                q(i,n) = q(ihi-(i-ihi)+1,n)
             end if

          end do

       else if (bc(1,2,n) .eq. REFLECT_ODD) then

          do i = blo(1), bhi(1)

             if (i > ihi) then
                q(i,n) = -q(ihi-(i-ihi)+1,n)
             end if

          end do

       end if

    end do

  end subroutine filccn

end module filcc_module
