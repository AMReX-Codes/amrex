
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
