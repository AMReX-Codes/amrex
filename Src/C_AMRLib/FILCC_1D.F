c
c $Id: FILCC_1D.F,v 1.1 2001-08-16 23:06:48 lijewski Exp $
c

#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"

#define SDIM 1

c ::: -----------------------------------------------------------
c ::: This routine is intended to be a generic fill function
c ::: for cell-centered data.  It knows how to extrapolate
c ::: and reflect data and is used to supplement the problem-specific
c ::: fill functions which call it.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: q           <=  array to fill
c ::: lo,hi        => index extent of q array
c ::: domlo,domhi  => index extent of problem domain
c ::: dx           => cell spacing
c ::: xlo          => physical location of lower left hand
c :::	              corner of q array
c ::: bc	   => array of boundary flags bc(SPACEDIM,lo:hi)
c ::: 
c ::: NOTE: all corner as well as edge data is filled if not EXT_DIR
c ::: -----------------------------------------------------------

      subroutine filcc(q,DIMS(q),domlo,domhi,dx,xlo,bc)

      integer    DIMDEC(q)
      integer    domlo(SDIM), domhi(SDIM)
      integer    bc(SDIM,2)
      REAL_T     xlo(SDIM), dx(SDIM)
      REAL_T     q(DIMV(q))

      integer    nlft, nrgt
      integer    ilo, ihi
      integer    i
      integer    is, ie

      nlft = max(0,domlo(1)-ARG_L1(q))
      nrgt = max(0,ARG_H1(q)-domhi(1))

      is = max(ARG_L1(q),domlo(1))
      ie = min(ARG_H1(q),domhi(1))

c     ::::: first fill sides
      if (nlft .gt. 0) then
         ilo = domlo(1)

	 if (bc(1,1) .eq. FOEXTRAP) then
	    do i = 1, nlft
	       q(ilo-i) = q(ilo)
	    end do
	 else if (bc(1,1) .eq. HOEXTRAP) then
	    do i = 2, nlft
	       q(ilo-i) = q(ilo) 
	    end do 
            if (ilo+2 .le. ie) then 
               q(ilo-1) = (fifteen*q(ilo) - ten*q(ilo+1) + 
     $                        three*q(ilo+2)) * eighth
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

	 if (bc(1,2) .eq. FOEXTRAP) then
	    do i = 1, nrgt
	       q(ihi+i) = q(ihi)
	    end do
         else if (bc(1,2) .eq. HOEXTRAP) then
            do i = 2, nrgt
               q(ihi+i) = q(ihi)
            end do
            if (ihi-2 .ge. is) then
	       q(ihi+1) = (fifteen*q(ihi) - ten*q(ihi-1) + 
     $                        three*q(ihi-2)) * eighth
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
      end
