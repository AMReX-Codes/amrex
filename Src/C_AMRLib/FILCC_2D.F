#undef  BL_LANG_CC
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "DIMS.H"

#define SDIM 2

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

      integer    nlft, nrgt, nbot, ntop
      integer    ilo, ihi, jlo, jhi
      integer    i, j
      integer    is, ie, js, je

      nlft = max(0,domlo(1)-ARG_L1(q))
      nrgt = max(0,ARG_H1(q)-domhi(1))
      nbot = max(0,domlo(2)-ARG_L2(q))
      ntop = max(0,ARG_H2(q)-domhi(2))

      is = max(ARG_L1(q),domlo(1))
      ie = min(ARG_H1(q),domhi(1))
      js = max(ARG_L2(q),domlo(2))
      je = min(ARG_H2(q),domhi(2))

c     ::::: first fill sides
      if (nlft .gt. 0) then
         ilo = domlo(1)

	 if (bc(1,1) .eq. EXTRAP) then
	    do i = 1, nlft
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ilo-i,j) = q(ilo,j)
	    enddo
	    enddo
	 elseif (bc(1,1) .eq. HOEXTRAP) then
	    do i = 2, nlft
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ilo-i,j) = q(ilo,j) 
	    enddo 
	    enddo 
            if (ilo+2 .le. ie) then 
	     do j = ARG_L2(q), ARG_H2(q)
		q(ilo-1,j) = (fifteen*q(ilo,j) - ten*q(ilo+1,j) + 
     $                        three*q(ilo+2,j))/eight 
	     enddo 
            else 
	     do j = ARG_L2(q), ARG_H2(q)
	       q(ilo-1,j) = two*q(ilo,j) - q(ilo+1,j) 
	     enddo
            endif
	 elseif (bc(1,1) .eq. REFLECT_EVEN) then
	    do i = 1, nlft
	     do j = ARG_L2(q), ARG_H2(q)
	       q(ilo-i,j) = q(ilo+i-1,j)
	    enddo
	    enddo
	 elseif (bc(1,1) .eq. REFLECT_ODD) then
	    do i = 1, nlft
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ilo-i,j) = -q(ilo+i-1,j)
	    enddo
	    enddo
	 endif
      endif

      if (nrgt .gt. 0) then
         ihi = domhi(1)

	 if (bc(1,2) .eq. EXTRAP) then
	    do i = 1, nrgt
	    do j = ARG_L2(q), ARG_H2(q)
	       q(ihi+i,j) = q(ihi,j)
	    enddo
	    enddo
         elseif (bc(1,2) .eq. HOEXTRAP) then
            do i = 2, nrgt
	    do j = ARG_L2(q), ARG_H2(q)
               q(ihi+i,j) = q(ihi,j)
            enddo
            enddo
            if (ihi-2 .ge. is) then
	     do j = ARG_L2(q), ARG_H2(q)
	       q(ihi+1,j) = (fifteen*q(ihi,j) - ten*q(ihi-1,j) + 
     $                        three*q(ihi-2,j))/eight
             enddo
            else
	     do j = ARG_L2(q), ARG_H2(q)
	       q(ihi+1,j) = two*q(ihi,j) - q(ihi-1,j) 
             enddo
            endif
	 elseif (bc(1,2) .eq. REFLECT_EVEN) then
	    do i = 1, nrgt
            do j = ARG_L2(q), ARG_H2(q)
	       q(ihi+i,j) = q(ihi-i+1,j)
	    enddo
	    enddo
	 elseif (bc(1,2) .eq. REFLECT_ODD) then
	    do i = 1, nrgt
            do j = ARG_L2(q), ARG_H2(q)
	       q(ihi+i,j) = -q(ihi-i+1,j)
	    enddo
	    enddo
	 endif
      endif

      if (nbot .gt. 0) then
         jlo = domlo(2)

	 if (bc(2,1) .eq. EXTRAP) then
	    do j = 1, nbot
	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jlo-j) = q(i,jlo)
	    enddo
	    enddo
         elseif (bc(2,1) .eq. HOEXTRAP) then
            do j = 2, nbot
	    do i = ARG_L1(q), ARG_H1(q)
               q(i,jlo-j) = q(i,jlo)
            enddo
            enddo
            if (jlo+2 .le. je) then
 	     do i = ARG_L1(q), ARG_H1(q)
               q(i,jlo-1) = (fifteen*q(i,jlo) - ten*q(i,jlo+1) + 
     $                        three*q(i,jlo+2))/eight
             enddo
            else
 	     do i = ARG_L1(q), ARG_H1(q)
               q(i,jlo-1) = two*q(i,jlo) - q(i,jlo+1)
             enddo
            endif
	 elseif (bc(2,1) .eq. REFLECT_EVEN) then
	    do j = 1, nbot
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jlo-j) = q(i,jlo+j-1)
	    enddo
	    enddo
	 elseif (bc(2,1) .eq. REFLECT_ODD) then
	    do j = 1, nbot
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jlo-j) = -q(i,jlo+j-1)
	    enddo
	    enddo
	 endif
      endif

      if (ntop .gt. 0) then
         jhi = domhi(2)

	 if (bc(2,2) .eq. EXTRAP) then
	    do j = 1, ntop
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jhi+j) = q(i,jhi)
	    enddo
	    enddo
         elseif (bc(2,2) .eq. HOEXTRAP) then
            do j = 2, ntop
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jhi+j) = q(i,jhi)
            enddo
            enddo
            if (jhi-2 .ge. js) then
 	    do i = ARG_L1(q), ARG_H1(q)
               q(i,jhi+1) = (fifteen*q(i,jhi) - ten*q(i,jhi-1) + 
     $                         three*q(i,jhi-2))/eight
             enddo
            else
 	     do i = ARG_L1(q), ARG_H1(q)
               q(i,jhi+1) = two*q(i,jhi) - q(i,jhi-1)
             enddo
            endif
	 elseif (bc(2,2) .eq. REFLECT_EVEN) then
	    do j = 1, ntop
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jhi+j) = q(i,jhi-j+1)
	    enddo
	    enddo
	 elseif (bc(2,2) .eq. REFLECT_ODD) then
	    do j = 1, ntop
 	    do i = ARG_L1(q), ARG_H1(q)
	       q(i,jhi+j) = -q(i,jhi-j+1)
	    enddo
	    enddo
	 endif
      endif

      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and.
     $    (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) ) then
        q(ilo-1,jlo-1) = half /eight * (
     $    (fifteen*q(ilo-1,jlo) - ten*q(ilo-1,jlo+1) + three*q(ilo-1,jlo+2)) +
     $    (fifteen*q(ilo,jlo-1) - ten*q(ilo+1,jlo-1) + three*q(ilo+2,jlo-1)) )
      endif

      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and.
     $    (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) ) then
        q(ilo-1,jhi+1) = half / eight * (
     $    (fifteen*q(ilo-1,jhi) - ten*q(ilo-1,jhi-1) + three*q(ilo-1,jhi-2)) +
     $    (fifteen*q(ilo,jhi+1) - ten*q(ilo+1,jhi+1) + three*q(ilo+2,jhi+1)) )
      endif

      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and.
     $    (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) ) then
        q(ihi+1,jlo-1) = half / eight * (
     $    (fifteen*q(ihi+1,jlo) - ten*q(ihi+1,jlo+1) + three*q(ihi+1,jlo+2)) +
     $    (fifteen*q(ihi,jlo-1) - ten*q(ihi-1,jlo-1) + three*q(ihi-2,jlo-1)) )
      endif

      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and.
     $    (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) ) then
        q(ihi+1,jhi+1) = half / eight * (
     $    (fifteen*q(ihi+1,jhi) - ten*q(ihi+1,jhi-1) + three*q(ihi+1,jhi-2)) +
     $    (fifteen*q(ihi,jhi+1) - ten*q(ihi-1,jhi+1) + three*q(ihi-2,jhi+1)) )
      endif

      return
      end
