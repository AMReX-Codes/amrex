c
c $Id: FILCC_3D.F,v 1.8 2002-07-31 22:32:02 lijewski Exp $
c

#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "ArrayLim.H"

#define SDIM 3

c ::: -----------------------------------------------------------
c ::: This routine is intended to be a generic fill function
c ::: for cell centered data.  It knows how to exrapolate,
c ::: and reflect data and can be used to suppliment problem
c ::: specific fill functions (ie. EXT_DIR).
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: q        <=  array to fill
c ::: DIMS(q)   => index extent of q array
c ::: domlo,hi  => index extent of problem domain
c ::: dx        => cell spacing
c ::: xlo       => physical location of lower left hand
c :::	           corner of q array
c ::: bc	=> array of boundary flags bc(SPACEDIM,lo:hi)
c ::: 
c ::: NOTE: corner data not used in computing soln but must have
c :::       reasonable values for arithmetic to live
c ::: -----------------------------------------------------------

      subroutine filcc(q,DIMS(q),domlo,domhi,dx,xlo,bc)

      integer    DIMDEC(q)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     xlo(SDIM), dx(SDIM)
      REAL_T     q(DIMV(q))
      integer    bc(SDIM,2)

      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    is,  ie,  js,  je,  ks,  ke
      integer    i, j, k

      is = max(ARG_L1(q),domlo(1))
      ie = min(ARG_H1(q),domhi(1))
      js = max(ARG_L2(q),domlo(2))
      je = min(ARG_H2(q),domhi(2))
      ks = max(ARG_L3(q),domlo(3))
      ke = min(ARG_H3(q),domhi(3))

      nlft = max(0,domlo(1)-ARG_L1(q))
      nrgt = max(0,ARG_H1(q)-domhi(1))
      nbot = max(0,domlo(2)-ARG_L2(q))
      ntop = max(0,ARG_H2(q)-domhi(2))
      ndwn = max(0,domlo(3)-ARG_L3(q))
      nup  = max(0,ARG_H3(q)-domhi(3))
c
c     ::::: first fill sides
c
      if (nlft .gt. 0) then
         ilo = domlo(1)

	 if (bc(1,1) .eq. FOEXTRAP) then
	    do i = 1, nlft
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ilo-i,j,k) = q(ilo,j,k)
                  end do
               end do
	    end do
	 else if (bc(1,1) .eq. HOEXTRAP) then
	    do i = 2, nlft
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ilo-i,j,k) = q(ilo,j,k)
                  end do
               end do
	    end do
            if (ilo+2 .le. ie) then
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ilo-1,j,k) = (15*q(ilo,j,k) - 10*q(ilo+1,j,k) + 
     $                    3*q(ilo+2,j,k)) * eighth
                  end do
               end do
            else  
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ilo-1,j,k) = half*(3*q(ilo,j,k) - q(ilo+1,j,k))
                  end do
               end do
            end if
	 else if (bc(1,1) .eq. REFLECT_EVEN) then
	    do i = 1, nlft
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ilo-i,j,k) = q(ilo+i-1,j,k)
                  end do
               end do
	    end do
	 else if (bc(1,1) .eq. REFLECT_ODD) then
	    do i = 1, nlft
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ilo-i,j,k) = -q(ilo+i-1,j,k)
                  end do
               end do
	    end do
	 end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)

	 if (bc(1,2) .eq. FOEXTRAP) then
	    do i = 1, nrgt
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ihi+i,j,k) = q(ihi,j,k)
                  end do
               end do
	    end do
         else if (bc(1,2) .eq. HOEXTRAP) then
            do i = 2, nrgt
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ihi+i,j,k) = q(ihi,j,k)
                  end do
               end do
            end do
            if (ihi-2 .ge. is) then
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ihi+1,j,k) = (15*q(ihi,j,k) - 10*q(ihi-1,j,k) + 
     $                    3*q(ihi-2,j,k)) * eighth
                  end do
               end do
            else
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ihi+1,j,k) = half*(3*q(ihi,j,k) - q(ihi-1,j,k))
                  end do
               end do
            end if
	 else if (bc(1,2) .eq. REFLECT_EVEN) then
	    do i = 1, nrgt
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ihi+i,j,k) = q(ihi-i+1,j,k)
                  end do
               end do
	    end do
	 else if (bc(1,2) .eq. REFLECT_ODD) then
	    do i = 1, nrgt
               do k = ARG_L3(q),ARG_H3(q)
                  do j = ARG_L2(q),ARG_H2(q)
                     q(ihi+i,j,k) = -q(ihi-i+1,j,k)
                  end do
               end do
	    end do
	 end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)
         
	 if (bc(2,1) .eq. FOEXTRAP) then
	    do j = 1, nbot
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jlo-j,k) = q(i,jlo,k)
                  end do
               end do
	    end do
         else if (bc(2,1) .eq. HOEXTRAP) then
            do j = 2, nbot
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jlo-j,k) = q(i,jlo,k)
                  end do
               end do
            end do
            if (jlo+2 .le. je) then
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jlo-1,k) = (15*q(i,jlo,k) - 10*q(i,jlo+1,k) + 
     $                    3*q(i,jlo+2,k)) * eighth
                  end do
               end do
	    else
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jlo-1,k) = half*(3*q(i,jlo,k) - q(i,jlo+1,k))
                  end do
               end do
	    end if
	 else if (bc(2,1) .eq. REFLECT_EVEN) then
	    do j = 1, nbot 
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jlo-j,k) = q(i,jlo+j-1,k)
                  end do
               end do
	    end do
	 else if (bc(2,1) .eq. REFLECT_ODD) then
	    do j = 1, nbot
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jlo-j,k) = -q(i,jlo+j-1,k)
                  end do
               end do
	    end do
	 end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)

	 if (bc(2,2) .eq. FOEXTRAP) then
	    do j = 1, ntop
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jhi+j,k) = q(i,jhi,k)
                  end do
               end do
	    end do
         else if (bc(2,2) .eq. HOEXTRAP) then
            do j = 2, ntop
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jhi+j,k) = q(i,jhi,k)
                  end do
               end do
            end do
            if (jhi-2 .ge. js) then
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jhi+1,k) = (15*q(i,jhi,k) - 10*q(i,jhi-1,k) + 
     $                    3*q(i,jhi-2,k)) * eighth
                  end do
               end do
	    else
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jhi+1,k) = half*(3*q(i,jhi,k) - q(i,jhi-1,k))
                  end do
               end do
	    end if
	 else if (bc(2,2) .eq. REFLECT_EVEN) then
	    do j = 1, ntop
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jhi+j,k) = q(i,jhi-j+1,k)
                  end do
               end do
	    end do
	 else if (bc(2,2) .eq. REFLECT_ODD) then
	    do j = 1, ntop
               do k = ARG_L3(q),ARG_H3(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,jhi+j,k) = -q(i,jhi-j+1,k)
                  end do
               end do
	    end do
	 end if
      end if

      if (ndwn .gt. 0) then
         klo = domlo(3)

	 if (bc(3,1) .eq. FOEXTRAP) then
	    do k = 1, ndwn
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,klo-k) = q(i,j,klo)
                  end do
               end do
	    end do
         else if (bc(3,1) .eq. HOEXTRAP) then
            do k = 2, ndwn
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,klo-k) = q(i,j,klo)
                  end do
               end do
            end do
            if (klo+2 .le. ke) then
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,klo-1) = (15*q(i,j,klo) - 10*q(i,j,klo+1) + 
     $                    3*q(i,j,klo+2)) * eighth
                  end do
               end do
	    else
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,klo-1) = half*(3*q(i,j,klo) - q(i,j,klo+1))
                  end do
               end do
	    end if
	 else if (bc(3,1) .eq. REFLECT_EVEN) then
	    do k = 1, ndwn
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,klo-k) = q(i,j,klo+k-1)
                  end do
               end do
	    end do
	 else if (bc(3,1) .eq. REFLECT_ODD) then
	    do k = 1, ndwn
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,klo-k) = -q(i,j,klo+k-1)
                  end do
               end do
	    end do
	 end if
      end if

      if (nup .gt. 0) then
         khi = domhi(3)

	 if (bc(3,2) .eq. FOEXTRAP) then
	    do k = 1, nup
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,khi+k) = q(i,j,khi)
                  end do
               end do
	    end do
         else if (bc(3,2) .eq. HOEXTRAP) then
            do k = 2, nup
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,khi+k) = q(i,j,khi)
                  end do
               end do
            end do
            if (khi-2 .ge. ks) then
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,khi+1) = (15*q(i,j,khi) - 10*q(i,j,khi-1) + 
     $                    3*q(i,j,khi-2)) * eighth
                  end do
               end do
	    else
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,khi+1) = half*(3*q(i,j,khi) - q(i,j,khi-1))
                  end do
               end do
	    end if
	 else if (bc(3,2) .eq. REFLECT_EVEN) then
	    do k = 1, nup
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,khi+k) = q(i,j,khi-k+1)
                  end do
               end do
	    end do
	 else if (bc(3,2) .eq. REFLECT_ODD) then
	    do k = 1, nup
               do j = ARG_L2(q),ARG_H2(q)
                  do i = ARG_L1(q),ARG_H1(q)
                     q(i,j,khi+k) = -q(i,j,khi-k+1)
                  end do
               end do
	    end do
	 end if
      end if
c
c    First correct the i-j edges and all corners
c
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and.
     $     (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) ) then
         if (jlo+2 .le. je) then
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jlo-1,k) = half * eighth * 
     $              (15*q(ilo-1,jlo,k) - 10*q(ilo-1,jlo+1,k) + 
     $              3*q(ilo-1,jlo+2,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jlo-1,k) = half * half * 
     $              (3*q(ilo-1,jlo,k) - q(ilo-1,jlo+1,k))
            end do
         end if

         if (ilo+2 .le. ie) then
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + half * eighth * 
     $              (15*q(ilo,jlo-1,k) - 10*q(ilo+1,jlo-1,k) + 
     $              3*q(ilo+2,jlo-1,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + half * half * 
     $              (3*q(ilo,jlo-1,k) - q(ilo+1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
            if (klo+2 .le. ke) then
               q(ilo-1,jlo-1,klo-1) = eighth * (
     $              (15*q(ilo-1,jlo-1,klo) - 10*q(ilo-1,jlo-1,klo+1) +
     $              3*q(ilo-1,jlo-1,klo+2)) )
            else
               q(ilo-1,jlo-1,klo-1) = half * 
     $              (3*q(ilo-1,jlo-1,klo) - q(ilo-1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
            if (khi-2 .ge. ks) then
               q(ilo-1,jlo-1,khi+1) = eighth * (
     $              (15*q(ilo-1,jlo-1,khi) - 10*q(ilo-1,jlo-1,khi-1) +
     $              3*q(ilo-1,jlo-1,khi-2)) )
            else
               q(ilo-1,jlo-1,khi+1) = half * 
     $              (3*q(ilo-1,jlo-1,khi) - q(ilo-1,jlo-1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and.
     $     (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) ) then
         if (jhi-2 .ge. js) then 
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jhi+1,k) = half * eighth * 
     $              (15*q(ilo-1,jhi,k) - 10*q(ilo-1,jhi-1,k) + 
     $              3*q(ilo-1,jhi-2,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jhi+1,k) = half * half * 
     $              (3*q(ilo-1,jhi,k) - q(ilo-1,jhi-1,k))
            end do
         end if

         if (ilo+2 .le. ie) then 
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + half * eighth * 
     $              (15*q(ilo,jhi+1,k) - 10*q(ilo+1,jhi+1,k) + 
     $              3*q(ilo+2,jhi+1,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + half * half * 
     $              (3*q(ilo,jhi+1,k) - q(ilo+1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
            if (klo+2 .le. ke) then
               q(ilo-1,jhi+1,klo-1) = eighth * (
     $              (15*q(ilo-1,jhi+1,klo) - 10*q(ilo-1,jhi+1,klo+1) +
     $              3*q(ilo-1,jhi+1,klo+2)) )
            else
               q(ilo-1,jhi+1,klo-1) = half * 
     $              (3*q(ilo-1,jhi+1,klo) - q(ilo-1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
            if (khi-2 .ge. ks) then
               q(ilo-1,jhi+1,khi+1) = eighth * (
     $              (15*q(ilo-1,jhi+1,khi) - 10*q(ilo-1,jhi+1,khi-1) +
     $              3*q(ilo-1,jhi+1,khi-2)) )
            else
               q(ilo-1,jhi+1,khi+1) = half * 
     $              (3*q(ilo-1,jhi+1,khi) - q(ilo-1,jhi+1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and.
     $     (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) ) then
         if (jlo+2 .le. je) then 
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jlo-1,k) = half * eighth * 
     $              (15*q(ihi+1,jlo,k) - 10*q(ihi+1,jlo+1,k) + 
     $              3*q(ihi+1,jlo+2,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jlo-1,k) = half * half * 
     $              (3*q(ihi+1,jlo,k) - q(ihi+1,jlo+1,k))
            end do
         end if

         if (ihi-2 .ge. is) then 
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + half * eighth * 
     $              (15*q(ihi,jlo-1,k) - 10*q(ihi-1,jlo-1,k) + 
     $              3*q(ihi-2,jlo-1,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + half * half * 
     $              (3*q(ihi,jlo-1,k) - q(ihi-1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
            if (klo+2 .le. ke) then
               q(ihi+1,jlo-1,klo-1) = eighth * 
     $              (15*q(ihi+1,jlo-1,klo) - 10*q(ihi+1,jlo-1,klo+1) +
     $              3*q(ihi+1,jlo-1,klo+2))
            else
               q(ihi+1,jlo-1,klo-1) = half * 
     $              (3*q(ihi+1,jlo-1,klo) - q(ihi+1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
            if (khi-2 .ge. ks) then
               q(ihi+1,jlo-1,khi+1) = eighth * 
     $              (15*q(ihi+1,jlo-1,khi) - 10*q(ihi+1,jlo-1,khi-1) +
     $              3*q(ihi+1,jlo-1,khi-2))
            else
               q(ihi+1,jlo-1,khi+1) = half * 
     $              (3*q(ihi+1,jlo-1,khi) - q(ihi+1,jlo-1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and.
     $     (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) ) then
         if (jhi-2 .ge. js) then
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jhi+1,k) = half * eighth * 
     $              (15*q(ihi+1,jhi,k) - 10*q(ihi+1,jhi-1,k) + 
     $              3*q(ihi+1,jhi-2,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jhi+1,k) = half * half * 
     $              (3*q(ihi+1,jhi,k) - q(ihi+1,jhi-1,k))
            end do
         end if

         if (ihi-2 .ge. is) then
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + half * eighth * 
     $              (15*q(ihi,jhi+1,k) - 10*q(ihi-1,jhi+1,k) + 
     $              3*q(ihi-2,jhi+1,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + half * half * 
     $              (3*q(ihi,jhi+1,k) - q(ihi-1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
            if (klo+2 .le. ke) then
               q(ihi+1,jhi+1,klo-1) = eighth *
     $              (15*q(ihi+1,jhi+1,klo) - 10*q(ihi+1,jhi+1,klo+1) +
     $              3*q(ihi+1,jhi+1,klo+2))
            else
               q(ihi+1,jhi+1,klo-1) = half * 
     $              (3*q(ihi+1,jhi+1,klo) - q(ihi+1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
            if (khi-2 .ge. ks) then
               q(ihi+1,jhi+1,khi+1) = eighth *
     $              (15*q(ihi+1,jhi+1,khi) - 10*q(ihi+1,jhi+1,khi-1) +
     $              3*q(ihi+1,jhi+1,khi-2))
            else
               q(ihi+1,jhi+1,khi+1) = half * 
     $              (3*q(ihi+1,jhi+1,khi) - q(ihi+1,jhi+1,khi-1))
            end if
         end if

      end if
c
c    Next correct the i-k edges
c
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and.
     $     (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,klo-1) = half * eighth * 
     $              (15*q(ilo-1,j,klo) - 10*q(ilo-1,j,klo+1) + 
     $              3*q(ilo-1,j,klo+2))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,klo-1) = half * half * 
     $              (3*q(ilo-1,j,klo) - q(ilo-1,j,klo+1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + half * eighth * 
     $              (15*q(ilo,j,klo-1) - 10*q(ilo+1,j,klo-1) + 
     $              3*q(ilo+2,j,klo-1))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + half * half * 
     $              (3*q(ilo,j,klo-1) - q(ilo+1,j,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and.
     $     (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,khi+1) = half * eighth * 
     $              (15*q(ilo-1,j,khi) - 10*q(ilo-1,j,khi-1) + 
     $              3*q(ilo-1,j,khi-2))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,khi+1) = half * half * 
     $              (3*q(ilo-1,j,khi) - q(ilo-1,j,khi-1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + half * eighth * 
     $              (15*q(ilo,j,khi+1) - 10*q(ilo+1,j,khi+1) + 
     $              3*q(ilo+2,j,khi+1))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + half * half * 
     $              (3*q(ilo,j,khi+1) - q(ilo+1,j,khi+1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and.
     $     (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,klo-1) = half * eighth *
     $              (15*q(ihi+1,j,klo) - 10*q(ihi+1,j,klo+1) + 
     $              3*q(ihi+1,j,klo+2))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,klo-1) = half * half * 
     $              (3*q(ihi+1,j,klo) - q(ihi+1,j,klo+1))
            end do
         end if

         if (ihi-2 .ge. is) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + half * eighth *
     $              (15*q(ihi,j,klo-1) - 10*q(ihi-1,j,klo-1) + 
     $              3*q(ihi-2,j,klo-1))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + half * half * 
     $              (3*q(ihi,j,klo-1) - q(ihi-1,j,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and.
     $     (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,khi+1) = half * eighth * 
     $              (15*q(ihi+1,j,khi) - 10*q(ihi+1,j,khi-1) + 
     $              3*q(ihi+1,j,khi-2))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,khi+1) = half * half * 
     $              (3*q(ihi+1,j,khi) - q(ihi+1,j,khi-1))
            end do
         end if
         if (ihi-2 .ge. is) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + half * eighth * 
     $              (15*q(ihi,j,khi+1) - 10*q(ihi-1,j,khi+1) + 
     $              3*q(ihi-2,j,khi+1))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + half * half * 
     $              (3*q(ihi,j,khi+1) - q(ihi-1,j,khi+1))
            end do
         end if
      end if
c
c    Next correct the j-k edges
c
      if ((nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) .and.
     $     (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,klo-1) = half * eighth *
     $              (15*q(i,jlo-1,klo) - 10*q(i,jlo-1,klo+1) + 
     $              3*q(i,jlo-1,klo+2))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,klo-1) = half * half * 
     $              (3*q(i,jlo-1,klo) - q(i,jlo-1,klo+1))
            end do
         end if
         if (jlo+2 .le. je) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + half * eighth * 
     $              (15*q(i,jlo,klo-1) - 10*q(i,jlo+1,klo-1) + 
     $              3*q(i,jlo+2,klo-1))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + half * half * 
     $              (3*q(i,jlo,klo-1) - q(i,jlo+1,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) .and.
     $     (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,khi+1) = half * eighth * 
     $              (15*q(i,jlo-1,khi) - 10*q(i,jlo-1,khi-1) + 
     $              3*q(i,jlo-1,khi-2))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,khi+1) = half * half *
     $              (3*q(i,jlo-1,khi) - q(i,jlo-1,khi-1))
            end do
         end if

         if (jlo+2 .le. je) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + half * eighth * 
     $              (15*q(i,jlo,khi+1) - 10*q(i,jlo+1,khi+1) + 
     $              3*q(i,jlo+2,khi+1))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + half * half *
     $              (3*q(i,jlo,khi+1) - q(i,jlo+1,khi+1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) .and.
     $     (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
         if (klo+2 .le. ke) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,klo-1) = half * eighth * 
     $              (15*q(i,jhi+1,klo) - 10*q(i,jhi+1,klo+1) + 
     $              3*q(i,jhi+1,klo+2))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,klo-1) = half * half * 
     $              (3*q(i,jhi+1,klo) - q(i,jhi+1,klo+1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + half * eighth * 
     $              (15*q(i,jhi,klo-1) - 10*q(i,jhi-1,klo-1) + 
     $              3*q(i,jhi-2,klo-1))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + half * half * 
     $              (3*q(i,jhi,klo-1) - q(i,jhi-1,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) .and.
     $     (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
         if (khi-2 .ge. ks) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,khi+1) = half * eighth * 
     $              (15*q(i,jhi+1,khi) - 10*q(i,jhi+1,khi-1) + 
     $              3*q(i,jhi+1,khi-2))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,khi+1) = half * half * 
     $              (3*q(i,jhi+1,khi) - q(i,jhi+1,khi-1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + half * eighth * 
     $              (15*q(i,jhi,khi+1) - 10*q(i,jhi-1,khi+1) + 
     $              3*q(i,jhi-2,khi+1))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + half * half * 
     $              (3*q(i,jhi,khi+1) - q(i,jhi-1,khi+1))
            end do
         end if
      end if

      end

      subroutine hoextraptocc(q,DIMS(q),domlo,domhi,dx,xlo)

      integer    DIMDEC(q)
      integer    domlo(SDIM), domhi(SDIM)
      REAL_T     xlo(SDIM), dx(SDIM)
      REAL_T     q(DIMV(q))

      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    is,  ie,  js,  je,  ks,  ke
      integer    i, j, k

      is = max(ARG_L1(q),domlo(1))
      ie = min(ARG_H1(q),domhi(1))
      js = max(ARG_L2(q),domlo(2))
      je = min(ARG_H2(q),domhi(2))
      ks = max(ARG_L3(q),domlo(3))
      ke = min(ARG_H3(q),domhi(3))

      nlft = max(0,domlo(1)-ARG_L1(q))
      nrgt = max(0,ARG_H1(q)-domhi(1))
      nbot = max(0,domlo(2)-ARG_L2(q))
      ntop = max(0,ARG_H2(q)-domhi(2))
      ndwn = max(0,domlo(3)-ARG_L3(q))
      nup  = max(0,ARG_H3(q)-domhi(3))
c
c     First fill sides.
c
      if (nlft .gt. 0) then
         ilo = domlo(1)

         do i = 2, nlft
            do k = ARG_L3(q),ARG_H3(q)
               do j = ARG_L2(q),ARG_H2(q)
                  q(ilo-i,j,k) = q(ilo,j,k)
               end do
            end do
         end do
         if (ilo+2 .le. ie) then
            do k = ARG_L3(q),ARG_H3(q)
               do j = ARG_L2(q),ARG_H2(q)
                  q(ilo-1,j,k) = 3*q(ilo,j,k) - 3*q(ilo+1,j,k) +
     $                           q(ilo+2,j,k)
               end do
            end do
         else  
            do k = ARG_L3(q),ARG_H3(q)
               do j = ARG_L2(q),ARG_H2(q)
                  q(ilo-1,j,k) = 2*q(ilo,j,k) - q(ilo+1,j,k)
               end do
            end do
         end if
      end if

      if (nrgt .gt. 0) then
         ihi = domhi(1)

         do i = 2, nrgt
            do k = ARG_L3(q),ARG_H3(q)
               do j = ARG_L2(q),ARG_H2(q)
                  q(ihi+i,j,k) = q(ihi,j,k)
               end do
            end do
         end do
         if (ihi-2 .ge. is) then
            do k = ARG_L3(q),ARG_H3(q)
               do j = ARG_L2(q),ARG_H2(q)
                  q(ihi+1,j,k) = 3*q(ihi,j,k) - 3*q(ihi-1,j,k) +
     $                           q(ihi-2,j,k)
               end do
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               do j = ARG_L2(q),ARG_H2(q)
                  q(ihi+1,j,k) = 2*q(ihi,j,k) - q(ihi-1,j,k)
               end do
            end do
         end if
      end if

      if (nbot .gt. 0) then
         jlo = domlo(2)
         
         do j = 2, nbot
            do k = ARG_L3(q),ARG_H3(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,jlo-j,k) = q(i,jlo,k)
               end do
            end do
         end do
         if (jlo+2 .le. je) then
            do k = ARG_L3(q),ARG_H3(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,jlo-1,k) = 3*q(i,jlo,k) - 3*q(i,jlo+1,k) +
     $                           q(i,jlo+2,k)
               end do
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,jlo-1,k) = 2*q(i,jlo,k) - q(i,jlo+1,k)
               end do
            end do
         end if
      end if

      if (ntop .gt. 0) then
         jhi = domhi(2)

         do j = 2, ntop
            do k = ARG_L3(q),ARG_H3(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,jhi+j,k) = q(i,jhi,k)
               end do
            end do
         end do
         if (jhi-2 .ge. js) then
            do k = ARG_L3(q),ARG_H3(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,jhi+1,k) = 3*q(i,jhi,k) - 3*q(i,jhi-1,k) +
     $                           q(i,jhi-2,k)
               end do
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,jhi+1,k) = 2*q(i,jhi,k) - q(i,jhi-1,k)
               end do
            end do
         end if
      end if

      if (ndwn .gt. 0) then
         klo = domlo(3)

         do k = 2, ndwn
            do j = ARG_L2(q),ARG_H2(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,j,klo-k) = q(i,j,klo)
               end do
            end do
         end do
         if (klo+2 .le. ke) then
            do j = ARG_L2(q),ARG_H2(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,j,klo-1) = 3*q(i,j,klo) - 3*q(i,j,klo+1) +
     $                           q(i,j,klo+2)
               end do
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,j,klo-1) = 2*q(i,j,klo) - q(i,j,klo+1)
               end do
            end do
         end if
      end if

      if (nup .gt. 0) then
         khi = domhi(3)

         do k = 2, nup
            do j = ARG_L2(q),ARG_H2(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,j,khi+k) = q(i,j,khi)
               end do
            end do
         end do
         if (khi-2 .ge. ks) then
            do j = ARG_L2(q),ARG_H2(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,j,khi+1) = 3*q(i,j,khi) - 3*q(i,j,khi-1) +
     $                           q(i,j,khi-2)
               end do
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               do i = ARG_L1(q),ARG_H1(q)
                  q(i,j,khi+1) = 2*q(i,j,khi) - q(i,j,khi-1)
               end do
            end do
         end if
      end if
c
c    First correct the i-j edges and all corners
c
      if ((nlft .gt. 0) .and. (nbot .gt. 0)) then
         if (jlo+2 .le. je) then
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jlo-1,k) = half *
     $              (3*q(ilo-1,jlo,k) - 3*q(ilo-1,jlo+1,k) +
     $              q(ilo-1,jlo+2,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jlo-1,k) = half *
     $              (2*q(ilo-1,jlo,k) - q(ilo-1,jlo+1,k))
            end do
         end if

         if (ilo+2 .le. ie) then
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + half *
     $              (3*q(ilo,jlo-1,k) - 3*q(ilo+1,jlo-1,k) + 
     $              q(ilo+2,jlo-1,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + half *
     $              (2*q(ilo,jlo-1,k) - q(ilo+1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               q(ilo-1,jlo-1,klo-1) = 
     $              (3*q(ilo-1,jlo-1,klo) - 3*q(ilo-1,jlo-1,klo+1) +
     $              q(ilo-1,jlo-1,klo+2))
            else
               q(ilo-1,jlo-1,klo-1) = 
     $              (2*q(ilo-1,jlo-1,klo) - q(ilo-1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               q(ilo-1,jlo-1,khi+1) = 
     $              (3*q(ilo-1,jlo-1,khi) - 3*q(ilo-1,jlo-1,khi-1) +
     $              q(ilo-1,jlo-1,khi-2))
            else
               q(ilo-1,jlo-1,khi+1) =
     $              (2*q(ilo-1,jlo-1,khi) - q(ilo-1,jlo-1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nlft .gt. 0) .and. (ntop .gt. 0)) then
         if (jhi-2 .ge. js) then 
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jhi+1,k) = half *
     $              (3*q(ilo-1,jhi,k) - 3*q(ilo-1,jhi-1,k) + 
     $              q(ilo-1,jhi-2,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jhi+1,k) = half *
     $              (2*q(ilo-1,jhi,k) - q(ilo-1,jhi-1,k))
            end do
         end if

         if (ilo+2 .le. ie) then 
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + half *
     $              (3*q(ilo,jhi+1,k) - 3*q(ilo+1,jhi+1,k) + 
     $              q(ilo+2,jhi+1,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + half *
     $              (2*q(ilo,jhi+1,k) - q(ilo+1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               q(ilo-1,jhi+1,klo-1) = 
     $              (3*q(ilo-1,jhi+1,klo) - 3*q(ilo-1,jhi+1,klo+1) +
     $              q(ilo-1,jhi+1,klo+2))
            else
               q(ilo-1,jhi+1,klo-1) =
     $              (2*q(ilo-1,jhi+1,klo) - q(ilo-1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               q(ilo-1,jhi+1,khi+1) =
     $              (3*q(ilo-1,jhi+1,khi) - 3*q(ilo-1,jhi+1,khi-1) +
     $              q(ilo-1,jhi+1,khi-2))
            else
               q(ilo-1,jhi+1,khi+1) =
     $              (2*q(ilo-1,jhi+1,khi) - q(ilo-1,jhi+1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0) .and. (nbot .gt. 0)) then
         if (jlo+2 .le. je) then 
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jlo-1,k) = half *
     $              (3*q(ihi+1,jlo,k) - 3*q(ihi+1,jlo+1,k) + 
     $              q(ihi+1,jlo+2,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jlo-1,k) = half *
     $              (2*q(ihi+1,jlo,k) - q(ihi+1,jlo+1,k))
            end do
         end if

         if (ihi-2 .ge. is) then 
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + half *
     $              (3*q(ihi,jlo-1,k) - 3*q(ihi-1,jlo-1,k) + 
     $              q(ihi-2,jlo-1,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + half *
     $              (2*q(ihi,jlo-1,k) - q(ihi-1,jlo-1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               q(ihi+1,jlo-1,klo-1) =
     $              (3*q(ihi+1,jlo-1,klo) - 3*q(ihi+1,jlo-1,klo+1) +
     $              q(ihi+1,jlo-1,klo+2))
            else
               q(ihi+1,jlo-1,klo-1) =
     $              (2*q(ihi+1,jlo-1,klo) - q(ihi+1,jlo-1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               q(ihi+1,jlo-1,khi+1) =
     $              (3*q(ihi+1,jlo-1,khi) - 3*q(ihi+1,jlo-1,khi-1) +
     $              q(ihi+1,jlo-1,khi-2))
            else
               q(ihi+1,jlo-1,khi+1) =
     $              (2*q(ihi+1,jlo-1,khi) - q(ihi+1,jlo-1,khi-1))
            end if
         end if

      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0) .and. (ntop .gt. 0)) then
         if (jhi-2 .ge. js) then
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jhi+1,k) = half *
     $              (3*q(ihi+1,jhi,k) - 3*q(ihi+1,jhi-1,k) + 
     $              q(ihi+1,jhi-2,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jhi+1,k) = half *
     $              (2*q(ihi+1,jhi,k) - q(ihi+1,jhi-1,k))
            end do
         end if

         if (ihi-2 .ge. is) then
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + half *
     $              (3*q(ihi,jhi+1,k) - 3*q(ihi-1,jhi+1,k) + 
     $              q(ihi-2,jhi+1,k))
            end do
         else
            do k = ARG_L3(q),ARG_H3(q)
               q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + half *
     $              (2*q(ihi,jhi+1,k) - q(ihi-1,jhi+1,k))
            end do
         end if

         if (ndwn .gt. 0) then
            if (klo+2 .le. ke) then
               q(ihi+1,jhi+1,klo-1) =
     $              (3*q(ihi+1,jhi+1,klo) - 3*q(ihi+1,jhi+1,klo+1) +
     $              q(ihi+1,jhi+1,klo+2))
            else
               q(ihi+1,jhi+1,klo-1) =
     $              (2*q(ihi+1,jhi+1,klo) - q(ihi+1,jhi+1,klo+1))
            end if
         end if

         if (nup .gt. 0) then
            if (khi-2 .ge. ks) then
               q(ihi+1,jhi+1,khi+1) =
     $              (3*q(ihi+1,jhi+1,khi) - 3*q(ihi+1,jhi+1,khi-1) +
     $              q(ihi+1,jhi+1,khi-2))
            else
               q(ihi+1,jhi+1,khi+1) =
     $              (2*q(ihi+1,jhi+1,khi) - q(ihi+1,jhi+1,khi-1))
            end if
         end if

      end if
c
c    Next correct the i-k edges
c
      if ((nlft .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,klo-1) = half *
     $              (3*q(ilo-1,j,klo) - 3*q(ilo-1,j,klo+1) + 
     $              q(ilo-1,j,klo+2))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,klo-1) = half *
     $              (2*q(ilo-1,j,klo) - q(ilo-1,j,klo+1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + half *
     $              (3*q(ilo,j,klo-1) - 3*q(ilo+1,j,klo-1) + 
     $              q(ilo+2,j,klo-1))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + half *
     $              (2*q(ilo,j,klo-1) - q(ilo+1,j,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nlft .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,khi+1) = half *
     $              (3*q(ilo-1,j,khi) - 3*q(ilo-1,j,khi-1) + 
     $              q(ilo-1,j,khi-2))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,khi+1) = half *
     $              (2*q(ilo-1,j,khi) - q(ilo-1,j,khi-1))
            end do
         end if

         if (ilo+2 .le. ie) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + half *
     $              (3*q(ilo,j,khi+1) - 3*q(ilo+1,j,khi+1) + 
     $              q(ilo+2,j,khi+1))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + half *
     $              (2*q(ilo,j,khi+1) - q(ilo+1,j,khi+1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,klo-1) = half *
     $              (3*q(ihi+1,j,klo) - 3*q(ihi+1,j,klo+1) + 
     $              q(ihi+1,j,klo+2))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,klo-1) = half *
     $              (2*q(ihi+1,j,klo) - q(ihi+1,j,klo+1))
            end do
         end if

         if (ihi-2 .ge. is) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + half *
     $              (3*q(ihi,j,klo-1) - 3*q(ihi-1,j,klo-1) + 
     $              q(ihi-2,j,klo-1))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + half *
     $              (2*q(ihi,j,klo-1) - q(ihi-1,j,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nrgt .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,khi+1) = half *
     $              (3*q(ihi+1,j,khi) - 3*q(ihi+1,j,khi-1) + 
     $              q(ihi+1,j,khi-2))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,khi+1) = half *
     $              (2*q(ihi+1,j,khi) - q(ihi+1,j,khi-1))
            end do
         end if
         if (ihi-2 .ge. is) then
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + half *
     $              (3*q(ihi,j,khi+1) - 3*q(ihi-1,j,khi+1) + 
     $              q(ihi-2,j,khi+1))
            end do
         else
            do j = ARG_L2(q),ARG_H2(q)
               q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + half *
     $              (2*q(ihi,j,khi+1) - q(ihi-1,j,khi+1))
            end do
         end if
      end if
c
c    Next correct the j-k edges
c
      if ((nbot .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,klo-1) = half *
     $              (3*q(i,jlo-1,klo) - 3*q(i,jlo-1,klo+1) + 
     $              q(i,jlo-1,klo+2))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,klo-1) = half *
     $              (2*q(i,jlo-1,klo) - q(i,jlo-1,klo+1))
            end do
         end if
         if (jlo+2 .le. je) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + half *
     $              (3*q(i,jlo,klo-1) - 3*q(i,jlo+1,klo-1) + 
     $              q(i,jlo+2,klo-1))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + half *
     $              (2*q(i,jlo,klo-1) - q(i,jlo+1,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((nbot .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,khi+1) = half *
     $              (3*q(i,jlo-1,khi) - 3*q(i,jlo-1,khi-1) + 
     $              q(i,jlo-1,khi-2))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,khi+1) = half *
     $              (2*q(i,jlo-1,khi) - q(i,jlo-1,khi-1))
            end do
         end if

         if (jlo+2 .le. je) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + half *
     $              (3*q(i,jlo,khi+1) - 3*q(i,jlo+1,khi+1) + 
     $              q(i,jlo+2,khi+1))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + half *
     $              (2*q(i,jlo,khi+1) - q(i,jlo+1,khi+1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((ntop .gt. 0) .and. (ndwn .gt. 0)) then
         if (klo+2 .le. ke) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,klo-1) = half *
     $              (3*q(i,jhi+1,klo) - 3*q(i,jhi+1,klo+1) + 
     $              q(i,jhi+1,klo+2))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,klo-1) = half *
     $              (2*q(i,jhi+1,klo) - q(i,jhi+1,klo+1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + half *
     $              (3*q(i,jhi,klo-1) - 3*q(i,jhi-1,klo-1) + 
     $              q(i,jhi-2,klo-1))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + half *
     $              (2*q(i,jhi,klo-1) - q(i,jhi-1,klo-1))
            end do
         end if
      end if
c
c ****************************************************************************
c
      if ((ntop .gt. 0) .and. (nup .gt. 0)) then
         if (khi-2 .ge. ks) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,khi+1) = half *
     $              (3*q(i,jhi+1,khi) - 3*q(i,jhi+1,khi-1) + 
     $              q(i,jhi+1,khi-2))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,khi+1) = half *
     $              (2*q(i,jhi+1,khi) - q(i,jhi+1,khi-1))
            end do
         end if
         if (jhi-2 .ge. js) then
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + half *
     $              (3*q(i,jhi,khi+1) - 3*q(i,jhi-1,khi+1) + 
     $              q(i,jhi-2,khi+1))
            end do
         else
            do i = ARG_L1(q),ARG_H1(q)
               q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + half *
     $              (2*q(i,jhi,khi+1) - q(i,jhi-1,khi+1))
            end do
         end if
      end if

      end
