
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_BC_TYPES.H"

#define SDIM 2

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

subroutine filcc(q,q_l1,q_l2,q_h1,q_h2,domlo,domhi,dx,xlo,bc)

  implicit none

  integer    q_l1, q_l2, q_h1, q_h2
  integer    domlo(SDIM), domhi(SDIM)
  integer    bc(SDIM,2)
  REAL_T     xlo(SDIM), dx(SDIM)
  REAL_T     q(q_l1:q_h1,q_l2:q_h2)

  integer    nlft, nrgt, nbot, ntop
  integer    ilo, ihi, jlo, jhi
  integer    i, j
  integer    is, ie, js, je

  nlft = max(0,domlo(1)-q_l1)
  nrgt = max(0,q_h1-domhi(1))
  nbot = max(0,domlo(2)-q_l2)
  ntop = max(0,q_h2-domhi(2))

  is = max(q_l1,domlo(1))
  ie = min(q_h1,domhi(1))
  js = max(q_l2,domlo(2))
  je = min(q_h2,domhi(2))

  !     ::::: first fill sides
  if (nlft .gt. 0) then
     ilo = domlo(1)

     if (bc(1,1) .eq. EXT_DIR) then
        !     set all ghost cell values to a prescribed dirichlet
        !     value; in this example, we have chosen 1
        !	    do i = 1, nlft
        !	    do j = q_l2, q_h2
        !	       q(ilo-i,j) = 1.d0
        !	    end do
        !	    end do
     else if (bc(1,1) .eq. FOEXTRAP) then
        do i = 1, nlft
           do j = q_l2, q_h2
              q(ilo-i,j) = q(ilo,j)
           end do
        end do
     else if (bc(1,1) .eq. HOEXTRAP) then
        do i = 2, nlft
           do j = q_l2, q_h2
              q(ilo-i,j) = q(ilo,j) 
           end do
        end do
        if (ilo+2 .le. ie) then 
           do j = q_l2, q_h2
              q(ilo-1,j) = (15*q(ilo,j) - 10*q(ilo+1,j) + &
                   3*q(ilo+2,j)) * eighth
           end do
        else 
           do j = q_l2, q_h2
              q(ilo-1,j) = half*(3*q(ilo,j) - q(ilo+1,j))
           end do
        end if
     else if (bc(1,1) .eq. REFLECT_EVEN) then
        do i = 1, nlft
           do j = q_l2, q_h2
              q(ilo-i,j) = q(ilo+i-1,j)
           end do
        end do
     else if (bc(1,1) .eq. REFLECT_ODD) then
        do i = 1, nlft
           do j = q_l2, q_h2
              q(ilo-i,j) = -q(ilo+i-1,j)
           end do
        end do
     end if
  end if

  if (nrgt .gt. 0) then
     ihi = domhi(1)

     if (bc(1,2) .eq. EXT_DIR) then
        !	    do i = 1, nrgt
        !	    do j = q_l2, q_h2
        !	       q(ihi+i,j) = 1.d0
        !	    end do
        !	    end do
     else if (bc(1,2) .eq. FOEXTRAP) then
        do i = 1, nrgt
           do j = q_l2, q_h2
              q(ihi+i,j) = q(ihi,j)
           end do
        end do
     else if (bc(1,2) .eq. HOEXTRAP) then
        do i = 2, nrgt
           do j = q_l2, q_h2
              q(ihi+i,j) = q(ihi,j)
           end do
        end do
        if (ihi-2 .ge. is) then
           do j = q_l2, q_h2
              q(ihi+1,j) = (15*q(ihi,j) - 10*q(ihi-1,j) + &
                   3*q(ihi-2,j)) * eighth
           end do
        else
           do j = q_l2, q_h2
              q(ihi+1,j) = half*(3*q(ihi,j) - q(ihi-1,j))
           end do
        end if
     else if (bc(1,2) .eq. REFLECT_EVEN) then
        do i = 1, nrgt
           do j = q_l2, q_h2
              q(ihi+i,j) = q(ihi-i+1,j)
           end do
        end do
     else if (bc(1,2) .eq. REFLECT_ODD) then
        do i = 1, nrgt
           do j = q_l2, q_h2
              q(ihi+i,j) = -q(ihi-i+1,j)
           end do
        end do
     end if
  end if

  if (nbot .gt. 0) then
     jlo = domlo(2)

     if (bc(2,1) .eq. EXT_DIR) then
        !	    do j = 1, nbot
        !	    do i = q_l1, q_h1
        !	       q(i,jlo-j) = 1.d0
        !	    end do
        !	    end do
     else if (bc(2,1) .eq. FOEXTRAP) then
        do j = 1, nbot
           do i = q_l1, q_h1
              q(i,jlo-j) = q(i,jlo)
           end do
        end do
     else if (bc(2,1) .eq. HOEXTRAP) then
        do j = 2, nbot
           do i = q_l1, q_h1
              q(i,jlo-j) = q(i,jlo)
           end do
        end do
        if (jlo+2 .le. je) then
           do i = q_l1, q_h1
              q(i,jlo-1) = (15*q(i,jlo) - 10*q(i,jlo+1) + &
                   3*q(i,jlo+2)) * eighth
           end do
        else
           do i = q_l1, q_h1
              q(i,jlo-1) = half*(3*q(i,jlo) - q(i,jlo+1))
           end do
        end if
     else if (bc(2,1) .eq. REFLECT_EVEN) then
        do j = 1, nbot
           do i = q_l1, q_h1
              q(i,jlo-j) = q(i,jlo+j-1)
           end do
        end do
     else if (bc(2,1) .eq. REFLECT_ODD) then
        do j = 1, nbot
           do i = q_l1, q_h1
              q(i,jlo-j) = -q(i,jlo+j-1)
           end do
        end do
     end if
  end if

  if (ntop .gt. 0) then
     jhi = domhi(2)

     if (bc(2,2) .eq. EXT_DIR) then
        !	    do j = 1, ntop
        ! 	    do i = q_l1, q_h1
        !	       q(i,jhi+j) = 1.d0
        !	    end do
        !	    end do
     else if (bc(2,2) .eq. FOEXTRAP) then
        do j = 1, ntop
           do i = q_l1, q_h1
              q(i,jhi+j) = q(i,jhi)
           end do
        end do
     else if (bc(2,2) .eq. HOEXTRAP) then
        do j = 2, ntop
           do i = q_l1, q_h1
              q(i,jhi+j) = q(i,jhi)
           end do
        end do
        if (jhi-2 .ge. js) then
           do i = q_l1, q_h1
              q(i,jhi+1) = (15*q(i,jhi) - 10*q(i,jhi-1) + &
                   3*q(i,jhi-2)) * eighth
           end do
        else
           do i = q_l1, q_h1
              q(i,jhi+1) = half*(3*q(i,jhi) - q(i,jhi-1))
           end do
        end if
     else if (bc(2,2) .eq. REFLECT_EVEN) then
        do j = 1, ntop
           do i = q_l1, q_h1
              q(i,jhi+j) = q(i,jhi-j+1)
           end do
        end do
     else if (bc(2,2) .eq. REFLECT_ODD) then
        do j = 1, ntop
           do i = q_l1, q_h1
              q(i,jhi+j) = -q(i,jhi-j+1)
           end do
        end do
     end if
  end if

  if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and. &
       (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) ) then

     if (jlo+2 .le. je) then 
        q(ilo-1,jlo-1) = half * eighth * &
             (15*q(ilo-1,jlo) - 10*q(ilo-1,jlo+1) + 3*q(ilo-1,jlo+2))
     else
        q(ilo-1,jlo-1) = half * half * &
             (3*q(ilo-1,jlo) - q(ilo-1,jlo+1))
     end if

     if (ilo+2 .le. ie) then 
        q(ilo-1,jlo-1) =  q(ilo-1,jlo-1) + half * eighth * &
             (15*q(ilo,jlo-1) - 10*q(ilo+1,jlo-1) + 3*q(ilo+2,jlo-1)) 
     else
        q(ilo-1,jlo-1) =  q(ilo-1,jlo-1) + half * half * &
             (3*q(ilo,jlo-1) - q(ilo+1,jlo-1))
     end if


  end if

  if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and. &
       (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) ) then

     if (jhi-2 .ge. js) then 
        q(ilo-1,jhi+1) = half * eighth * &
             (15*q(ilo-1,jhi) - 10*q(ilo-1,jhi-1) + 3*q(ilo-1,jhi-2))
     else
        q(ilo-1,jhi+1) = half * half * &
             (3*q(ilo-1,jhi) - q(ilo-1,jhi-1))
     end if

     if (ilo+2 .le. ie) then 
        q(ilo-1,jhi+1) = q(ilo-1,jhi+1) + half * eighth * &
             (15*q(ilo,jhi+1) - 10*q(ilo+1,jhi+1) + 3*q(ilo+2,jhi+1))
     else
        q(ilo-1,jhi+1) = q(ilo-1,jhi+1) + half * half * &
             (3*q(ilo,jhi+1) - q(ilo+1,jhi+1))
     end if
  end if

  if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and. &
       (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) ) then
     if (jlo+2 .le. je) then 
        q(ihi+1,jlo-1) = half * eighth * &
             (15*q(ihi+1,jlo) - 10*q(ihi+1,jlo+1) + 3*q(ihi+1,jlo+2))
     else
        q(ihi+1,jlo-1) = half * half * &
             (3*q(ihi+1,jlo) - q(ihi+1,jlo+1))
     end if

     if (ihi-2 .ge. is) then 
        q(ihi+1,jlo-1) = q(ihi+1,jlo-1) + half * eighth * &
             (15*q(ihi,jlo-1) - 10*q(ihi-1,jlo-1) + 3*q(ihi-2,jlo-1))
     else
        q(ihi+1,jlo-1) = q(ihi+1,jlo-1) + half * half * &
             (3*q(ihi,jlo-1) - q(ihi-1,jlo-1))
     end if
  end if

  if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and. &
       (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) ) then

     if (jhi-2 .ge. js) then 
        q(ihi+1,jhi+1) = half * eighth * &
             (15*q(ihi+1,jhi) - 10*q(ihi+1,jhi-1) + 3*q(ihi+1,jhi-2))
     else
        q(ihi+1,jhi+1) = half * half * &
             (3*q(ihi+1,jhi) - q(ihi+1,jhi-1))
     end if

     if (ihi-2 .ge. is) then 
        q(ihi+1,jhi+1) = q(ihi+1,jhi+1) + half * eighth * &
             (15*q(ihi,jhi+1) - 10*q(ihi-1,jhi+1) + 3*q(ihi-2,jhi+1))
     else
        q(ihi+1,jhi+1) = q(ihi+1,jhi+1) + half * half * &
             (3*q(ihi,jhi+1) - q(ihi-1,jhi+1))
     end if

  end if

  return
end subroutine filcc

subroutine hoextraptocc(q,q_l1,q_l2,q_h1,q_h2,domlo,domhi,dx,xlo)

  implicit none

  integer    q_l1, q_l2, q_h1, q_h2
  integer    domlo(SDIM), domhi(SDIM)
  REAL_T     xlo(SDIM), dx(SDIM)
  REAL_T     q(q_l1:q_h1,q_l2:q_h2)

  integer    nlft, nrgt, nbot, ntop
  integer    ilo, ihi, jlo, jhi
  integer    i, j
  integer    is, ie, js, je

  nlft = max(0,domlo(1)-q_l1)
  nrgt = max(0,q_h1-domhi(1))
  nbot = max(0,domlo(2)-q_l2)
  ntop = max(0,q_h2-domhi(2))

  is = max(q_l1,domlo(1))
  ie = min(q_h1,domhi(1))
  js = max(q_l2,domlo(2))
  je = min(q_h2,domhi(2))

  !
  !     Set these to invalid values, they shouldn't be used if not reset
  !
  ilo = -10
  jlo = -10
  ihi = 100000000
  jhi = 100000000

  !
  !     First fill sides.
  !
  if (nlft .gt. 0) then
     ilo = domlo(1)
     do i = 2, nlft
        do j = q_l2, q_h2
           q(ilo-i,j) = q(ilo,j) 
        end do
     end do
     if (ilo+2 .le. ie) then 
        do j = q_l2, q_h2
           q(ilo-1,j) = 3*q(ilo,j) - 3*q(ilo+1,j) + q(ilo+2,j)
        end do
     else 
        do j = q_l2, q_h2
           q(ilo-1,j) = 2*q(ilo,j) - q(ilo+1,j)
        end do
     end if
  end if

  if (nrgt .gt. 0) then
     ihi = domhi(1)
     do i = 2, nrgt
        do j = q_l2, q_h2
           q(ihi+i,j) = q(ihi,j)
        end do
     end do
     if (ihi-2 .ge. is) then
        do j = q_l2, q_h2
           q(ihi+1,j) = 3*q(ihi,j) - 3*q(ihi-1,j) + q(ihi-2,j)
        end do
     else
        do j = q_l2, q_h2
           q(ihi+1,j) = 2*q(ihi,j) - q(ihi-1,j)
        end do
     end if
  end if

  if (nbot .gt. 0) then
     jlo = domlo(2)
     do j = 2, nbot
        do i = q_l1, q_h1
           q(i,jlo-j) = q(i,jlo)
        end do
     end do
     if (jlo+2 .le. je) then
        do i = q_l1, q_h1
           q(i,jlo-1) = 3*q(i,jlo) - 3*q(i,jlo+1) + q(i,jlo+2)
        end do
     else
        do i = q_l1, q_h1
           q(i,jlo-1) = 2*q(i,jlo) - q(i,jlo+1)
        end do
     end if
  end if

  if (ntop .gt. 0) then
     jhi = domhi(2)
     do j = 2, ntop
        do i = q_l1, q_h1
           q(i,jhi+j) = q(i,jhi)
        end do
     end do
     if (jhi-2 .ge. js) then
        do i = q_l1, q_h1
           q(i,jhi+1) = 3*q(i,jhi) - 3*q(i,jhi-1) + q(i,jhi-2)
        end do
     else
        do i = q_l1, q_h1
           q(i,jhi+1) = 2*q(i,jhi) - q(i,jhi-1)
        end do
     end if
  end if

  if ((nlft .gt. 0) .and. (nbot .gt. 0)) then
     if (jlo+2 .le. je) then
        q(ilo-1,jlo-1) = half * &
             (3*q(ilo-1,jlo) - 3*q(ilo-1,jlo+1) + q(ilo-1,jlo+2))
     else
        q(ilo-1,jlo-1) = half * (2*q(ilo-1,jlo) - q(ilo-1,jlo+1))
     end if

     if (ilo+2 .le. ie) then 
        q(ilo-1,jlo-1) =  q(ilo-1,jlo-1) + half * &
             (3*q(ilo,jlo-1) - 3*q(ilo+1,jlo-1) + q(ilo+2,jlo-1)) 
     else
        q(ilo-1,jlo-1) =  q(ilo-1,jlo-1) + half * &
             (2*q(ilo,jlo-1) - q(ilo+1,jlo-1))
     end if
  end if

  if ((nlft .gt. 0) .and. (ntop .gt. 0)) then 
     if (jhi-2 .ge. js) then 
        q(ilo-1,jhi+1) = half * &
             (3*q(ilo-1,jhi) - 3*q(ilo-1,jhi-1) + q(ilo-1,jhi-2))
     else
        q(ilo-1,jhi+1) = half * (2*q(ilo-1,jhi) - q(ilo-1,jhi-1))
     end if

     if (ilo+2 .le. ie) then 
        q(ilo-1,jhi+1) = q(ilo-1,jhi+1) + half * &
             (3*q(ilo,jhi+1) - 3*q(ilo+1,jhi+1) + q(ilo+2,jhi+1))
     else
        q(ilo-1,jhi+1) = q(ilo-1,jhi+1) + half * &
             (2*q(ilo,jhi+1) - q(ilo+1,jhi+1))
     end if
  end if

  if ((nrgt .gt. 0) .and. (nbot .gt. 0)) then 
     if (jlo+2 .le. je) then 
        q(ihi+1,jlo-1) = half * &
             (3*q(ihi+1,jlo) - 3*q(ihi+1,jlo+1) + q(ihi+1,jlo+2))
     else
        q(ihi+1,jlo-1) = half * (2*q(ihi+1,jlo) - q(ihi+1,jlo+1))
     end if

     if (ihi-2 .ge. is) then 
        q(ihi+1,jlo-1) = q(ihi+1,jlo-1) + half * &
             (3*q(ihi,jlo-1) - 3*q(ihi-1,jlo-1) + q(ihi-2,jlo-1))
     else
        q(ihi+1,jlo-1) = q(ihi+1,jlo-1) + half * &
             (2*q(ihi,jlo-1) - q(ihi-1,jlo-1))
     end if
  end if

  if ((nrgt .gt. 0) .and. (ntop .gt. 0)) then 
     if (jhi-2 .ge. js) then 
        q(ihi+1,jhi+1) = half * &
             (3*q(ihi+1,jhi) - 3*q(ihi+1,jhi-1) + q(ihi+1,jhi-2))
     else
        q(ihi+1,jhi+1) = half * (2*q(ihi+1,jhi) - q(ihi+1,jhi-1))
     end if

     if (ihi-2 .ge. is) then 
        q(ihi+1,jhi+1) = q(ihi+1,jhi+1) + half * &
             (3*q(ihi,jhi+1) - 3*q(ihi-1,jhi+1) + q(ihi-2,jhi+1))
     else
        q(ihi+1,jhi+1) = q(ihi+1,jhi+1) + half * &
             (2*q(ihi,jhi+1) - q(ihi-1,jhi+1))
     end if
  end if

end subroutine hoextraptocc





module filcc_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

#ifdef AMREX_USE_CUDA
  attributes(device) &
#endif
  subroutine filccn(blo, bhi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)

    implicit none

    integer,  intent(in   ) :: blo(3), bhi(3)
    integer,  intent(in   ) :: q_lo(3), q_hi(3)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    integer,  intent(in   ) :: bc(2,2,ncomp)
    real(rt), intent(in   ) :: xlo(2), dx(2)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),ncomp)
    integer,  intent(in   ) :: ncomp

    integer :: ilo, ihi, jlo, jhi
    integer :: is, ie, js, je
    integer :: i, j, n

    is = max(q_lo(1), domlo(1))
    ie = min(q_hi(1), domhi(1))
    js = max(q_lo(2), domlo(2))
    je = min(q_hi(2), domhi(2))

    ilo = domlo(1)
    ihi = domhi(1)
    jlo = domlo(2)
    jhi = domhi(2)

    do n = 1, ncomp
       do j = blo(2), bhi(2)
          do i = blo(1), bhi(1)

             if (bc(1,1,n) .eq. EXT_DIR) then

                ! Do nothing.

             else if (bc(1,1,n) .eq. FOEXTRAP) then

                if (i < ilo) then
                   q(i,j,n) = q(ilo,j,n)
                end if

             else if (bc(1,1,n) .eq. HOEXTRAP) then

                if (i < ilo - 1) then
                   q(i,j,n) = q(ilo,j,n)
                else if (i == ilo - 1) then
                   if (ilo+2 <= ie) then
                      q(i,j,n) = eighth * (15*q(ilo,j,n) - 10*q(ilo+1,j,n) + 3*q(ilo+2,j,n))
                   else
                      q(i,j,n) = half * (3*q(ilo,j,n) - q(ilo+1,j,n))
                   end if
                end if

             else if (bc(1,1,n) .eq. REFLECT_EVEN) then

                if (i < ilo) then
                   q(i,j,n) = q(ilo+(ilo-i)-1,j,n)
                end if

             else if (bc(1,1,n) .eq. REFLECT_ODD) then

                if (i < ilo) then
                   q(i,j,n) = -q(ilo+(ilo-i)-1,j,n)
                end if

             end if



             if (bc(1,2,n) .eq. EXT_DIR) then

                ! Do nothing.

             else if (bc(1,2,n) .eq. FOEXTRAP) then

                if (i > ihi) then
                   q(i,j,n) = q(ihi,j,n)
                end if

             else if (bc(1,2,n) .eq. HOEXTRAP) then

                if (i > ihi + 1) then
                   q(i,j,n) = q(ihi,j,n)
                else if (i == ihi + 1) then
                   if (ihi-2 >= is) then
                      q(i,j,n) = eighth * (15*q(ihi,j,n) - 10*q(ihi-1,j,n) + 3*q(ihi-2,j,n))
                   else
                      q(i,j,n) = half * (3*q(ihi,j,n) - q(ihi-1,j,n))
                   end if
                end if

             else if (bc(1,2,n) .eq. REFLECT_EVEN) then

                if (i > ihi) then
                   q(i,j,n) = q(ihi-(i-ihi)+1,j,n)
                end if

             else if (bc(1,2,n) .eq. REFLECT_ODD) then

                if (i > ihi) then
                   q(i,j,n) = -q(ihi-(i-ihi)+1,j,n)
                end if

             end if



             if (bc(2,1,n) .eq. EXT_DIR) then

                ! Do nothing.

             else if (bc(2,1,n) .eq. FOEXTRAP) then

                if (j < jlo) then
                   q(i,j,n) = q(i,jlo,n)
                end if

             else if (bc(2,1,n) .eq. HOEXTRAP) then

                if (j < jlo - 1) then
                   q(i,j,n) = q(i,jlo,n)
                else if (j == jlo - 1) then
                   if (jlo+2 <= je) then
                      q(i,j,n) = eighth * (15*q(i,jlo,n) - 10*q(i,jlo+1,n) + 3*q(i,jlo+2,n))
                   else
                      q(i,j,n) = half * (3*q(i,jlo,n) - q(i,jlo+1,n))
                   end if
                end if

             else if (bc(2,1,n) .eq. REFLECT_EVEN) then

                if (j < jlo) then
                   q(i,j,n) = q(i,jlo+(jlo-j)-1,n)
                end if

             else if (bc(2,1,n) .eq. REFLECT_ODD) then

                if (j < jlo) then
                   q(i,j,n) = -q(i,jlo+(jlo-j)-1,n)
                end if

             end if



             if (bc(2,2,n) .eq. EXT_DIR) then

                ! Do nothing.

             else if (bc(2,2,n) .eq. FOEXTRAP) then

                if (j > jhi) then
                   q(i,j,n) = q(i,jhi,n)
                end if

             else if (bc(2,2,n) .eq. HOEXTRAP) then

                if (j > jhi + 1) then
                   q(i,j,n) = q(i,jhi,n)
                else if (j == jhi + 1) then
                   if (jhi-2 >= js) then
                      q(i,j,n) = eighth * (15*q(i,jhi,n) - 10*q(i,jhi-1,n) + 3*q(i,jhi-2,n))
                   else
                      q(i,j,n) = half * (3*q(i,jhi,n) - q(i,jhi-1,n))
                   end if
                end if

             else if (bc(2,2,n) .eq. REFLECT_EVEN) then

                if (j > jhi) then
                   q(i,j,n) = q(i,jhi-(j-jhi)+1,n)
                end if

             else if (bc(2,2,n) .eq. REFLECT_ODD) then

                if (j > jhi) then
                   q(i,j,n) = -q(i,jhi-(j-jhi)+1,n)
                end if

             end if


             if (i == ilo-1 .and. bc(1,1,n) .eq. HOEXTRAP .and. &
                 j == jlo-1 .and. bc(2,1,n) .eq. HOEXTRAP) then

                if (jlo+2 <= je) then
                   q(i,j,n) = half * eighth * (15*q(ilo-1,jlo,n) - 10*q(ilo-1,jlo+1,n) + 3*q(ilo-1,jlo+2,n))
                else
                   q(i,j,n) = half * half * (3*q(ilo-1,jlo,n) - q(ilo-1,jlo+1,n))
                end if

                if (ilo+2 <= ie) then
                   q(i,j,n) = q(ilo-1,jlo-1,n) + &
                              half * eighth * (15*q(ilo,jlo-1,n) - 10*q(ilo+1,jlo-1,n) + 3*q(ilo+2,jlo-1,n))
                else
                   q(i,j,n) = q(ilo-1,jlo-1,n) + half * half * (3*q(ilo,jlo-1,n) - q(ilo+1,jlo-1,n))
                end if

             end if

             if (i == ilo-1 .and. bc(1,1,n) .eq. HOEXTRAP .and. &
                 j == jhi+1 .and. bc(2,2,n) .eq. HOEXTRAP) then

                if (jhi-2 >= js) then
                   q(i,j,n) = half * eighth * (15*q(ilo-1,jhi,n) - 10*q(ilo-1,jhi-1,n) + 3*q(ilo-1,jhi-2,n))
                else
                   q(i,j,n) = half * half * (3*q(ilo-1,jhi,n) - q(ilo-1,jhi-1,n))
                end if

                if (ilo+2 <= ie) then
                   q(i,j,n) = q(ilo-1,jhi+1,n) + &
                              half * eighth * (15*q(ilo,jhi+1,n) - 10*q(ilo+1,jhi+1,n) + 3*q(ilo+2,jhi+1,n))
                else
                   q(i,j,n) = q(ilo-1,jhi+1,n) + half * half * (3*q(ilo,jhi+1,n) - q(ilo+1,jhi+1,n))
                end if

             end if

             if (i == ihi+1 .and. bc(1,2,n) .eq. HOEXTRAP .and. &
                 j == jlo-1 .and. bc(2,1,n) .eq. HOEXTRAP) then

                if (jlo+2 <= je) then
                   q(i,j,n) = half * eighth * (15*q(ihi+1,jlo,n) - 10*q(ihi+1,jlo+1,n) + 3*q(ihi+1,jlo+2,n))
                else
                   q(i,j,n) = half * half * (3*q(ihi+1,jlo,n) - q(ihi+1,jlo+1,n))
                end if

                if (ihi-2 >= is) then
                   q(i,j,n) = q(ihi+1,jlo-1,n) + &
                              half * eighth * (15*q(ihi,jlo-1,n) - 10*q(ihi-1,jlo-1,n) + 3*q(ihi-2,jlo-1,n))
                else
                   q(i,j,n) = q(ihi+1,jlo-1,n) + half * half * (3*q(ihi,jlo-1,n) - q(ihi-1,jlo-1,n))
                end if

             end if

             if (i == ihi+1 .and. bc(1,2,n) .eq. HOEXTRAP .and. &
                 j == jhi+1 .and. bc(2,2,n) .eq. HOEXTRAP) then

                if (jhi-2 >= js) then
                   q(i,j,n) = half * eighth * (15*q(ihi+1,jhi,n) - 10*q(ihi+1,jhi-1,n) + 3*q(ihi+1,jhi-2,n))
                else
                   q(i,j,n) = half * half * (3*q(ihi+1,jhi,n) - q(ihi+1,jhi-1,n))
                end if

                if (ihi-2 >= is) then
                   q(i,j,n) = q(ihi+1,jhi+1,n) + &
                              half * eighth * (15*q(ihi,jhi+1,n) - 10*q(ihi-1,jhi+1,n) + 3*q(ihi-2,jhi+1,n))
                else
                   q(i,j,n) = q(ihi+1,jhi+1,n) + half * half * (3*q(ihi,jhi+1,n) - q(ihi-1,jhi+1,n))
                end if

             end if

          end do
       end do
    end do
  end subroutine filccn

end module filcc_module
