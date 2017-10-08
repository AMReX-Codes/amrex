
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_CONSTANTS.H"
#include "AMReX_BC_TYPES.H"

! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for cell centered data.  It knows how to exrapolate,
! ::: and reflect data and can be used to suppliment problem
! ::: specific fill functions (ie. EXT_DIR).
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: q        <=  array to fill
! ::: DIMS(q)   => index extent of q array
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	           corner of q array
! ::: bc	=> array of boundary flags bc(SPACEDIM,lo:hi)
! ::: 
! ::: NOTE: corner data not used in computing soln but must have
! :::       reasonable values for arithmetic to live
! ::: -----------------------------------------------------------

subroutine filcc(q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,domlo,domhi,dx,xlo,bc)

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in   ) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer,  intent(in   ) :: domlo(3), domhi(3)
  real(rt), intent(in   ) :: xlo(3), dx(3)
  real(rt), intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
  integer,  intent(in   ) :: bc(3,2)

  integer :: nlft, nrgt, nbot, ntop, nup, ndwn
  integer :: ilo, ihi, jlo, jhi, klo, khi
  integer :: is,  ie,  js,  je,  ks,  ke
  integer :: i, j, k

  is = max(q_l1,domlo(1))
  ie = min(q_h1,domhi(1))
  js = max(q_l2,domlo(2))
  je = min(q_h2,domhi(2))
  ks = max(q_l3,domlo(3))
  ke = min(q_h3,domhi(3))

  nlft = max(0,domlo(1)-q_l1)
  nrgt = max(0,q_h1-domhi(1))
  nbot = max(0,domlo(2)-q_l2)
  ntop = max(0,q_h2-domhi(2))
  ndwn = max(0,domlo(3)-q_l3)
  nup  = max(0,q_h3-domhi(3))

  !
  !     ::::: first fill sides
  !
  if (nlft .gt. 0) then
     ilo = domlo(1)

     if (bc(1,1) .eq. EXT_DIR) then
        !     set all ghost cell values to a prescribed dirichlet
        !     value; in this example, we have chosen 1
        !	    do i = 1, nlft
        !           do k = q_l3,q_h3
        !           do j = q_l2,q_h2
        !              q(ilo-i,j,k) = 1.d0
        !           end do
        !           end do
        !	    end do
     else if (bc(1,1) .eq. FOEXTRAP) then
        do i = 1, nlft
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ilo-i,j,k) = q(ilo,j,k)
              end do
           end do
        end do
     else if (bc(1,1) .eq. HOEXTRAP) then
        do i = 2, nlft
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ilo-i,j,k) = q(ilo,j,k)
              end do
           end do
        end do
        if (ilo+2 .le. ie) then
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ilo-1,j,k) = (15*q(ilo,j,k) - 10*q(ilo+1,j,k) + &
                      3*q(ilo+2,j,k)) * eighth
              end do
           end do
        else  
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ilo-1,j,k) = half*(3*q(ilo,j,k) - q(ilo+1,j,k))
              end do
           end do
        end if
     else if (bc(1,1) .eq. REFLECT_EVEN) then
        do i = 1, nlft
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ilo-i,j,k) = q(ilo+i-1,j,k)
              end do
           end do
        end do
     else if (bc(1,1) .eq. REFLECT_ODD) then
        do i = 1, nlft
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ilo-i,j,k) = -q(ilo+i-1,j,k)
              end do
           end do
        end do
     end if
  end if

  if (nrgt .gt. 0) then
     ihi = domhi(1)

     if (bc(1,2) .eq. EXT_DIR) then
        !	    do i = 1, nrgt
        !           do k = q_l3,q_h3
        !           do j = q_l2,q_h2
        !              q(ihi+i,j,k) = 1.d0
        !           end do
        !           end do
        !	    end do
     else if (bc(1,2) .eq. FOEXTRAP) then
        do i = 1, nrgt
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ihi+i,j,k) = q(ihi,j,k)
              end do
           end do
        end do
     else if (bc(1,2) .eq. HOEXTRAP) then
        do i = 2, nrgt
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ihi+i,j,k) = q(ihi,j,k)
              end do
           end do
        end do
        if (ihi-2 .ge. is) then
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ihi+1,j,k) = (15*q(ihi,j,k) - 10*q(ihi-1,j,k) + &
                      3*q(ihi-2,j,k)) * eighth
              end do
           end do
        else
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ihi+1,j,k) = half*(3*q(ihi,j,k) - q(ihi-1,j,k))
              end do
           end do
        end if
     else if (bc(1,2) .eq. REFLECT_EVEN) then
        do i = 1, nrgt
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ihi+i,j,k) = q(ihi-i+1,j,k)
              end do
           end do
        end do
     else if (bc(1,2) .eq. REFLECT_ODD) then
        do i = 1, nrgt
           do k = q_l3,q_h3
              do j = q_l2,q_h2
                 q(ihi+i,j,k) = -q(ihi-i+1,j,k)
              end do
           end do
        end do
     end if
  end if

  if (nbot .gt. 0) then
     jlo = domlo(2)

     if (bc(2,1) .eq. EXT_DIR) then
        !	    do j = 1, nbot
        !           do k = q_l3,q_h3
        !           do i = q_l1,q_h1
        !              q(i,jlo-j,k) = 1.d0
        !           end do
        !           end do
        !	    end do
     else if (bc(2,1) .eq. FOEXTRAP) then
        do j = 1, nbot
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jlo-j,k) = q(i,jlo,k)
              end do
           end do
        end do
     else if (bc(2,1) .eq. HOEXTRAP) then
        do j = 2, nbot
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jlo-j,k) = q(i,jlo,k)
              end do
           end do
        end do
        if (jlo+2 .le. je) then
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jlo-1,k) = (15*q(i,jlo,k) - 10*q(i,jlo+1,k) + &
                      3*q(i,jlo+2,k)) * eighth
              end do
           end do
        else
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jlo-1,k) = half*(3*q(i,jlo,k) - q(i,jlo+1,k))
              end do
           end do
        end if
     else if (bc(2,1) .eq. REFLECT_EVEN) then
        do j = 1, nbot 
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jlo-j,k) = q(i,jlo+j-1,k)
              end do
           end do
        end do
     else if (bc(2,1) .eq. REFLECT_ODD) then
        do j = 1, nbot
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jlo-j,k) = -q(i,jlo+j-1,k)
              end do
           end do
        end do
     end if
  end if

  if (ntop .gt. 0) then
     jhi = domhi(2)

     if (bc(2,2) .eq. EXT_DIR) then
        !	    do j = 1, ntop
        !           do k = q_l3,q_h3
        !           do i = q_l1,q_h1
        !              q(i,jhi+j,k) = 1.d0
        !           end do
        !           end do
        !	    end do
     else if (bc(2,2) .eq. FOEXTRAP) then
        do j = 1, ntop
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jhi+j,k) = q(i,jhi,k)
              end do
           end do
        end do
     else if (bc(2,2) .eq. HOEXTRAP) then
        do j = 2, ntop
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jhi+j,k) = q(i,jhi,k)
              end do
           end do
        end do
        if (jhi-2 .ge. js) then
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jhi+1,k) = (15*q(i,jhi,k) - 10*q(i,jhi-1,k) + &
                      3*q(i,jhi-2,k)) * eighth
              end do
           end do
        else
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jhi+1,k) = half*(3*q(i,jhi,k) - q(i,jhi-1,k))
              end do
           end do
        end if
     else if (bc(2,2) .eq. REFLECT_EVEN) then
        do j = 1, ntop
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jhi+j,k) = q(i,jhi-j+1,k)
              end do
           end do
        end do
     else if (bc(2,2) .eq. REFLECT_ODD) then
        do j = 1, ntop
           do k = q_l3,q_h3
              do i = q_l1,q_h1
                 q(i,jhi+j,k) = -q(i,jhi-j+1,k)
              end do
           end do
        end do
     end if
  end if

  if (ndwn .gt. 0) then
     klo = domlo(3)

     if (bc(3,1) .eq. EXT_DIR) then
        !	    do k = 1, ndwn
        !           do j = q_l2,q_h2
        !           do i = q_l1,q_h1
        !              q(i,j,klo-k) = 1.d0
        !           end do
        !           end do
        !	    end do
     else if (bc(3,1) .eq. FOEXTRAP) then
        do k = 1, ndwn
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,klo-k) = q(i,j,klo)
              end do
           end do
        end do
     else if (bc(3,1) .eq. HOEXTRAP) then
        do k = 2, ndwn
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,klo-k) = q(i,j,klo)
              end do
           end do
        end do
        if (klo+2 .le. ke) then
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,klo-1) = (15*q(i,j,klo) - 10*q(i,j,klo+1) + &
                      3*q(i,j,klo+2)) * eighth
              end do
           end do
        else
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,klo-1) = half*(3*q(i,j,klo) - q(i,j,klo+1))
              end do
           end do
        end if
     else if (bc(3,1) .eq. REFLECT_EVEN) then
        do k = 1, ndwn
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,klo-k) = q(i,j,klo+k-1)
              end do
           end do
        end do
     else if (bc(3,1) .eq. REFLECT_ODD) then
        do k = 1, ndwn
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,klo-k) = -q(i,j,klo+k-1)
              end do
           end do
        end do
     end if
  end if

  if (nup .gt. 0) then
     khi = domhi(3)

     if (bc(3,2) .eq. EXT_DIR) then
        !	    do k = 1, nup
        !           do j = q_l2,q_h2
        !           do i = q_l1,q_h1
        !              q(i,j,khi+k) = 1.d0
        !           end do
        !           end do
        !	    end do
     else if (bc(3,2) .eq. FOEXTRAP) then
        do k = 1, nup
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,khi+k) = q(i,j,khi)
              end do
           end do
        end do
     else if (bc(3,2) .eq. HOEXTRAP) then
        do k = 2, nup
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,khi+k) = q(i,j,khi)
              end do
           end do
        end do
        if (khi-2 .ge. ks) then
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,khi+1) = (15*q(i,j,khi) - 10*q(i,j,khi-1) + &
                      3*q(i,j,khi-2)) * eighth
              end do
           end do
        else
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,khi+1) = half*(3*q(i,j,khi) - q(i,j,khi-1))
              end do
           end do
        end if
     else if (bc(3,2) .eq. REFLECT_EVEN) then
        do k = 1, nup
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,khi+k) = q(i,j,khi-k+1)
              end do
           end do
        end do
     else if (bc(3,2) .eq. REFLECT_ODD) then
        do k = 1, nup
           do j = q_l2,q_h2
              do i = q_l1,q_h1
                 q(i,j,khi+k) = -q(i,j,khi-k+1)
              end do
           end do
        end do
     end if
  end if
  !
  !    First correct the i-j edges and all corners
  !
  if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and. &
       (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) ) then
     if (jlo+2 .le. je) then
        do k = q_l3,q_h3
           q(ilo-1,jlo-1,k) = half * eighth * &
                (15*q(ilo-1,jlo,k) - 10*q(ilo-1,jlo+1,k) + &
                3*q(ilo-1,jlo+2,k))
        end do
     else
        do k = q_l3,q_h3
           q(ilo-1,jlo-1,k) = half * half * &
                (3*q(ilo-1,jlo,k) - q(ilo-1,jlo+1,k))
        end do
     end if

     if (ilo+2 .le. ie) then
        do k = q_l3,q_h3
           q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + half * eighth * &
                (15*q(ilo,jlo-1,k) - 10*q(ilo+1,jlo-1,k) + &
                3*q(ilo+2,jlo-1,k))
        end do
     else
        do k = q_l3,q_h3
           q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + half * half * &
                (3*q(ilo,jlo-1,k) - q(ilo+1,jlo-1,k))
        end do
     end if

     if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
        if (klo+2 .le. ke) then
           q(ilo-1,jlo-1,klo-1) = eighth * ( &
                (15*q(ilo-1,jlo-1,klo) - 10*q(ilo-1,jlo-1,klo+1) + &
                3*q(ilo-1,jlo-1,klo+2)) )
        else
           q(ilo-1,jlo-1,klo-1) = half * &
                (3*q(ilo-1,jlo-1,klo) - q(ilo-1,jlo-1,klo+1))
        end if
     end if

     if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
        if (khi-2 .ge. ks) then
           q(ilo-1,jlo-1,khi+1) = eighth * ( &
                (15*q(ilo-1,jlo-1,khi) - 10*q(ilo-1,jlo-1,khi-1) + &
                3*q(ilo-1,jlo-1,khi-2)) )
        else
           q(ilo-1,jlo-1,khi+1) = half * &
                (3*q(ilo-1,jlo-1,khi) - q(ilo-1,jlo-1,khi-1))
        end if
     end if

  end if
  !
  ! ****************************************************************************
  !
  if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and. &
       (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) ) then
     if (jhi-2 .ge. js) then 
        do k = q_l3,q_h3
           q(ilo-1,jhi+1,k) = half * eighth * &
                (15*q(ilo-1,jhi,k) - 10*q(ilo-1,jhi-1,k) + &
                3*q(ilo-1,jhi-2,k))
        end do
     else
        do k = q_l3,q_h3
           q(ilo-1,jhi+1,k) = half * half * &
                (3*q(ilo-1,jhi,k) - q(ilo-1,jhi-1,k))
        end do
     end if

     if (ilo+2 .le. ie) then 
        do k = q_l3,q_h3
           q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + half * eighth * &
                (15*q(ilo,jhi+1,k) - 10*q(ilo+1,jhi+1,k) + &
                3*q(ilo+2,jhi+1,k))
        end do
     else
        do k = q_l3,q_h3
           q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + half * half * &
                (3*q(ilo,jhi+1,k) - q(ilo+1,jhi+1,k))
        end do
     end if

     if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
        if (klo+2 .le. ke) then
           q(ilo-1,jhi+1,klo-1) = eighth * ( &
                (15*q(ilo-1,jhi+1,klo) - 10*q(ilo-1,jhi+1,klo+1) + &
                3*q(ilo-1,jhi+1,klo+2)) )
        else
           q(ilo-1,jhi+1,klo-1) = half * &
                (3*q(ilo-1,jhi+1,klo) - q(ilo-1,jhi+1,klo+1))
        end if
     end if

     if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
        if (khi-2 .ge. ks) then
           q(ilo-1,jhi+1,khi+1) = eighth * ( &
                (15*q(ilo-1,jhi+1,khi) - 10*q(ilo-1,jhi+1,khi-1) + &
                3*q(ilo-1,jhi+1,khi-2)) )
        else
           q(ilo-1,jhi+1,khi+1) = half * &
                (3*q(ilo-1,jhi+1,khi) - q(ilo-1,jhi+1,khi-1))
        end if
     end if

  end if
  !
  ! ****************************************************************************
  !
  if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and. &
       (nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) ) then
     if (jlo+2 .le. je) then 
        do k = q_l3,q_h3
           q(ihi+1,jlo-1,k) = half * eighth * &
                (15*q(ihi+1,jlo,k) - 10*q(ihi+1,jlo+1,k) + &
                3*q(ihi+1,jlo+2,k))
        end do
     else
        do k = q_l3,q_h3
           q(ihi+1,jlo-1,k) = half * half * &
                (3*q(ihi+1,jlo,k) - q(ihi+1,jlo+1,k))
        end do
     end if

     if (ihi-2 .ge. is) then 
        do k = q_l3,q_h3
           q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + half * eighth * &
                (15*q(ihi,jlo-1,k) - 10*q(ihi-1,jlo-1,k) + &
                3*q(ihi-2,jlo-1,k))
        end do
     else
        do k = q_l3,q_h3
           q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + half * half * &
                (3*q(ihi,jlo-1,k) - q(ihi-1,jlo-1,k))
        end do
     end if

     if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
        if (klo+2 .le. ke) then
           q(ihi+1,jlo-1,klo-1) = eighth * &
                (15*q(ihi+1,jlo-1,klo) - 10*q(ihi+1,jlo-1,klo+1) + &
                3*q(ihi+1,jlo-1,klo+2))
        else
           q(ihi+1,jlo-1,klo-1) = half * &
                (3*q(ihi+1,jlo-1,klo) - q(ihi+1,jlo-1,klo+1))
        end if
     end if

     if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
        if (khi-2 .ge. ks) then
           q(ihi+1,jlo-1,khi+1) = eighth * &
                (15*q(ihi+1,jlo-1,khi) - 10*q(ihi+1,jlo-1,khi-1) + &
                3*q(ihi+1,jlo-1,khi-2))
        else
           q(ihi+1,jlo-1,khi+1) = half * &
                (3*q(ihi+1,jlo-1,khi) - q(ihi+1,jlo-1,khi-1))
        end if
     end if

  end if
  !
  ! ****************************************************************************
  !
  if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and. &
       (ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) ) then
     if (jhi-2 .ge. js) then
        do k = q_l3,q_h3
           q(ihi+1,jhi+1,k) = half * eighth * &
                (15*q(ihi+1,jhi,k) - 10*q(ihi+1,jhi-1,k) + &
                3*q(ihi+1,jhi-2,k))
        end do
     else
        do k = q_l3,q_h3
           q(ihi+1,jhi+1,k) = half * half * &
                (3*q(ihi+1,jhi,k) - q(ihi+1,jhi-1,k))
        end do
     end if

     if (ihi-2 .ge. is) then
        do k = q_l3,q_h3
           q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + half * eighth * &
                (15*q(ihi,jhi+1,k) - 10*q(ihi-1,jhi+1,k) + &
                3*q(ihi-2,jhi+1,k))
        end do
     else
        do k = q_l3,q_h3
           q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + half * half * &
                (3*q(ihi,jhi+1,k) - q(ihi-1,jhi+1,k))
        end do
     end if

     if (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) then
        if (klo+2 .le. ke) then
           q(ihi+1,jhi+1,klo-1) = eighth * &
                (15*q(ihi+1,jhi+1,klo) - 10*q(ihi+1,jhi+1,klo+1) + &
                3*q(ihi+1,jhi+1,klo+2))
        else
           q(ihi+1,jhi+1,klo-1) = half * &
                (3*q(ihi+1,jhi+1,klo) - q(ihi+1,jhi+1,klo+1))
        end if
     end if

     if (nup .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) then
        if (khi-2 .ge. ks) then
           q(ihi+1,jhi+1,khi+1) = eighth * &
                (15*q(ihi+1,jhi+1,khi) - 10*q(ihi+1,jhi+1,khi-1) + &
                3*q(ihi+1,jhi+1,khi-2))
        else
           q(ihi+1,jhi+1,khi+1) = half * &
                (3*q(ihi+1,jhi+1,khi) - q(ihi+1,jhi+1,khi-1))
        end if
     end if

  end if
  !
  !    Next correct the i-k edges
  !
  if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and. &
       (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
     if (klo+2 .le. ke) then
        do j = q_l2,q_h2
           q(ilo-1,j,klo-1) = half * eighth * &
                (15*q(ilo-1,j,klo) - 10*q(ilo-1,j,klo+1) + &
                3*q(ilo-1,j,klo+2))
        end do
     else
        do j = q_l2,q_h2
           q(ilo-1,j,klo-1) = half * half * &
                (3*q(ilo-1,j,klo) - q(ilo-1,j,klo+1))
        end do
     end if

     if (ilo+2 .le. ie) then
        do j = q_l2,q_h2
           q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + half * eighth * &
                (15*q(ilo,j,klo-1) - 10*q(ilo+1,j,klo-1) + &
                3*q(ilo+2,j,klo-1))
        end do
     else
        do j = q_l2,q_h2
           q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + half * half * &
                (3*q(ilo,j,klo-1) - q(ilo+1,j,klo-1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((nlft .gt. 0 .and. bc(1,1) .eq. HOEXTRAP) .and. &
       (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
     if (khi-2 .ge. ks) then
        do j = q_l2,q_h2
           q(ilo-1,j,khi+1) = half * eighth * &
                (15*q(ilo-1,j,khi) - 10*q(ilo-1,j,khi-1) + &
                3*q(ilo-1,j,khi-2))
        end do
     else
        do j = q_l2,q_h2
           q(ilo-1,j,khi+1) = half * half * &
                (3*q(ilo-1,j,khi) - q(ilo-1,j,khi-1))
        end do
     end if

     if (ilo+2 .le. ie) then
        do j = q_l2,q_h2
           q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + half * eighth * &
                (15*q(ilo,j,khi+1) - 10*q(ilo+1,j,khi+1) + &
                3*q(ilo+2,j,khi+1))
        end do
     else
        do j = q_l2,q_h2
           q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + half * half * &
                (3*q(ilo,j,khi+1) - q(ilo+1,j,khi+1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and. &
       (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
     if (klo+2 .le. ke) then
        do j = q_l2,q_h2
           q(ihi+1,j,klo-1) = half * eighth * &
                (15*q(ihi+1,j,klo) - 10*q(ihi+1,j,klo+1) + &
                3*q(ihi+1,j,klo+2))
        end do
     else
        do j = q_l2,q_h2
           q(ihi+1,j,klo-1) = half * half * &
                (3*q(ihi+1,j,klo) - q(ihi+1,j,klo+1))
        end do
     end if

     if (ihi-2 .ge. is) then
        do j = q_l2,q_h2
           q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + half * eighth * &
                (15*q(ihi,j,klo-1) - 10*q(ihi-1,j,klo-1) + &
                3*q(ihi-2,j,klo-1))
        end do
     else
        do j = q_l2,q_h2
           q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + half * half * &
                (3*q(ihi,j,klo-1) - q(ihi-1,j,klo-1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((nrgt .gt. 0 .and. bc(1,2) .eq. HOEXTRAP) .and. &
       (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
     if (khi-2 .ge. ks) then
        do j = q_l2,q_h2
           q(ihi+1,j,khi+1) = half * eighth * &
                (15*q(ihi+1,j,khi) - 10*q(ihi+1,j,khi-1) + &
                3*q(ihi+1,j,khi-2))
        end do
     else
        do j = q_l2,q_h2
           q(ihi+1,j,khi+1) = half * half * &
                (3*q(ihi+1,j,khi) - q(ihi+1,j,khi-1))
        end do
     end if
     if (ihi-2 .ge. is) then
        do j = q_l2,q_h2
           q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + half * eighth * &
                (15*q(ihi,j,khi+1) - 10*q(ihi-1,j,khi+1) + &
                3*q(ihi-2,j,khi+1))
        end do
     else
        do j = q_l2,q_h2
           q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + half * half * &
                (3*q(ihi,j,khi+1) - q(ihi-1,j,khi+1))
        end do
     end if
  end if
  !
  !    Next correct the j-k edges
  !
  if ((nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) .and. &
       (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
     if (klo+2 .le. ke) then
        do i = q_l1,q_h1
           q(i,jlo-1,klo-1) = half * eighth * &
                (15*q(i,jlo-1,klo) - 10*q(i,jlo-1,klo+1) + &
                3*q(i,jlo-1,klo+2))
        end do
     else
        do i = q_l1,q_h1
           q(i,jlo-1,klo-1) = half * half * &
                (3*q(i,jlo-1,klo) - q(i,jlo-1,klo+1))
        end do
     end if
     if (jlo+2 .le. je) then
        do i = q_l1,q_h1
           q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + half * eighth * &
                (15*q(i,jlo,klo-1) - 10*q(i,jlo+1,klo-1) + &
                3*q(i,jlo+2,klo-1))
        end do
     else
        do i = q_l1,q_h1
           q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + half * half * &
                (3*q(i,jlo,klo-1) - q(i,jlo+1,klo-1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((nbot .gt. 0 .and. bc(2,1) .eq. HOEXTRAP) .and. &
       (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
     if (khi-2 .ge. ks) then
        do i = q_l1,q_h1
           q(i,jlo-1,khi+1) = half * eighth * &
                (15*q(i,jlo-1,khi) - 10*q(i,jlo-1,khi-1) + &
                3*q(i,jlo-1,khi-2))
        end do
     else
        do i = q_l1,q_h1
           q(i,jlo-1,khi+1) = half * half * &
                (3*q(i,jlo-1,khi) - q(i,jlo-1,khi-1))
        end do
     end if

     if (jlo+2 .le. je) then
        do i = q_l1,q_h1
           q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + half * eighth * &
                (15*q(i,jlo,khi+1) - 10*q(i,jlo+1,khi+1) + &
                3*q(i,jlo+2,khi+1))
        end do
     else
        do i = q_l1,q_h1
           q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + half * half * &
                (3*q(i,jlo,khi+1) - q(i,jlo+1,khi+1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) .and. &
       (ndwn .gt. 0 .and. bc(3,1) .eq. HOEXTRAP) ) then
     if (klo+2 .le. ke) then
        do i = q_l1,q_h1
           q(i,jhi+1,klo-1) = half * eighth * &
                (15*q(i,jhi+1,klo) - 10*q(i,jhi+1,klo+1) + &
                3*q(i,jhi+1,klo+2))
        end do
     else
        do i = q_l1,q_h1
           q(i,jhi+1,klo-1) = half * half * &
                (3*q(i,jhi+1,klo) - q(i,jhi+1,klo+1))
        end do
     end if
     if (jhi-2 .ge. js) then
        do i = q_l1,q_h1
           q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + half * eighth * &
                (15*q(i,jhi,klo-1) - 10*q(i,jhi-1,klo-1) + &
                3*q(i,jhi-2,klo-1))
        end do
     else
        do i = q_l1,q_h1
           q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + half * half * &
                (3*q(i,jhi,klo-1) - q(i,jhi-1,klo-1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((ntop .gt. 0 .and. bc(2,2) .eq. HOEXTRAP) .and. &
       (nup  .gt. 0 .and. bc(3,2) .eq. HOEXTRAP) ) then
     if (khi-2 .ge. ks) then
        do i = q_l1,q_h1
           q(i,jhi+1,khi+1) = half * eighth * &
                (15*q(i,jhi+1,khi) - 10*q(i,jhi+1,khi-1) + &
                3*q(i,jhi+1,khi-2))
        end do
     else
        do i = q_l1,q_h1
           q(i,jhi+1,khi+1) = half * half * &
                (3*q(i,jhi+1,khi) - q(i,jhi+1,khi-1))
        end do
     end if
     if (jhi-2 .ge. js) then
        do i = q_l1,q_h1
           q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + half * eighth * &
                (15*q(i,jhi,khi+1) - 10*q(i,jhi-1,khi+1) + &
                3*q(i,jhi-2,khi+1))
        end do
     else
        do i = q_l1,q_h1
           q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + half * half * &
                (3*q(i,jhi,khi+1) - q(i,jhi-1,khi+1))
        end do
     end if
  end if

end subroutine filcc



subroutine hoextraptocc(q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,domlo,domhi,dx,xlo)

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in   ) :: q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
  integer,  intent(in   ) :: domlo(3), domhi(3)
  real(rt), intent(in   ) :: xlo(3), dx(3)
  real(rt), intent(inout) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)

  integer :: nlft, nrgt, nbot, ntop, nup, ndwn
  integer :: ilo, ihi, jlo, jhi, klo, khi
  integer :: is,  ie,  js,  je,  ks,  ke
  integer :: i, j, k

  is = max(q_l1,domlo(1))
  ie = min(q_h1,domhi(1))
  js = max(q_l2,domlo(2))
  je = min(q_h2,domhi(2))
  ks = max(q_l3,domlo(3))
  ke = min(q_h3,domhi(3))

  nlft = max(0,domlo(1)-q_l1)
  nrgt = max(0,q_h1-domhi(1))
  nbot = max(0,domlo(2)-q_l2)
  ntop = max(0,q_h2-domhi(2))
  ndwn = max(0,domlo(3)-q_l3)
  nup  = max(0,q_h3-domhi(3))
  !
  !     First fill sides.
  !
  if (nlft .gt. 0) then
     ilo = domlo(1)

     do i = 2, nlft
        do k = q_l3,q_h3
           do j = q_l2,q_h2
              q(ilo-i,j,k) = q(ilo,j,k)
           end do
        end do
     end do
     if (ilo+2 .le. ie) then
        do k = q_l3,q_h3
           do j = q_l2,q_h2
              q(ilo-1,j,k) = 3*q(ilo,j,k) - 3*q(ilo+1,j,k) + &
                   q(ilo+2,j,k)
           end do
        end do
     else  
        do k = q_l3,q_h3
           do j = q_l2,q_h2
              q(ilo-1,j,k) = 2*q(ilo,j,k) - q(ilo+1,j,k)
           end do
        end do
     end if
  end if

  if (nrgt .gt. 0) then
     ihi = domhi(1)

     do i = 2, nrgt
        do k = q_l3,q_h3
           do j = q_l2,q_h2
              q(ihi+i,j,k) = q(ihi,j,k)
           end do
        end do
     end do
     if (ihi-2 .ge. is) then
        do k = q_l3,q_h3
           do j = q_l2,q_h2
              q(ihi+1,j,k) = 3*q(ihi,j,k) - 3*q(ihi-1,j,k) + &
                   q(ihi-2,j,k)
           end do
        end do
     else
        do k = q_l3,q_h3
           do j = q_l2,q_h2
              q(ihi+1,j,k) = 2*q(ihi,j,k) - q(ihi-1,j,k)
           end do
        end do
     end if
  end if

  if (nbot .gt. 0) then
     jlo = domlo(2)

     do j = 2, nbot
        do k = q_l3,q_h3
           do i = q_l1,q_h1
              q(i,jlo-j,k) = q(i,jlo,k)
           end do
        end do
     end do
     if (jlo+2 .le. je) then
        do k = q_l3,q_h3
           do i = q_l1,q_h1
              q(i,jlo-1,k) = 3*q(i,jlo,k) - 3*q(i,jlo+1,k) + &
                   q(i,jlo+2,k)
           end do
        end do
     else
        do k = q_l3,q_h3
           do i = q_l1,q_h1
              q(i,jlo-1,k) = 2*q(i,jlo,k) - q(i,jlo+1,k)
           end do
        end do
     end if
  end if

  if (ntop .gt. 0) then
     jhi = domhi(2)

     do j = 2, ntop
        do k = q_l3,q_h3
           do i = q_l1,q_h1
              q(i,jhi+j,k) = q(i,jhi,k)
           end do
        end do
     end do
     if (jhi-2 .ge. js) then
        do k = q_l3,q_h3
           do i = q_l1,q_h1
              q(i,jhi+1,k) = 3*q(i,jhi,k) - 3*q(i,jhi-1,k) + &
                   q(i,jhi-2,k)
           end do
        end do
     else
        do k = q_l3,q_h3
           do i = q_l1,q_h1
              q(i,jhi+1,k) = 2*q(i,jhi,k) - q(i,jhi-1,k)
           end do
        end do
     end if
  end if

  if (ndwn .gt. 0) then
     klo = domlo(3)

     do k = 2, ndwn
        do j = q_l2,q_h2
           do i = q_l1,q_h1
              q(i,j,klo-k) = q(i,j,klo)
           end do
        end do
     end do
     if (klo+2 .le. ke) then
        do j = q_l2,q_h2
           do i = q_l1,q_h1
              q(i,j,klo-1) = 3*q(i,j,klo) - 3*q(i,j,klo+1) + &
                   q(i,j,klo+2)
           end do
        end do
     else
        do j = q_l2,q_h2
           do i = q_l1,q_h1
              q(i,j,klo-1) = 2*q(i,j,klo) - q(i,j,klo+1)
           end do
        end do
     end if
  end if

  if (nup .gt. 0) then
     khi = domhi(3)

     do k = 2, nup
        do j = q_l2,q_h2
           do i = q_l1,q_h1
              q(i,j,khi+k) = q(i,j,khi)
           end do
        end do
     end do
     if (khi-2 .ge. ks) then
        do j = q_l2,q_h2
           do i = q_l1,q_h1
              q(i,j,khi+1) = 3*q(i,j,khi) - 3*q(i,j,khi-1) + &
                   q(i,j,khi-2)
           end do
        end do
     else
        do j = q_l2,q_h2
           do i = q_l1,q_h1
              q(i,j,khi+1) = 2*q(i,j,khi) - q(i,j,khi-1)
           end do
        end do
     end if
  end if
  !
  !    First correct the i-j edges and all corners
  !
  if ((nlft .gt. 0) .and. (nbot .gt. 0)) then
     if (jlo+2 .le. je) then
        do k = q_l3,q_h3
           q(ilo-1,jlo-1,k) = half * &
                (3*q(ilo-1,jlo,k) - 3*q(ilo-1,jlo+1,k) + &
                q(ilo-1,jlo+2,k))
        end do
     else
        do k = q_l3,q_h3
           q(ilo-1,jlo-1,k) = half * &
                (2*q(ilo-1,jlo,k) - q(ilo-1,jlo+1,k))
        end do
     end if

     if (ilo+2 .le. ie) then
        do k = q_l3,q_h3
           q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + half * &
                (3*q(ilo,jlo-1,k) - 3*q(ilo+1,jlo-1,k) + &
                q(ilo+2,jlo-1,k))
        end do
     else
        do k = q_l3,q_h3
           q(ilo-1,jlo-1,k) = q(ilo-1,jlo-1,k) + half * &
                (2*q(ilo,jlo-1,k) - q(ilo+1,jlo-1,k))
        end do
     end if

     if (ndwn .gt. 0) then
        if (klo+2 .le. ke) then
           q(ilo-1,jlo-1,klo-1) = &
                (3*q(ilo-1,jlo-1,klo) - 3*q(ilo-1,jlo-1,klo+1) + &
                q(ilo-1,jlo-1,klo+2))
        else
           q(ilo-1,jlo-1,klo-1) = &
                (2*q(ilo-1,jlo-1,klo) - q(ilo-1,jlo-1,klo+1))
        end if
     end if

     if (nup .gt. 0) then
        if (khi-2 .ge. ks) then
           q(ilo-1,jlo-1,khi+1) = &
                (3*q(ilo-1,jlo-1,khi) - 3*q(ilo-1,jlo-1,khi-1) + &
                q(ilo-1,jlo-1,khi-2))
        else
           q(ilo-1,jlo-1,khi+1) = &
                (2*q(ilo-1,jlo-1,khi) - q(ilo-1,jlo-1,khi-1))
        end if
     end if

  end if
  !
  ! ****************************************************************************
  !
  if ((nlft .gt. 0) .and. (ntop .gt. 0)) then
     if (jhi-2 .ge. js) then 
        do k = q_l3,q_h3
           q(ilo-1,jhi+1,k) = half * &
                (3*q(ilo-1,jhi,k) - 3*q(ilo-1,jhi-1,k) + &
                q(ilo-1,jhi-2,k))
        end do
     else
        do k = q_l3,q_h3
           q(ilo-1,jhi+1,k) = half * &
                (2*q(ilo-1,jhi,k) - q(ilo-1,jhi-1,k))
        end do
     end if

     if (ilo+2 .le. ie) then 
        do k = q_l3,q_h3
           q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + half * &
                (3*q(ilo,jhi+1,k) - 3*q(ilo+1,jhi+1,k) + &
                q(ilo+2,jhi+1,k))
        end do
     else
        do k = q_l3,q_h3
           q(ilo-1,jhi+1,k) = q(ilo-1,jhi+1,k) + half * &
                (2*q(ilo,jhi+1,k) - q(ilo+1,jhi+1,k))
        end do
     end if

     if (ndwn .gt. 0) then
        if (klo+2 .le. ke) then
           q(ilo-1,jhi+1,klo-1) = &
                (3*q(ilo-1,jhi+1,klo) - 3*q(ilo-1,jhi+1,klo+1) + &
                q(ilo-1,jhi+1,klo+2))
        else
           q(ilo-1,jhi+1,klo-1) = &
                (2*q(ilo-1,jhi+1,klo) - q(ilo-1,jhi+1,klo+1))
        end if
     end if

     if (nup .gt. 0) then
        if (khi-2 .ge. ks) then
           q(ilo-1,jhi+1,khi+1) = &
                (3*q(ilo-1,jhi+1,khi) - 3*q(ilo-1,jhi+1,khi-1) + &
                q(ilo-1,jhi+1,khi-2))
        else
           q(ilo-1,jhi+1,khi+1) = &
                (2*q(ilo-1,jhi+1,khi) - q(ilo-1,jhi+1,khi-1))
        end if
     end if

  end if
  !
  ! ****************************************************************************
  !
  if ((nrgt .gt. 0) .and. (nbot .gt. 0)) then
     if (jlo+2 .le. je) then 
        do k = q_l3,q_h3
           q(ihi+1,jlo-1,k) = half * &
                (3*q(ihi+1,jlo,k) - 3*q(ihi+1,jlo+1,k) + &
                q(ihi+1,jlo+2,k))
        end do
     else
        do k = q_l3,q_h3
           q(ihi+1,jlo-1,k) = half * &
                (2*q(ihi+1,jlo,k) - q(ihi+1,jlo+1,k))
        end do
     end if

     if (ihi-2 .ge. is) then 
        do k = q_l3,q_h3
           q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + half * &
                (3*q(ihi,jlo-1,k) - 3*q(ihi-1,jlo-1,k) + &
                q(ihi-2,jlo-1,k))
        end do
     else
        do k = q_l3,q_h3
           q(ihi+1,jlo-1,k) = q(ihi+1,jlo-1,k) + half * &
                (2*q(ihi,jlo-1,k) - q(ihi-1,jlo-1,k))
        end do
     end if

     if (ndwn .gt. 0) then
        if (klo+2 .le. ke) then
           q(ihi+1,jlo-1,klo-1) = &
                (3*q(ihi+1,jlo-1,klo) - 3*q(ihi+1,jlo-1,klo+1) + &
                q(ihi+1,jlo-1,klo+2))
        else
           q(ihi+1,jlo-1,klo-1) = &
                (2*q(ihi+1,jlo-1,klo) - q(ihi+1,jlo-1,klo+1))
        end if
     end if

     if (nup .gt. 0) then
        if (khi-2 .ge. ks) then
           q(ihi+1,jlo-1,khi+1) = &
                (3*q(ihi+1,jlo-1,khi) - 3*q(ihi+1,jlo-1,khi-1) + &
                q(ihi+1,jlo-1,khi-2))
        else
           q(ihi+1,jlo-1,khi+1) = &
                (2*q(ihi+1,jlo-1,khi) - q(ihi+1,jlo-1,khi-1))
        end if
     end if

  end if
  !
  ! ****************************************************************************
  !
  if ((nrgt .gt. 0) .and. (ntop .gt. 0)) then
     if (jhi-2 .ge. js) then
        do k = q_l3,q_h3
           q(ihi+1,jhi+1,k) = half * &
                (3*q(ihi+1,jhi,k) - 3*q(ihi+1,jhi-1,k) + &
                q(ihi+1,jhi-2,k))
        end do
     else
        do k = q_l3,q_h3
           q(ihi+1,jhi+1,k) = half * &
                (2*q(ihi+1,jhi,k) - q(ihi+1,jhi-1,k))
        end do
     end if

     if (ihi-2 .ge. is) then
        do k = q_l3,q_h3
           q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + half * &
                (3*q(ihi,jhi+1,k) - 3*q(ihi-1,jhi+1,k) + &
                q(ihi-2,jhi+1,k))
        end do
     else
        do k = q_l3,q_h3
           q(ihi+1,jhi+1,k) = q(ihi+1,jhi+1,k) + half * &
                (2*q(ihi,jhi+1,k) - q(ihi-1,jhi+1,k))
        end do
     end if

     if (ndwn .gt. 0) then
        if (klo+2 .le. ke) then
           q(ihi+1,jhi+1,klo-1) = &
                (3*q(ihi+1,jhi+1,klo) - 3*q(ihi+1,jhi+1,klo+1) + &
                q(ihi+1,jhi+1,klo+2))
        else
           q(ihi+1,jhi+1,klo-1) = &
                (2*q(ihi+1,jhi+1,klo) - q(ihi+1,jhi+1,klo+1))
        end if
     end if

     if (nup .gt. 0) then
        if (khi-2 .ge. ks) then
           q(ihi+1,jhi+1,khi+1) = &
                (3*q(ihi+1,jhi+1,khi) - 3*q(ihi+1,jhi+1,khi-1) + &
                q(ihi+1,jhi+1,khi-2))
        else
           q(ihi+1,jhi+1,khi+1) = &
                (2*q(ihi+1,jhi+1,khi) - q(ihi+1,jhi+1,khi-1))
        end if
     end if

  end if
  !
  !    Next correct the i-k edges
  !
  if ((nlft .gt. 0) .and. (ndwn .gt. 0)) then
     if (klo+2 .le. ke) then
        do j = q_l2,q_h2
           q(ilo-1,j,klo-1) = half * &
                (3*q(ilo-1,j,klo) - 3*q(ilo-1,j,klo+1) + &
                q(ilo-1,j,klo+2))
        end do
     else
        do j = q_l2,q_h2
           q(ilo-1,j,klo-1) = half * &
                (2*q(ilo-1,j,klo) - q(ilo-1,j,klo+1))
        end do
     end if

     if (ilo+2 .le. ie) then
        do j = q_l2,q_h2
           q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + half * &
                (3*q(ilo,j,klo-1) - 3*q(ilo+1,j,klo-1) + &
                q(ilo+2,j,klo-1))
        end do
     else
        do j = q_l2,q_h2
           q(ilo-1,j,klo-1) = q(ilo-1,j,klo-1) + half * &
                (2*q(ilo,j,klo-1) - q(ilo+1,j,klo-1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((nlft .gt. 0) .and. (nup .gt. 0)) then
     if (khi-2 .ge. ks) then
        do j = q_l2,q_h2
           q(ilo-1,j,khi+1) = half * &
                (3*q(ilo-1,j,khi) - 3*q(ilo-1,j,khi-1) + &
                q(ilo-1,j,khi-2))
        end do
     else
        do j = q_l2,q_h2
           q(ilo-1,j,khi+1) = half * &
                (2*q(ilo-1,j,khi) - q(ilo-1,j,khi-1))
        end do
     end if

     if (ilo+2 .le. ie) then
        do j = q_l2,q_h2
           q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + half * &
                (3*q(ilo,j,khi+1) - 3*q(ilo+1,j,khi+1) + &
                q(ilo+2,j,khi+1))
        end do
     else
        do j = q_l2,q_h2
           q(ilo-1,j,khi+1) = q(ilo-1,j,khi+1) + half * &
                (2*q(ilo,j,khi+1) - q(ilo+1,j,khi+1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((nrgt .gt. 0) .and. (ndwn .gt. 0)) then
     if (klo+2 .le. ke) then
        do j = q_l2,q_h2
           q(ihi+1,j,klo-1) = half * &
                (3*q(ihi+1,j,klo) - 3*q(ihi+1,j,klo+1) + &
                q(ihi+1,j,klo+2))
        end do
     else
        do j = q_l2,q_h2
           q(ihi+1,j,klo-1) = half * &
                (2*q(ihi+1,j,klo) - q(ihi+1,j,klo+1))
        end do
     end if

     if (ihi-2 .ge. is) then
        do j = q_l2,q_h2
           q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + half * &
                (3*q(ihi,j,klo-1) - 3*q(ihi-1,j,klo-1) + &
                q(ihi-2,j,klo-1))
        end do
     else
        do j = q_l2,q_h2
           q(ihi+1,j,klo-1) = q(ihi+1,j,klo-1) + half * &
                (2*q(ihi,j,klo-1) - q(ihi-1,j,klo-1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((nrgt .gt. 0) .and. (nup .gt. 0)) then
     if (khi-2 .ge. ks) then
        do j = q_l2,q_h2
           q(ihi+1,j,khi+1) = half * &
                (3*q(ihi+1,j,khi) - 3*q(ihi+1,j,khi-1) + &
                q(ihi+1,j,khi-2))
        end do
     else
        do j = q_l2,q_h2
           q(ihi+1,j,khi+1) = half * &
                (2*q(ihi+1,j,khi) - q(ihi+1,j,khi-1))
        end do
     end if
     if (ihi-2 .ge. is) then
        do j = q_l2,q_h2
           q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + half * &
                (3*q(ihi,j,khi+1) - 3*q(ihi-1,j,khi+1) + &
                q(ihi-2,j,khi+1))
        end do
     else
        do j = q_l2,q_h2
           q(ihi+1,j,khi+1) = q(ihi+1,j,khi+1) + half * &
                (2*q(ihi,j,khi+1) - q(ihi-1,j,khi+1))
        end do
     end if
  end if
  !
  !    Next correct the j-k edges
  !
  if ((nbot .gt. 0) .and. (ndwn .gt. 0)) then
     if (klo+2 .le. ke) then
        do i = q_l1,q_h1
           q(i,jlo-1,klo-1) = half * &
                (3*q(i,jlo-1,klo) - 3*q(i,jlo-1,klo+1) + &
                q(i,jlo-1,klo+2))
        end do
     else
        do i = q_l1,q_h1
           q(i,jlo-1,klo-1) = half * &
                (2*q(i,jlo-1,klo) - q(i,jlo-1,klo+1))
        end do
     end if
     if (jlo+2 .le. je) then
        do i = q_l1,q_h1
           q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + half * &
                (3*q(i,jlo,klo-1) - 3*q(i,jlo+1,klo-1) + &
                q(i,jlo+2,klo-1))
        end do
     else
        do i = q_l1,q_h1
           q(i,jlo-1,klo-1) = q(i,jlo-1,klo-1) + half * &
                (2*q(i,jlo,klo-1) - q(i,jlo+1,klo-1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((nbot .gt. 0) .and. (nup .gt. 0)) then
     if (khi-2 .ge. ks) then
        do i = q_l1,q_h1
           q(i,jlo-1,khi+1) = half * &
                (3*q(i,jlo-1,khi) - 3*q(i,jlo-1,khi-1) + &
                q(i,jlo-1,khi-2))
        end do
     else
        do i = q_l1,q_h1
           q(i,jlo-1,khi+1) = half * &
                (2*q(i,jlo-1,khi) - q(i,jlo-1,khi-1))
        end do
     end if

     if (jlo+2 .le. je) then
        do i = q_l1,q_h1
           q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + half * &
                (3*q(i,jlo,khi+1) - 3*q(i,jlo+1,khi+1) + &
                q(i,jlo+2,khi+1))
        end do
     else
        do i = q_l1,q_h1
           q(i,jlo-1,khi+1) = q(i,jlo-1,khi+1) + half * &
                (2*q(i,jlo,khi+1) - q(i,jlo+1,khi+1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((ntop .gt. 0) .and. (ndwn .gt. 0)) then
     if (klo+2 .le. ke) then
        do i = q_l1,q_h1
           q(i,jhi+1,klo-1) = half * &
                (3*q(i,jhi+1,klo) - 3*q(i,jhi+1,klo+1) + &
                q(i,jhi+1,klo+2))
        end do
     else
        do i = q_l1,q_h1
           q(i,jhi+1,klo-1) = half * &
                (2*q(i,jhi+1,klo) - q(i,jhi+1,klo+1))
        end do
     end if
     if (jhi-2 .ge. js) then
        do i = q_l1,q_h1
           q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + half * &
                (3*q(i,jhi,klo-1) - 3*q(i,jhi-1,klo-1) + &
                q(i,jhi-2,klo-1))
        end do
     else
        do i = q_l1,q_h1
           q(i,jhi+1,klo-1) = q(i,jhi+1,klo-1) + half * &
                (2*q(i,jhi,klo-1) - q(i,jhi-1,klo-1))
        end do
     end if
  end if
  !
  ! ****************************************************************************
  !
  if ((ntop .gt. 0) .and. (nup .gt. 0)) then
     if (khi-2 .ge. ks) then
        do i = q_l1,q_h1
           q(i,jhi+1,khi+1) = half * &
                (3*q(i,jhi+1,khi) - 3*q(i,jhi+1,khi-1) + &
                q(i,jhi+1,khi-2))
        end do
     else
        do i = q_l1,q_h1
           q(i,jhi+1,khi+1) = half * &
                (2*q(i,jhi+1,khi) - q(i,jhi+1,khi-1))
        end do
     end if
     if (jhi-2 .ge. js) then
        do i = q_l1,q_h1
           q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + half * &
                (3*q(i,jhi,khi+1) - 3*q(i,jhi-1,khi+1) + &
                q(i,jhi-2,khi+1))
        end do
     else
        do i = q_l1,q_h1
           q(i,jhi+1,khi+1) = q(i,jhi+1,khi+1) + half * &
                (2*q(i,jhi,khi+1) - q(i,jhi-1,khi+1))
        end do
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
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: xlo(3), dx(3)
    real(rt), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),ncomp)
    integer,  intent(in   ) :: bc(3,2,ncomp)
    integer,  intent(in   ) :: ncomp

    integer :: ilo, ihi, jlo, jhi, klo, khi
    integer :: is, ie, js, je, ks, ke
    integer :: i, j, k, n

    is = max(q_lo(1), domlo(1))
    ie = min(q_hi(1), domhi(1))
    js = max(q_lo(2), domlo(2))
    je = min(q_hi(2), domhi(2))
    ks = max(q_lo(3), domlo(3))
    ke = min(q_hi(3), domhi(3))

    ilo = domlo(1)
    ihi = domhi(1)
    jlo = domlo(2)
    jhi = domhi(2)
    klo = domlo(3)
    khi = domhi(3)

    do n = 1, ncomp
       do k = blo(3), bhi(3)
          do j = blo(2), bhi(2)
             do i = blo(1), bhi(1)

                if (bc(1,1,n) .eq. EXT_DIR) then

                   ! Do nothing.

                else if (bc(1,1,n) .eq. FOEXTRAP) then

                   if (i < ilo) then
                      q(i,j,k,n) = q(ilo,j,k,n)
                   end if

                else if (bc(1,1,n) .eq. HOEXTRAP) then

                   if (i < ilo - 1) then
                      q(i,j,k,n) = q(ilo,j,k,n)
                   else if (i == ilo - 1) then
                      if (ilo+2 <= ie) then
                         q(i,j,k,n) = eighth * (15*q(ilo,j,k,n) - 10*q(ilo+1,j,k,n) + 3*q(ilo+2,j,k,n))
                      else
                         q(i,j,k,n) = half * (3*q(ilo,j,k,n) - q(ilo+1,j,k,n))
                      end if
                   end if

                else if (bc(1,1,n) .eq. REFLECT_EVEN) then

                   if (i < ilo) then
                      q(i,j,k,n) = q(ilo+(ilo-i)-1,j,k,n)
                   end if

                else if (bc(1,1,n) .eq. REFLECT_ODD) then

                   if (i < ilo) then
                      q(i,j,k,n) = -q(ilo+(ilo-i)-1,j,k,n)
                   end if

                end if



                if (bc(1,2,n) .eq. EXT_DIR) then

                   ! Do nothing.

                else if (bc(1,2,n) .eq. FOEXTRAP) then

                   if (i > ihi) then
                      q(i,j,k,n) = q(ihi,j,k,n)
                   end if

                else if (bc(1,2,n) .eq. HOEXTRAP) then

                   if (i > ihi + 1) then
                      q(i,j,k,n) = q(ihi,j,k,n)
                   else if (i == ihi + 1) then
                      if (ihi-2 >= is) then
                         q(i,j,k,n) = eighth * (15*q(ihi,j,k,n) - 10*q(ihi-1,j,k,n) + 3*q(ihi-2,j,k,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi,j,k,n) - q(ihi-1,j,k,n))
                      end if
                   end if

                else if (bc(1,2,n) .eq. REFLECT_EVEN) then

                   if (i > ihi) then
                      q(i,j,k,n) = q(ihi-(i-ihi)+1,j,k,n)
                   end if

                else if (bc(1,2,n) .eq. REFLECT_ODD) then

                   if (i > ihi) then
                      q(i,j,k,n) = -q(ihi-(i-ihi)+1,j,k,n)
                   end if

                end if



                if (bc(2,1,n) .eq. EXT_DIR) then

                   ! Do nothing.

                else if (bc(2,1,n) .eq. FOEXTRAP) then

                   if (j < jlo) then
                      q(i,j,k,n) = q(i,jlo,k,n)
                   end if

                else if (bc(2,1,n) .eq. HOEXTRAP) then

                   if (j < jlo - 1) then
                      q(i,j,k,n) = q(i,jlo,k,n)
                   else if (j == jlo - 1) then
                      if (jlo+2 <= je) then
                         q(i,j,k,n) = eighth * (15*q(i,jlo,k,n) - 10*q(i,jlo+1,k,n) + 3*q(i,jlo+2,k,n))
                      else
                         q(i,j,k,n) = half * (3*q(i,jlo,k,n) - q(i,jlo+1,k,n))
                      end if
                   end if

                else if (bc(2,1,n) .eq. REFLECT_EVEN) then

                   if (j < jlo) then
                      q(i,j,k,n) = q(i,jlo+(jlo-j)-1,k,n)
                   end if

                else if (bc(2,1,n) .eq. REFLECT_ODD) then

                   if (j < jlo) then
                      q(i,j,k,n) = -q(i,jlo+(jlo-j)-1,k,n)
                   end if

                end if



                if (bc(2,2,n) .eq. EXT_DIR) then

                   ! Do nothing.

                else if (bc(2,2,n) .eq. FOEXTRAP) then

                   if (j > jhi) then
                      q(i,j,k,n) = q(i,jhi,k,n)
                   end if

                else if (bc(2,2,n) .eq. HOEXTRAP) then

                   if (j > jhi + 1) then
                      q(i,j,k,n) = q(i,jhi,k,n)
                   else if (j == jhi + 1) then
                      if (jhi-2 >= js) then
                         q(i,j,k,n) = eighth * (15*q(i,jhi,k,n) - 10*q(i,jhi-1,k,n) + 3*q(i,jhi-2,k,n))
                      else
                         q(i,j,k,n) = half * (3*q(i,jhi,k,n) - q(i,jhi-1,k,n))
                      end if
                   end if

                else if (bc(2,2,n) .eq. REFLECT_EVEN) then

                   if (j > jhi) then
                      q(i,j,k,n) = q(i,jhi-(j-jhi)+1,k,n)
                   end if

                else if (bc(2,2,n) .eq. REFLECT_ODD) then

                   if (j > jhi) then
                      q(i,j,k,n) = -q(i,jhi-(j-jhi)+1,k,n)
                   end if

                end if



                if (bc(3,1,n) .eq. EXT_DIR) then

                   ! Do nothing.

                else if (bc(3,1,n) .eq. FOEXTRAP) then

                   if (k < klo) then
                      q(i,j,k,n) = q(i,j,klo,n)
                   end if

                else if (bc(3,1,n) .eq. HOEXTRAP) then

                   if (k < klo - 1) then
                      q(i,j,k,n) = q(i,j,klo,n)
                   else if (k == klo - 1) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * (15*q(i,j,klo,n) - 10*q(i,j,klo+1,n) + 3*q(i,j,klo+2,n))
                      else
                         q(i,j,k,n) = half * (3*q(i,j,klo,n) - q(i,j,klo+1,n))
                      end if
                   end if

                else if (bc(3,1,n) .eq. REFLECT_EVEN) then

                   if (k < klo) then
                      q(i,j,k,n) = q(i,j,klo+(klo-k)-1,n)
                   end if

                else if (bc(3,1,n) .eq. REFLECT_ODD) then

                   if (k < klo) then
                      q(i,j,k,n) = -q(i,j,klo+(klo-k)-1,n)
                   end if

                end if



                if (bc(3,2,n) .eq. EXT_DIR) then

                   ! Do nothing.

                else if (bc(3,2,n) .eq. FOEXTRAP) then

                   if (k > khi) then
                      q(i,j,k,n) = q(i,j,khi,n)
                   end if

                else if (bc(3,2,n) .eq. HOEXTRAP) then

                   if (k > khi + 1) then
                      q(i,j,k,n) = q(i,j,khi,n)
                   else if (k == khi + 1) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * (15*q(i,j,khi,n) - 10*q(i,j,khi-1,n) + 3*q(i,j,khi-2,n))
                      else
                         q(i,j,k,n) = half * (3*q(i,j,khi,n) - q(i,j,khi-1,n))
                      end if
                   end if

                else if (bc(3,2,n) .eq. REFLECT_EVEN) then

                   if (k > khi) then
                      q(i,j,k,n) = q(i,j,khi-(k-khi)+1,n)
                   end if

                else if (bc(3,2,n) .eq. REFLECT_ODD) then

                   if (k > khi) then
                      q(i,j,k,n) = -q(i,j,khi-(k-khi)+1,n)
                   end if

                end if



                !
                ! First correct the i-j edges and all corners
                !

                if (i == ilo-1 .and. bc(1,1,n) .eq. HOEXTRAP .and. &
                    j == jlo-1 .and. bc(2,1,n) .eq. HOEXTRAP) then

                   if (jlo+2 <= je) then
                      q(i,j,k,n) = half * eighth * (15*q(ilo-1,jlo,k,n) - 10*q(ilo-1,jlo+1,k,n) + 3*q(ilo-1,jlo+2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ilo-1,jlo,k,n) - q(ilo-1,jlo+1,k,n))
                   end if

                   if (ilo+2 <= ie) then
                      q(i,j,k,n) = q(ilo-1,jlo-1,k,n) + &
                                   half * eighth * (15*q(ilo,jlo-1,k,n) - 10*q(ilo+1,jlo-1,k,n) + 3*q(ilo+2,jlo-1,k,n))
                   else
                      q(i,j,k,n) = q(ilo-1,jlo-1,k,n) + half * half * (3*q(ilo,jlo-1,k,n) - q(ilo+1,jlo-1,k,n))
                   end if

                   if (k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jlo-1,klo,n) - 10*q(ilo-1,jlo-1,klo+1,n) + &
                                                  3*q(ilo-1,jlo-1,klo+2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jlo-1,klo,n) - q(ilo-1,jlo-1,klo+1,n))
                      end if
                   end if

                   if (k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jlo-1,khi,n) - 10*q(ilo-1,jlo-1,khi-1,n) + &
                                                  3*q(ilo-1,jlo-1,khi-2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jlo-1,khi,n) - q(ilo-1,jlo-1,khi-1,n))
                      end if
                   end if

                end if

                !
                ! ****************************************************************************
                !

                if (i == ilo-1 .and. bc(1,1,n) .eq. HOEXTRAP .and. &
                    j == jhi+1 .and. bc(2,2,n) .eq. HOEXTRAP) then

                   if (jhi-2 >= js) then
                      q(i,j,k,n) = half * eighth * (15*q(ilo-1,jhi,k,n) - 10*q(ilo-1,jhi-1,k,n) + 3*q(ilo-1,jhi-2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ilo-1,jhi,k,n) - q(ilo-1,jhi-1,k,n))
                   end if

                   if (ilo+2 <= ie) then
                      q(i,j,k,n) = q(ilo-1,jhi+1,k,n) + &
                                   half * eighth * (15*q(ilo,jhi+1,k,n) - 10*q(ilo+1,jhi+1,k,n) + 3*q(ilo+2,jhi+1,k,n))
                   else
                      q(i,j,k,n) = q(ilo-1,jhi+1,k,n) + half * half * (3*q(ilo,jhi+1,k,n) - q(ilo+1,jhi+1,k,n))
                   end if

                   if (k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jhi+1,klo,n) - 10*q(ilo-1,jhi+1,klo+1,n) + &
                                                  3*q(ilo-1,jhi+1,klo+2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jhi+1,klo,n) - q(ilo-1,jhi+1,klo+1,n))
                      end if
                   end if

                   if (k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jhi+1,khi,n) - 10*q(ilo-1,jhi+1,khi-1,n) + &
                                                  3*q(ilo-1,jhi+1,khi-2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jhi+1,khi,n) - q(ilo-1,jhi+1,khi-1,n))
                      end if
                   end if

                end if

                !
                ! ****************************************************************************
                !

                if (i == ihi+1 .and. bc(1,2,n) .eq. HOEXTRAP .and. &
                    j == jlo-1 .and. bc(2,1,n) .eq. HOEXTRAP) then

                   if (jlo+2 <= je) then
                      q(i,j,k,n) = half * eighth * (15*q(ihi+1,jlo,k,n) - 10*q(ihi+1,jlo+1,k,n) + 3*q(ihi+1,jlo+2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ihi+1,jlo,k,n) - q(ihi+1,jlo+1,k,n))
                   end if

                   if (ihi-2 >= is) then
                      q(i,j,k,n) = q(ihi+1,jlo-1,k,n) + &
                                   half * eighth * (15*q(ihi,jlo-1,k,n) - 10*q(ihi-1,jlo-1,k,n) + 3*q(ihi-2,jlo-1,k,n))
                   else
                      q(i,j,k,n) = q(ihi+1,jlo-1,k,n) + half * half * (3*q(ihi,jlo-1,k,n) - q(ihi-1,jlo-1,k,n))
                   end if

                   if (k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jlo-1,klo,n) - 10*q(ihi+1,jlo-1,klo+1,n) + 3*q(ihi+1,jlo-1,klo+2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jlo-1,klo,n) - q(ihi+1,jlo-1,klo+1,n))
                      end if
                   end if

                   if (k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jlo-1,khi,n) - 10*q(ihi+1,jlo-1,khi-1,n) + 3*q(ihi+1,jlo-1,khi-2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jlo-1,khi,n) - q(ihi+1,jlo-1,khi-1,n))
                      end if
                   end if

                end if

                !
                ! ****************************************************************************
                !

                if (i == ihi+1 .and. bc(1,2,n) .eq. HOEXTRAP .and. &
                    j == jhi+1 .and. bc(2,2,n) .eq. HOEXTRAP) then

                   if (jhi-2 >= js) then
                      q(i,j,k,n) = half * eighth * (15*q(ihi+1,jhi,k,n) - 10*q(ihi+1,jhi-1,k,n) + 3*q(ihi+1,jhi-2,k,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ihi+1,jhi,k,n) - q(ihi+1,jhi-1,k,n))
                   end if

                   if (ihi-2 >= is) then
                      q(i,j,k,n) = q(ihi+1,jhi+1,k,n) + &
                                   half * eighth * (15*q(ihi,jhi+1,k,n) - 10*q(ihi-1,jhi+1,k,n) + 3*q(ihi-2,jhi+1,k,n))
                   else
                      q(i,j,k,n) = q(ihi+1,jhi+1,k,n) + half * half * (3*q(ihi,jhi+1,k,n) - q(ihi-1,jhi+1,k,n))
                   end if

                   if (k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jhi+1,klo,n) - 10*q(ihi+1,jhi+1,klo+1,n) + 3*q(ihi+1,jhi+1,klo+2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jhi+1,klo,n) - q(ihi+1,jhi+1,klo+1,n))
                      end if
                   end if

                   if (k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jhi+1,khi,n) - 10*q(ihi+1,jhi+1,khi-1,n) + 3*q(ihi+1,jhi+1,khi-2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jhi+1,khi,n) - q(ihi+1,jhi+1,khi-1,n))
                      end if
                   end if

                end if

                !
                ! Next correct the i-k edges
                !

                if (i == ilo-1 .and. bc(1,1,n) .eq. HOEXTRAP .and. &
                    k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then

                   if (klo+2 <= ke) then
                      q(i,j,k,n) = half * eighth * (15*q(ilo-1,j,klo,n) - 10*q(ilo-1,j,klo+1,n) + 3*q(ilo-1,j,klo+2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ilo-1,j,klo,n) - q(ilo-1,j,klo+1,n))
                   end if

                   if (ilo+2 <= ie) then
                      q(i,j,k,n) = q(ilo-1,j,klo-1,n) + &
                                   half * eighth * (15*q(ilo,j,klo-1,n) - 10*q(ilo+1,j,klo-1,n) + 3*q(ilo+2,j,klo-1,n))
                   else
                      q(i,j,k,n) = q(ilo-1,j,klo-1,n) + half * half * (3*q(ilo,j,klo-1,n) - q(ilo+1,j,klo-1,n))
                   end if

                end if

                !
                ! ****************************************************************************
                !

                if (i == ilo-1 .and. bc(1,1,n) .eq. HOEXTRAP .and. &
                    k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then

                   if (khi-2 >= ks) then
                      q(i,j,k,n) = half * eighth * (15*q(ilo-1,j,khi,n) - 10*q(ilo-1,j,khi-1,n) + 3*q(ilo-1,j,khi-2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ilo-1,j,khi,n) - q(ilo-1,j,khi-1,n))
                   end if

                   if (ilo+2 <= ie) then
                      q(i,j,k,n) = q(ilo-1,j,khi+1,n) + &
                                   half * eighth * (15*q(ilo,j,khi+1,n) - 10*q(ilo+1,j,khi+1,n) + 3*q(ilo+2,j,khi+1,n))
                   else
                      q(i,j,k,n) = q(ilo-1,j,khi+1,n) + half * half * (3*q(ilo,j,khi+1,n) - q(ilo+1,j,khi+1,n))
                   end if

                end if

                !
                ! ****************************************************************************
                !

                if (i == ihi+1 .and. bc(1,2,n) .eq. HOEXTRAP .and. &
                    k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then

                   if (klo+2 <= ke) then
                      q(i,j,k,n) = half * eighth * (15*q(ihi+1,j,klo,n) - 10*q(ihi+1,j,klo+1,n) + 3*q(ihi+1,j,klo+2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ihi+1,j,klo,n) - q(ihi+1,j,klo+1,n))
                   end if

                   if (ihi-2 >= is) then
                      q(i,j,k,n) = q(ihi+1,j,klo-1,n) + &
                                 half * eighth * (15*q(ihi,j,klo-1,n) - 10*q(ihi-1,j,klo-1,n) + 3*q(ihi-2,j,klo-1,n))
                   else
                      q(i,j,k,n) = q(ihi+1,j,klo-1,n) + half * half * (3*q(ihi,j,klo-1,n) - q(ihi-1,j,klo-1,n))
                   end if

                end if

                !
                ! ****************************************************************************
                !
                if (i == ihi+1 .and. bc(1,2,n) .eq. HOEXTRAP .and. &
                    k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then

                   if (khi-2 >= ks) then
                      q(i,j,k,n) = half * eighth * (15*q(ihi+1,j,khi,n) - 10*q(ihi+1,j,khi-1,n) + 3*q(ihi+1,j,khi-2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(ihi+1,j,khi,n) - q(ihi+1,j,khi-1,n))
                   end if

                   if (ihi-2 >= is) then
                      q(i,j,k,n) = q(ihi+1,j,khi+1,n) + &
                                   half * eighth * (15*q(ihi,j,khi+1,n) - 10*q(ihi-1,j,khi+1,n) + 3*q(ihi-2,j,khi+1,n))
                   else
                      q(i,j,k,n) = q(ihi+1,j,khi+1,n) + half * half * (3*q(ihi,j,khi+1,n) - q(ihi-1,j,khi+1,n))
                   end if

                end if

                !
                ! Next correct the j-k edges
                !

                if (j == jlo-1 .and. bc(2,1,n) .eq. HOEXTRAP .and. &
                    k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then

                   if (klo+2 <= ke) then
                      q(i,j,k,n) = half * eighth * (15*q(i,jlo-1,klo,n) - 10*q(i,jlo-1,klo+1,n) + 3*q(i,jlo-1,klo+2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(i,jlo-1,klo,n) - q(i,jlo-1,klo+1,n))
                   end if

                   if (jlo+2 <= je) then
                      q(i,j,k,n) = q(i,jlo-1,klo-1,n) + &
                                   half * eighth * (15*q(i,jlo,klo-1,n) - 10*q(i,jlo+1,klo-1,n) + 3*q(i,jlo+2,klo-1,n))
                   else
                      q(i,j,k,n) = q(i,jlo-1,klo-1,n) + half * half * (3*q(i,jlo,klo-1,n) - q(i,jlo+1,klo-1,n))
                   end if

                end if

                !
                ! ****************************************************************************
                !

                if (j == jlo-1 .and. bc(2,1,n) .eq. HOEXTRAP .and. &
                    k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then

                   if (khi-2 >= ks) then
                      q(i,j,k,n) = half * eighth * (15*q(i,jlo-1,khi,n) - 10*q(i,jlo-1,khi-1,n) + 3*q(i,jlo-1,khi-2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(i,jlo-1,khi,n) - q(i,jlo-1,khi-1,n))
                   end if

                   if (jlo+2 <= je) then
                      q(i,j,k,n) = q(i,jlo-1,khi+1,n) + &
                                    half * eighth * (15*q(i,jlo,khi+1,n) - 10*q(i,jlo+1,khi+1,n) + 3*q(i,jlo+2,khi+1,n))
                   else
                      q(i,j,k,n) = q(i,jlo-1,khi+1,n) + half * half * (3*q(i,jlo,khi+1,n) - q(i,jlo+1,khi+1,n))
                   end if

                end if

                !
                ! ****************************************************************************
                !

                if (j == jhi+1 .and. bc(2,2,n) .eq. HOEXTRAP .and. &
                    k == klo-1 .and. bc(3,1,n) .eq. HOEXTRAP) then

                   if (klo+2 <= ke) then
                      q(i,j,k,n) = half * eighth * (15*q(i,jhi+1,klo,n) - 10*q(i,jhi+1,klo+1,n) + 3*q(i,jhi+1,klo+2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(i,jhi+1,klo,n) - q(i,jhi+1,klo+1,n))
                   end if

                   if (jhi-2 >= js) then
                      q(i,j,k,n) = q(i,jhi+1,klo-1,n) + &
                                   half * eighth * (15*q(i,jhi,klo-1,n) - 10*q(i,jhi-1,klo-1,n) + 3*q(i,jhi-2,klo-1,n))
                   else
                      q(i,j,k,n) = q(i,jhi+1,klo-1,n) + half * half * (3*q(i,jhi,klo-1,n) - q(i,jhi-1,klo-1,n))
                   end if

                end if

                !
                ! ****************************************************************************
                !

                if (j == jhi+1 .and. bc(2,2,n) .eq. HOEXTRAP .and. &
                    k == khi+1 .and. bc(3,2,n) .eq. HOEXTRAP) then

                   if (khi-2 >= ks) then
                      q(i,j,k,n) = half * eighth * (15*q(i,jhi+1,khi,n) - 10*q(i,jhi+1,khi-1,n) + 3*q(i,jhi+1,khi-2,n))
                   else
                      q(i,j,k,n) = half * half * (3*q(i,jhi+1,khi,n) - q(i,jhi+1,khi-1,n))
                   end if

                   if (jhi-2 >= js) then
                      q(i,j,k,n) = q(i,jhi+1,khi+1,n) + &
                                   half * eighth * (15*q(i,jhi,khi+1,n) - 10*q(i,jhi-1,khi+1,n) + 3*q(i,jhi-2,khi+1,n))
                   else
                      q(i,j,k,n) = q(i,jhi+1,khi+1,n) + half * half * (3*q(i,jhi,khi+1,n) - q(i,jhi-1,khi+1,n))
                   end if

                end if

             end do
          end do
       end do
    end do

  end subroutine filccn

end module filcc_module
