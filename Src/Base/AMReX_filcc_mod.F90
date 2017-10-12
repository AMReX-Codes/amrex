
#include "AMReX_BC_TYPES.H"
#include "AMReX_CONSTANTS.H"

module amrex_filcc_module

  use amrex_fort_module, only : amrex_real, amrex_spacedim, get_loop_bounds
#ifdef AMREX_USE_CUDA
  use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

  implicit none

  interface amrex_filcc
     module procedure amrex_filcc_1
     module procedure amrex_filcc_n
  end interface amrex_filcc

  private
  public :: amrex_filcc, amrex_fab_filcc, filccn

contains

  subroutine amrex_filcc_n(q,qlo,qhi,domlo,domhi,dx,xlo,bclo,bchi)
    integer, intent(in) :: qlo(4), qhi(4)
    integer, dimension(amrex_spacedim), intent(in) :: domlo, domhi
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bclo(amrex_spacedim,*), bchi(amrex_spacedim,*)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qlo(4):qhi(4))
    integer :: i, bc(amrex_spacedim,2)
    do i = qlo(4), qhi(4)
       bc(:,1) = bclo(:,i)
       bc(:,2) = bchi(:,i)
#if (BL_SPACEDIM == 3)
       call filcc(q(:,:,:,i),qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3),domlo,domhi,dx,xlo,bc)
#elif (BL_SPACEDIM == 2)
       call filcc(q(:,:,:,i),qlo(1),qlo(2),       qhi(1),qhi(2),       domlo,domhi,dx,xlo,bc)
#else
       call filcc(q(:,:,:,i),qlo(1),              qhi(1),              domlo,domhi,dx,xlo,bc)
#endif
    end do
  end subroutine amrex_filcc_n

#if (BL_SPACEDIM == 3)

  subroutine amrex_filcc_1(q,qlo1,qlo2,qlo3,qhi1,qhi2,qhi3,domlo,domhi,dx,xlo,bc)
    integer, intent(in) :: qlo1,qlo2,qlo3,qhi1,qhi2,qhi3,domlo(amrex_spacedim),domhi(amrex_spacedim)
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim, 2)
    real(amrex_real), intent(inout) :: q(qlo1:qhi1,qlo2:qhi2,qlo3:qhi3)
    call filcc(q,qlo1,qlo2,qlo3,qhi1,qhi2,qhi3,domlo,domhi,dx,xlo,bc)
  end subroutine amrex_filcc_1

#elif (BL_SPACEDIM == 2)

  subroutine amrex_filcc_1(q,qlo1,qlo2,qhi1,qhi2,domlo,domhi,dx,xlo,bc)
    integer, intent(in) :: qlo1,qlo2,qhi1,qhi2,domlo(amrex_spacedim),domhi(amrex_spacedim)
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim, 2)
    real(amrex_real), intent(inout) :: q(qlo1:qhi1,qlo2:qhi2)
    call filcc(q,qlo1,qlo2,qhi1,qhi2,domlo,domhi,dx,xlo,bc)
  end subroutine amrex_filcc_1

  subroutine amrex_filcc_2(q,qlo1,qlo2,qhi1,qhi2,domlo,domhi,dx,xlo,bclo,bchi)
    integer, intent(in) :: qlo1,qlo2,qhi1,qhi2,domlo(amrex_spacedim),domhi(amrex_spacedim)
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bclo(amrex_spacedim), bchi(amrex_spacedim)
    real(amrex_real), intent(inout) :: q(qlo1:qhi1,qlo2:qhi2)
    integer :: bc(amrex_spacedim,2)
    bc(:,1) = bclo
    bc(:,2) = bchi
    call filcc(q,qlo1,qlo2,qhi1,qhi2,domlo,domhi,dx,xlo,bc)
  end subroutine amrex_filcc_2

#else

  subroutine amrex_filcc_1(q,qlo1,qhi1,domlo,domhi,dx,xlo,bc)
    integer, intent(in) :: qlo1,qhi1,domlo(amrex_spacedim),domhi(amrex_spacedim)
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim, 2)
    real(amrex_real), intent(inout) :: q(qlo1:qhi1)
    call filcc(q,qlo1,qhi1,domlo,domhi,dx,xlo,bc)
  end subroutine amrex_filcc_1

#endif

  AMREX_LAUNCH subroutine amrex_fab_filcc (q, qlo, qhi, nq, domlo, domhi, dx, xlo, bc) &
       bind(c, name='amrex_fab_filcc')

    implicit none

    integer, intent(in) :: qlo(3), qhi(3), nq
    integer, dimension(amrex_spacedim), intent(in) :: domlo, domhi
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim,2,nq)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nq)

    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, qlo, qhi)

    call filccn(blo, bhi, q, qlo, qhi, nq, domlo, domhi, dx, xlo, bc)

  end subroutine amrex_fab_filcc



  AMREX_DEVICE subroutine filccn(blo, bhi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)

    implicit none

    integer,          intent(in   ) :: blo(3), bhi(3)
    integer,          intent(in   ) :: q_lo(3), q_hi(3)
    integer,          intent(in   ) :: domlo(amrex_spacedim), domhi(amrex_spacedim)
    real(amrex_real), intent(in   ) :: xlo(amrex_spacedim), dx(amrex_spacedim)
    real(amrex_real), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),ncomp)
    integer,          intent(in   ) :: bc(amrex_spacedim,2,ncomp)
    integer,          intent(in   ) :: ncomp

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

       if (bc(1,1,n) .eq. EXT_DIR) then

          ! Do nothing.

       else if (bc(1,1,n) .eq. FOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i < ilo) then
                      q(i,j,k,n) = q(ilo,j,k,n)
                   end if

                end do
             end do
          end do

       else if (bc(1,1,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i < ilo - 1) then
                      q(i,j,k,n) = q(ilo,j,k,n)
                   else if (i == ilo - 1) then
                      if (ilo+2 <= ie) then
                         q(i,j,k,n) = eighth * (15*q(ilo,j,k,n) - 10*q(ilo+1,j,k,n) + 3*q(ilo+2,j,k,n))
                      else
                         q(i,j,k,n) = half * (3*q(ilo,j,k,n) - q(ilo+1,j,k,n))
                      end if
                   end if

                end do
             end do
          end do

       else if (bc(1,1,n) .eq. REFLECT_EVEN) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i < ilo) then
                      q(i,j,k,n) = q(ilo+(ilo-i)-1,j,k,n)
                   end if

                end do
             end do
          end do

       else if (bc(1,1,n) .eq. REFLECT_ODD) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i < ilo) then
                      q(i,j,k,n) = -q(ilo+(ilo-i)-1,j,k,n)
                   end if

                end do
             end do
          end do

       end if



       if (bc(1,2,n) .eq. EXT_DIR) then

          ! Do nothing.

       else if (bc(1,2,n) .eq. FOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i > ihi) then
                      q(i,j,k,n) = q(ihi,j,k,n)
                   end if

                end do
             end do
          end do

       else if (bc(1,2,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i > ihi + 1) then
                      q(i,j,k,n) = q(ihi,j,k,n)
                   else if (i == ihi + 1) then
                      if (ihi-2 >= is) then
                         q(i,j,k,n) = eighth * (15*q(ihi,j,k,n) - 10*q(ihi-1,j,k,n) + 3*q(ihi-2,j,k,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi,j,k,n) - q(ihi-1,j,k,n))
                      end if
                   end if

                end do
             end do
          end do

       else if (bc(1,2,n) .eq. REFLECT_EVEN) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i > ihi) then
                      q(i,j,k,n) = q(ihi-(i-ihi)+1,j,k,n)
                   end if

                end do
             end do
          end do

       else if (bc(1,2,n) .eq. REFLECT_ODD) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i > ihi) then
                      q(i,j,k,n) = -q(ihi-(i-ihi)+1,j,k,n)
                   end if

                end do
             end do
          end do

       end if



#if AMREX_SPACEDIM >= 2
       if (bc(2,1,n) .eq. EXT_DIR) then

          ! Do nothing.

       else if (bc(2,1,n) .eq. FOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j < jlo) then
                      q(i,j,k,n) = q(i,jlo,k,n)
                   end if

                end do
             end do
          end do

       else if (bc(2,1,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j < jlo - 1) then
                      q(i,j,k,n) = q(i,jlo,k,n)
                   else if (j == jlo - 1) then
                      if (jlo+2 <= je) then
                         q(i,j,k,n) = eighth * (15*q(i,jlo,k,n) - 10*q(i,jlo+1,k,n) + 3*q(i,jlo+2,k,n))
                      else
                         q(i,j,k,n) = half * (3*q(i,jlo,k,n) - q(i,jlo+1,k,n))
                      end if
                   end if

                end do
             end do
          end do

       else if (bc(2,1,n) .eq. REFLECT_EVEN) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j < jlo) then
                      q(i,j,k,n) = q(i,jlo+(jlo-j)-1,k,n)
                   end if

                end do
             end do
          end do

       else if (bc(2,1,n) .eq. REFLECT_ODD) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j < jlo) then
                      q(i,j,k,n) = -q(i,jlo+(jlo-j)-1,k,n)
                   end if

                end do
             end do
          end do

       end if



       if (bc(2,2,n) .eq. EXT_DIR) then

          ! Do nothing.

       else if (bc(2,2,n) .eq. FOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j > jhi) then
                      q(i,j,k,n) = q(i,jhi,k,n)
                   end if

                end do
             end do
          end do

       else if (bc(2,2,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j > jhi + 1) then
                      q(i,j,k,n) = q(i,jhi,k,n)
                   else if (j == jhi + 1) then
                      if (jhi-2 >= js) then
                         q(i,j,k,n) = eighth * (15*q(i,jhi,k,n) - 10*q(i,jhi-1,k,n) + 3*q(i,jhi-2,k,n))
                      else
                         q(i,j,k,n) = half * (3*q(i,jhi,k,n) - q(i,jhi-1,k,n))
                      end if
                   end if

                end do
             end do
          end do

       else if (bc(2,2,n) .eq. REFLECT_EVEN) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j > jhi) then
                      q(i,j,k,n) = q(i,jhi-(j-jhi)+1,k,n)
                   end if

                end do
             end do
          end do

       else if (bc(2,2,n) .eq. REFLECT_ODD) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j > jhi) then
                      q(i,j,k,n) = -q(i,jhi-(j-jhi)+1,k,n)
                   end if

                end do
             end do
          end do

       end if
#endif



#if AMREX_SPACEDIM == 3
       if (bc(3,1,n) .eq. EXT_DIR) then

          ! Do nothing.

       else if (bc(3,1,n) .eq. FOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (k < klo) then
                      q(i,j,k,n) = q(i,j,klo,n)
                   end if

                end do
             end do
          end do

       else if (bc(3,1,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (k < klo - 1) then
                      q(i,j,k,n) = q(i,j,klo,n)
                   else if (k == klo - 1) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * (15*q(i,j,klo,n) - 10*q(i,j,klo+1,n) + 3*q(i,j,klo+2,n))
                      else
                         q(i,j,k,n) = half * (3*q(i,j,klo,n) - q(i,j,klo+1,n))
                      end if
                   end if

                end do
             end do
          end do

       else if (bc(3,1,n) .eq. REFLECT_EVEN) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (k < klo) then
                      q(i,j,k,n) = q(i,j,klo+(klo-k)-1,n)
                   end if

                end do
             end do
          end do

       else if (bc(3,1,n) .eq. REFLECT_ODD) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (k < klo) then
                      q(i,j,k,n) = -q(i,j,klo+(klo-k)-1,n)
                   end if

                end do
             end do
          end do

       end if



       if (bc(3,2,n) .eq. EXT_DIR) then

          ! Do nothing.

       else if (bc(3,2,n) .eq. FOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (k > khi) then
                      q(i,j,k,n) = q(i,j,khi,n)
                   end if

                end do
             end do
          end do

       else if (bc(3,2,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (k > khi + 1) then
                      q(i,j,k,n) = q(i,j,khi,n)
                   else if (k == khi + 1) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * (15*q(i,j,khi,n) - 10*q(i,j,khi-1,n) + 3*q(i,j,khi-2,n))
                      else
                         q(i,j,k,n) = half * (3*q(i,j,khi,n) - q(i,j,khi-1,n))
                      end if
                   end if

                end do
             end do
          end do

       else if (bc(3,2,n) .eq. REFLECT_EVEN) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (k > khi) then
                      q(i,j,k,n) = q(i,j,khi-(k-khi)+1,n)
                   end if

                end do
             end do
          end do

       else if (bc(3,2,n) .eq. REFLECT_ODD) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (k > khi) then
                      q(i,j,k,n) = -q(i,j,khi-(k-khi)+1,n)
                   end if

                end do
             end do
          end do

       end if
#endif



#if AMREX_SPACEDIM >= 2
       ! Now take care of the higher contributions

       !
       ! First correct the i-j edges and all corners
       !

       if (bc(1,1,n) .eq. HOEXTRAP .and. bc(2,1,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i == ilo-1 .and. j == jlo-1) then

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

                end do
             end do
          end do

       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,1,n) .eq. HOEXTRAP .and. bc(2,2,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i == ilo-1 .and. j == jhi+1) then

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

                end do
             end do
          end do

       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. HOEXTRAP .and. bc(2,1,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i == ihi+1 .and. j == jlo-1) then

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

                end do
             end do
          end do

       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. HOEXTRAP .and. bc(2,2,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i == ihi+1 .and. j == jhi+1) then

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

                end do
             end do
          end do

       end if
#endif

#if AMREX_SPACEDIM == 3
       !
       ! Next correct the i-k edges
       !

       if (bc(1,1,n) .eq. HOEXTRAP .and. bc(3,1,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i == ilo-1 .and. k == klo-1) then

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

                end do
             end do
          end do

       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,1,n) .eq. HOEXTRAP .and. bc(3,2,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i == ilo-1 .and. k == khi+1) then

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

                end do
             end do
          end do

       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. HOEXTRAP .and. bc(3,1,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i == ihi+1 .and. k == klo-1) then

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

                end do
             end do
          end do

       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. HOEXTRAP .and. bc(3,2,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (i == ihi+1 .and. k == khi+1) then

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

                end do
             end do
          end do

       end if

       !
       ! Next correct the j-k edges
       !

       if (bc(2,1,n) .eq. HOEXTRAP .and. bc(3,1,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j == jlo-1 .and. k == klo-1) then

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

                end do
             end do
          end do

       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,1,n) .eq. HOEXTRAP .and. bc(3,2,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j == jlo-1 .and. k == khi+1) then

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

                end do
             end do
          end do

       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,2,n) .eq. HOEXTRAP .and. bc(3,1,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j == jhi+1 .and. k == klo-1) then

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

                end do
             end do
          end do

       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,2,n) .eq. HOEXTRAP .and. bc(3,2,n) .eq. HOEXTRAP) then

          do k = blo(3), bhi(3)
             do j = blo(2), bhi(2)
                do i = blo(1), bhi(1)

                   if (j == jhi+1 .and. k == khi+1) then

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

       end if
#endif

    end do

  end subroutine filccn

end module amrex_filcc_module
