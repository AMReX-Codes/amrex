
module amrex_filcc_module

  use amrex_fort_module, only : amrex_real, amrex_spacedim, amrex_get_loop_bounds
  use amrex_bc_types_module
  use amrex_constants_module

  implicit none

  interface amrex_filcc
     module procedure amrex_filcc_1
     module procedure amrex_filcc_n
  end interface amrex_filcc

  private
  public :: amrex_filcc, amrex_fab_filcc, amrex_filccn, amrex_hoextraptocc
#if (AMREX_SPACEDIM == 3)
  public :: amrex_hoextraptocc_3d
#endif
#if (AMREX_SPACEDIM == 2)
  public :: amrex_hoextraptocc_2d
#endif
#if defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA)
  public :: amrex_filccn_device
#endif

#ifndef AMREX_XSDK
  public :: filccn
#endif

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
       call amrex_filccn(qlo(1:3), qhi(1:3), q(:,:,:,i), qlo(1:3), qhi(1:3), 1, &
            domlo, domhi, dx, xlo, bc);
    end do
  end subroutine amrex_filcc_n

#if (AMREX_SPACEDIM == 3)

  subroutine amrex_filcc_1(q,qlo1,qlo2,qlo3,qhi1,qhi2,qhi3,domlo,domhi,dx,xlo,bc)
    integer, intent(in) :: qlo1,qlo2,qlo3,qhi1,qhi2,qhi3,domlo(amrex_spacedim),domhi(amrex_spacedim)
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim, 2)
    real(amrex_real), intent(inout) :: q(qlo1:qhi1,qlo2:qhi2,qlo3:qhi3)
    integer :: q_lo(3), q_hi(3)
    q_lo = [qlo1,qlo2,qlo3]
    q_hi = [qhi1,qhi2,qhi3]
    call amrex_filccn(q_lo, q_hi, q, q_lo, q_hi, 1, domlo, domhi, dx, xlo, bc);
  end subroutine amrex_filcc_1

#elif (AMREX_SPACEDIM == 2)

  subroutine amrex_filcc_1(q,qlo1,qlo2,qhi1,qhi2,domlo,domhi,dx,xlo,bc)
    integer, intent(in) :: qlo1,qlo2,qhi1,qhi2,domlo(amrex_spacedim),domhi(amrex_spacedim)
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim, 2)
    real(amrex_real), intent(inout) :: q(qlo1:qhi1,qlo2:qhi2)
    integer :: q_lo(3), q_hi(3)
    q_lo = [qlo1,qlo2,0]
    q_hi = [qhi1,qhi2,0]
    call amrex_filccn(q_lo, q_hi, q, q_lo, q_hi, 1, domlo, domhi, dx, xlo, bc);
  end subroutine amrex_filcc_1

#else

  subroutine amrex_filcc_1(q,qlo1,qhi1,domlo,domhi,dx,xlo,bc)
    integer, intent(in) :: qlo1,qhi1,domlo(amrex_spacedim),domhi(amrex_spacedim)
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim, 2)
    real(amrex_real), intent(inout) :: q(qlo1:qhi1)
    integer :: q_lo(3), q_hi(3)
    q_lo = [qlo1,0,0]
    q_hi = [qhi1,0,0]
    call amrex_filccn(q_lo, q_hi, q, q_lo, q_hi, 1, domlo, domhi, dx, xlo, bc);
  end subroutine amrex_filcc_1

#endif

  subroutine amrex_fab_filcc (q, qlo, qhi, nq, domlo, domhi, dx, xlo, bc) &
       bind(c, name='amrex_fab_filcc')

    implicit none

    integer, intent(in) :: qlo(3), qhi(3), nq
    integer, dimension(amrex_spacedim), intent(in) :: domlo, domhi
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim,2,nq)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nq)

    integer :: lo(3), hi(3)

    call amrex_get_loop_bounds(lo, hi, qlo, qhi)

    call amrex_filccn(lo, hi, q, qlo, qhi, nq, domlo, domhi, dx, xlo, bc)

  end subroutine amrex_fab_filcc

#ifndef AMREX_XSDK
  subroutine filccn(lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)
    implicit none
    integer,          intent(in   ) :: lo(3), hi(3)
    integer,          intent(in   ) :: q_lo(3), q_hi(3)
    integer,          intent(in   ) :: ncomp
    integer,          intent(in   ) :: domlo(amrex_spacedim), domhi(amrex_spacedim)
    real(amrex_real), intent(in   ) :: xlo(amrex_spacedim), dx(amrex_spacedim)
    real(amrex_real), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),ncomp)
    integer,          intent(in   ) :: bc(amrex_spacedim,2,ncomp)
    call amrex_filccn(lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)
  end subroutine filccn
#endif

  subroutine amrex_filccn(lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)

    implicit none

    integer,          intent(in   ) :: lo(3), hi(3)
    integer,          intent(in   ) :: q_lo(3), q_hi(3)
    integer,          intent(in   ) :: ncomp
    integer,          intent(in   ) :: domlo(amrex_spacedim), domhi(amrex_spacedim)
    real(amrex_real), intent(in   ) :: xlo(amrex_spacedim), dx(amrex_spacedim)
    real(amrex_real), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),ncomp)
    integer,          intent(in   ) :: bc(amrex_spacedim,2,ncomp)

    integer :: ilo, ihi, jlo, jhi, klo, khi
    integer :: is, ie, js, je, ks, ke
    integer :: i, j, k, n
    integer :: imin, imax, jmin, jmax, kmin, kmax

    is = max(q_lo(1), domlo(1))
    ie = min(q_hi(1), domhi(1))
    ilo = domlo(1)
    ihi = domhi(1)

#if AMREX_SPACEDIM >= 2
    js = max(q_lo(2), domlo(2))
    je = min(q_hi(2), domhi(2))
    jlo = domlo(2)
    jhi = domhi(2)
#endif

#if AMREX_SPACEDIM == 3
    ks = max(q_lo(3), domlo(3))
    ke = min(q_hi(3), domhi(3))
    klo = domlo(3)
    khi = domhi(3)
#endif

    do n = 1, ncomp

       if (lo(1) < ilo) then
          imin = lo(1)
          imax = ilo-1

          if (bc(1,1,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(1,1,n) .eq. amrex_bc_foextrap) then
             
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ilo,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,1,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax

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
             
          else if (bc(1,1,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ilo+(ilo-i)-1,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,1,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = -q(ilo+(ilo-i)-1,j,k,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(1) > ihi) then
          imin = ihi+1
          imax = hi(1)

          if (bc(1,2,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(1,2,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ihi,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,2,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax

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

          else if (bc(1,2,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ihi-(i-ihi)+1,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,2,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = -q(ihi-(i-ihi)+1,j,k,n)
                   end do
                end do
             end do
             
          end if

       end if

#if AMREX_SPACEDIM >= 2

       if (lo(2) < jlo) then
          jmin = lo(2)
          jmax = jlo-1

          if (bc(2,1,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(2,1,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jlo,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,1,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)

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

          else if (bc(2,1,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jlo+(jlo-j)-1,k,n)
                   end do
                end do
             end do

          else if (bc(2,1,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,jlo+(jlo-j)-1,k,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(2) > jhi) then
          jmin = jhi+1
          jmax = hi(2)

          if (bc(2,2,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(2,2,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jhi,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,2,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)

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

          else if (bc(2,2,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jhi-(j-jhi)+1,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,2,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,jhi-(j-jhi)+1,k,n)
                   end do
                end do
             end do
             
          end if

       end if
#endif




#if AMREX_SPACEDIM == 3

       if (lo(3) < klo) then
          kmin = lo(3)
          kmax = klo-1

          if (bc(3,1,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.
             
          else if (bc(3,1,n) .eq. amrex_bc_foextrap) then
             
             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,klo,n)
                   end do
                end do
             end do

          else if (bc(3,1,n) .eq. amrex_bc_hoextrap) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)

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
             
          else if (bc(3,1,n) .eq. amrex_bc_reflect_even) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,klo+(klo-k)-1,n)
                   end do
                end do
             end do
             
          else if (bc(3,1,n) .eq. amrex_bc_reflect_odd) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,j,klo+(klo-k)-1,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(3) > khi) then
          kmin = khi+1
          kmax = hi(3)

          if (bc(3,2,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.
             
          else if (bc(3,2,n) .eq. amrex_bc_foextrap) then
             
             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,khi,n)
                   end do
                end do
             end do

          else if (bc(3,2,n) .eq. amrex_bc_hoextrap) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)

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
             
          else if (bc(3,2,n) .eq. amrex_bc_reflect_even) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,khi-(k-khi)+1,n)
                   end do
                end do
             end do
             
          else if (bc(3,2,n) .eq. amrex_bc_reflect_odd) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,j,khi-(k-khi)+1,n)
                   end do
                end do
             end do

          end if

       end if
#endif



#if AMREX_SPACEDIM >= 2
       ! Now take care of the higher contributions

       !
       ! First correct the i-j edges and all corners
       !

       if (bc(1,1,n) .eq. amrex_bc_hoextrap .and. bc(2,1,n) .eq. amrex_bc_hoextrap) then

          if (lo(1) < ilo .and. lo(2) < jlo) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)

             i = ilo-1
             j = jlo-1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)
                      
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
                   
#if AMREX_SPACEDIM == 3
                   
                   if (k == klo-1 .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jlo-1,klo,n) - 10*q(ilo-1,jlo-1,klo+1,n) + &
                              3*q(ilo-1,jlo-1,klo+2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jlo-1,klo,n) - q(ilo-1,jlo-1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jlo-1,khi,n) - 10*q(ilo-1,jlo-1,khi-1,n) + &
                              3*q(ilo-1,jlo-1,khi-2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jlo-1,khi,n) - q(ilo-1,jlo-1,khi-1,n))
                      end if
                   end if
#endif
                   
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,1,n) .eq. amrex_bc_hoextrap .and. bc(2,2,n) .eq. amrex_bc_hoextrap) then

          if (lo(1) < ilo .and. hi(2) > jhi) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)

             i = ilo-1
             j = jhi+1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)

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

#if (AMREX_SPACEDIM == 3)
                   if (k == klo-1 .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jhi+1,klo,n) - 10*q(ilo-1,jhi+1,klo+1,n) + &
                              3*q(ilo-1,jhi+1,klo+2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jhi+1,klo,n) - q(ilo-1,jhi+1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jhi+1,khi,n) - 10*q(ilo-1,jhi+1,khi-1,n) + &
                              3*q(ilo-1,jhi+1,khi-2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jhi+1,khi,n) - q(ilo-1,jhi+1,khi-1,n))
                      end if
                   end if

#endif

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. amrex_bc_hoextrap .and. bc(2,1,n) .eq. amrex_bc_hoextrap) then

          if (hi(1) > ihi .and. lo(2) < jlo) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)

             i = ihi+1
             j = jlo-1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)

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
                   
#if (AMREX_SPACEDIM == 3)
                   if (k == klo-1 .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jlo-1,klo,n) - 10*q(ihi+1,jlo-1,klo+1,n) + 3*q(ihi+1,jlo-1,klo+2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jlo-1,klo,n) - q(ihi+1,jlo-1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jlo-1,khi,n) - 10*q(ihi+1,jlo-1,khi-1,n) + 3*q(ihi+1,jlo-1,khi-2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jlo-1,khi,n) - q(ihi+1,jlo-1,khi-1,n))
                      end if
                   end if
#endif

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. amrex_bc_hoextrap .and. bc(2,2,n) .eq. amrex_bc_hoextrap) then

          if (hi(1) > ihi .and. hi(2) > jhi) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)

             i = ihi+1
             j = jhi+1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)
                   
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
                   
#if (AMREX_SPACEDIM == 3)
                   if (k == klo-1 .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jhi+1,klo,n) - 10*q(ihi+1,jhi+1,klo+1,n) + 3*q(ihi+1,jhi+1,klo+2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jhi+1,klo,n) - q(ihi+1,jhi+1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jhi+1,khi,n) - 10*q(ihi+1,jhi+1,khi-1,n) + 3*q(ihi+1,jhi+1,khi-2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jhi+1,khi,n) - q(ihi+1,jhi+1,khi-1,n))
                      end if
                   end if
#endif
                   
                end do
             end if
          end if
       end if
#endif

#if AMREX_SPACEDIM == 3
       !
       ! Next correct the i-k edges
       !

       if (bc(1,1,n) .eq. amrex_bc_hoextrap .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then

          if (lo(1) < ilo .and. lo(3) < klo) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             i = ilo-1
             k = klo-1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then
             
                do j = lo(2), hi(2)

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

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,1,n) .eq. amrex_bc_hoextrap .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then

          if (lo(1) < ilo .and. hi(3) > khi) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             i = ilo-1
             k = khi+1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then

                do j = lo(2), hi(2)

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
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. amrex_bc_hoextrap .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then

          if (hi(1) > ihi .and. lo(3) < klo) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             i = ihi+1
             k = klo-1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then

                do j = lo(2), hi(2)

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
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. amrex_bc_hoextrap .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then
          
          if (hi(1) > ihi .and. hi(3) > khi) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             i = ihi+1
             k = khi+1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then

                do j = lo(2), hi(2)

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
                   
                end do
             end if
          end if
       end if

       !
       ! Next correct the j-k edges
       !

       if (bc(2,1,n) .eq. amrex_bc_hoextrap .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then

          if (lo(2) < jlo .and. lo(3) < klo) then
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             j = jlo-1
             k = klo-1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

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

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,1,n) .eq. amrex_bc_hoextrap .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then

          if (lo(2) < jlo .and. hi(3) > khi) then
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             j = jlo-1
             k = khi+1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

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
                   
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,2,n) .eq. amrex_bc_hoextrap .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then

          if (hi(2) > jhi .and. lo(3) < klo) then
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             j = jhi+1
             k = klo-1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

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
                   
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,2,n) .eq. amrex_bc_hoextrap .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then

          if (hi(2) > jhi .and. hi(3) > khi) then
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             j = jhi+1
             k = khi+1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

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
                   
                end do
             end if
          end if
       end if
#endif

    end do

  end subroutine amrex_filccn

#if defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) subroutine amrex_filccn_device(lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)

    implicit none

    integer,          intent(in   ) :: lo(3), hi(3)
    integer,          intent(in   ) :: q_lo(3), q_hi(3)
    integer,          intent(in   ) :: ncomp
    integer,          intent(in   ) :: domlo(amrex_spacedim), domhi(amrex_spacedim)
    real(amrex_real), intent(in   ) :: xlo(amrex_spacedim), dx(amrex_spacedim)
    real(amrex_real), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),ncomp)
    integer,          intent(in   ) :: bc(amrex_spacedim,2,ncomp)

    integer :: ilo, ihi, jlo, jhi, klo, khi
    integer :: is, ie, js, je, ks, ke
    integer :: i, j, k, n
    integer :: imin, imax, jmin, jmax, kmin, kmax

    is = max(q_lo(1), domlo(1))
    ie = min(q_hi(1), domhi(1))
    ilo = domlo(1)
    ihi = domhi(1)

#if AMREX_SPACEDIM >= 2
    js = max(q_lo(2), domlo(2))
    je = min(q_hi(2), domhi(2))
    jlo = domlo(2)
    jhi = domhi(2)
#endif

#if AMREX_SPACEDIM == 3
    ks = max(q_lo(3), domlo(3))
    ke = min(q_hi(3), domhi(3))
    klo = domlo(3)
    khi = domhi(3)
#endif

    do n = 1, ncomp

       if (lo(1) < ilo) then
          imin = lo(1)
          imax = min(hi(1),ilo-1)

          if (bc(1,1,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(1,1,n) .eq. amrex_bc_foextrap) then
             
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ilo,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,1,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax

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
             
          else if (bc(1,1,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ilo+(ilo-i)-1,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,1,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = -q(ilo+(ilo-i)-1,j,k,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(1) > ihi) then
          imin = max(lo(1),ihi+1)
          imax = hi(1)

          if (bc(1,2,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(1,2,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ihi,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,2,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax

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

          else if (bc(1,2,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ihi-(i-ihi)+1,j,k,n)
                   end do
                end do
             end do
             
          else if (bc(1,2,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = -q(ihi-(i-ihi)+1,j,k,n)
                   end do
                end do
             end do
             
          end if

       end if

       ! For CUDA we need to synchronize the threadblock after each
       ! dimension, since the results for the corners depend
       ! on the i, j, and k directions being done in order.
       ! Note: this will only work if the threadblock size is
       ! larger in each dimension than the number of ghost zones.

       call syncthreads()

#if AMREX_SPACEDIM >= 2

       if (lo(2) < jlo) then
          jmin = lo(2)
          jmax = min(hi(2),jlo-1)

          if (bc(2,1,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(2,1,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jlo,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,1,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)

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

          else if (bc(2,1,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jlo+(jlo-j)-1,k,n)
                   end do
                end do
             end do

          else if (bc(2,1,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,jlo+(jlo-j)-1,k,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(2) > jhi) then
          jmin = max(lo(2),jhi+1)
          jmax = hi(2)

          if (bc(2,2,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(2,2,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jhi,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,2,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)

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

          else if (bc(2,2,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,jhi-(j-jhi)+1,k,n)
                   end do
                end do
             end do
             
          else if (bc(2,2,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = jmin, jmax
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,jhi-(j-jhi)+1,k,n)
                   end do
                end do
             end do
             
          end if

       end if
#endif

       call syncthreads()



#if AMREX_SPACEDIM == 3

       if (lo(3) < klo) then
          kmin = lo(3)
          kmax = min(hi(3),klo-1)

          if (bc(3,1,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.
             
          else if (bc(3,1,n) .eq. amrex_bc_foextrap) then
             
             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,klo,n)
                   end do
                end do
             end do

          else if (bc(3,1,n) .eq. amrex_bc_hoextrap) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)

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
             
          else if (bc(3,1,n) .eq. amrex_bc_reflect_even) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,klo+(klo-k)-1,n)
                   end do
                end do
             end do
             
          else if (bc(3,1,n) .eq. amrex_bc_reflect_odd) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,j,klo+(klo-k)-1,n)
                   end do
                end do
             end do
             
          end if

       end if

       if (hi(3) > khi) then
          kmin = max(lo(3),khi+1)
          kmax = hi(3)

          if (bc(3,2,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.
             
          else if (bc(3,2,n) .eq. amrex_bc_foextrap) then
             
             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,khi,n)
                   end do
                end do
             end do

          else if (bc(3,2,n) .eq. amrex_bc_hoextrap) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)

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
             
          else if (bc(3,2,n) .eq. amrex_bc_reflect_even) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = q(i,j,khi-(k-khi)+1,n)
                   end do
                end do
             end do
             
          else if (bc(3,2,n) .eq. amrex_bc_reflect_odd) then

             do k = kmin, kmax
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      q(i,j,k,n) = -q(i,j,khi-(k-khi)+1,n)
                   end do
                end do
             end do

          end if

       end if
#endif

       call syncthreads()



#if AMREX_SPACEDIM >= 2
       ! Now take care of the higher contributions

       !
       ! First correct the i-j edges and all corners
       !

       if (bc(1,1,n) .eq. amrex_bc_hoextrap .and. bc(2,1,n) .eq. amrex_bc_hoextrap) then

          if (lo(1) < ilo .and. lo(2) < jlo) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)

             i = ilo-1
             j = jlo-1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)
                      
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
                   
#if AMREX_SPACEDIM == 3
                   
                   if (k == klo-1 .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jlo-1,klo,n) - 10*q(ilo-1,jlo-1,klo+1,n) + &
                              3*q(ilo-1,jlo-1,klo+2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jlo-1,klo,n) - q(ilo-1,jlo-1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jlo-1,khi,n) - 10*q(ilo-1,jlo-1,khi-1,n) + &
                              3*q(ilo-1,jlo-1,khi-2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jlo-1,khi,n) - q(ilo-1,jlo-1,khi-1,n))
                      end if
                   end if
#endif
                   
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,1,n) .eq. amrex_bc_hoextrap .and. bc(2,2,n) .eq. amrex_bc_hoextrap) then

          if (lo(1) < ilo .and. hi(2) > jhi) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)

             i = ilo-1
             j = jhi+1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)

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

#if (AMREX_SPACEDIM == 3)
                   if (k == klo-1 .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jhi+1,klo,n) - 10*q(ilo-1,jhi+1,klo+1,n) + &
                              3*q(ilo-1,jhi+1,klo+2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jhi+1,klo,n) - q(ilo-1,jhi+1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * ( (15*q(ilo-1,jhi+1,khi,n) - 10*q(ilo-1,jhi+1,khi-1,n) + &
                              3*q(ilo-1,jhi+1,khi-2,n)) )
                      else
                         q(i,j,k,n) = half * (3*q(ilo-1,jhi+1,khi,n) - q(ilo-1,jhi+1,khi-1,n))
                      end if
                   end if

#endif

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. amrex_bc_hoextrap .and. bc(2,1,n) .eq. amrex_bc_hoextrap) then

          if (hi(1) > ihi .and. lo(2) < jlo) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)

             i = ihi+1
             j = jlo-1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)

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
                   
#if (AMREX_SPACEDIM == 3)
                   if (k == klo-1 .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jlo-1,klo,n) - 10*q(ihi+1,jlo-1,klo+1,n) + 3*q(ihi+1,jlo-1,klo+2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jlo-1,klo,n) - q(ihi+1,jlo-1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jlo-1,khi,n) - 10*q(ihi+1,jlo-1,khi-1,n) + 3*q(ihi+1,jlo-1,khi-2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jlo-1,khi,n) - q(ihi+1,jlo-1,khi-1,n))
                      end if
                   end if
#endif

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. amrex_bc_hoextrap .and. bc(2,2,n) .eq. amrex_bc_hoextrap) then

          if (hi(1) > ihi .and. hi(2) > jhi) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)

             i = ihi+1
             j = jhi+1

             if (i.ge.imin .and. i.le.imax .and. j.ge.jmin .and. j.le.jmax) then

                do k = lo(3), hi(3)
                   
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
                   
#if (AMREX_SPACEDIM == 3)
                   if (k == klo-1 .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then
                      if (klo+2 <= ke) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jhi+1,klo,n) - 10*q(ihi+1,jhi+1,klo+1,n) + 3*q(ihi+1,jhi+1,klo+2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jhi+1,klo,n) - q(ihi+1,jhi+1,klo+1,n))
                      end if
                   end if
                   
                   if (k == khi+1 .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then
                      if (khi-2 >= ks) then
                         q(i,j,k,n) = eighth * (15*q(ihi+1,jhi+1,khi,n) - 10*q(ihi+1,jhi+1,khi-1,n) + 3*q(ihi+1,jhi+1,khi-2,n))
                      else
                         q(i,j,k,n) = half * (3*q(ihi+1,jhi+1,khi,n) - q(ihi+1,jhi+1,khi-1,n))
                      end if
                   end if
#endif
                   
                end do
             end if
          end if
       end if
#endif

#if AMREX_SPACEDIM == 3
       !
       ! Next correct the i-k edges
       !

       if (bc(1,1,n) .eq. amrex_bc_hoextrap .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then

          if (lo(1) < ilo .and. lo(3) < klo) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             i = ilo-1
             k = klo-1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then
             
                do j = lo(2), hi(2)

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

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,1,n) .eq. amrex_bc_hoextrap .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then

          if (lo(1) < ilo .and. hi(3) > khi) then
             imin = lo(1)
             imax = min(hi(1),ilo-1)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             i = ilo-1
             k = khi+1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then

                do j = lo(2), hi(2)

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
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. amrex_bc_hoextrap .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then

          if (hi(1) > ihi .and. lo(3) < klo) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             i = ihi+1
             k = klo-1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then

                do j = lo(2), hi(2)

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
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(1,2,n) .eq. amrex_bc_hoextrap .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then
          
          if (hi(1) > ihi .and. hi(3) > khi) then
             imin = max(lo(1),ihi+1)
             imax = hi(1)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             i = ihi+1
             k = khi+1

             if (i.ge.imin .and. i.le.imax .and. k.ge.kmin .and. k.le.kmax) then

                do j = lo(2), hi(2)

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
                   
                end do
             end if
          end if
       end if

       !
       ! Next correct the j-k edges
       !

       if (bc(2,1,n) .eq. amrex_bc_hoextrap .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then

          if (lo(2) < jlo .and. lo(3) < klo) then
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             j = jlo-1
             k = klo-1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

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

                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,1,n) .eq. amrex_bc_hoextrap .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then

          if (lo(2) < jlo .and. hi(3) > khi) then
             jmin = lo(2)
             jmax = min(hi(2),jlo-1)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             j = jlo-1
             k = khi+1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

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
                   
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,2,n) .eq. amrex_bc_hoextrap .and. bc(3,1,n) .eq. amrex_bc_hoextrap) then

          if (hi(2) > jhi .and. lo(3) < klo) then
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)
             kmin = lo(3)
             kmax = min(hi(3),klo-1)

             j = jhi+1
             k = klo-1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

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
                   
                end do
             end if
          end if
       end if

       !
       ! ****************************************************************************
       !

       if (bc(2,2,n) .eq. amrex_bc_hoextrap .and. bc(3,2,n) .eq. amrex_bc_hoextrap) then

          if (hi(2) > jhi .and. hi(3) > khi) then
             jmin = max(lo(2),jhi+1)
             jmax = hi(2)
             kmin = max(lo(3),khi+1)
             kmax = hi(3)

             j = jhi+1
             k = khi+1

             if (j.ge.jmin .and. j.le.jmax .and. k.ge.kmin .and. k.le.kmax) then

                do i = lo(1), hi(1)

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
                   
                end do
             end if
          end if
       end if
#endif

    end do

  end subroutine amrex_filccn_device
#endif

  subroutine amrex_hoextraptocc (q, qlo, qhi, domlo, domhi, dx, xlo) &
       bind(c,name='amrex_hoextraptocc')
    integer, intent(in) :: qlo(3), qhi(3), domlo(*), domhi(*)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    real(amrex_real), intent(in) :: dx(*), xlo(*)

#if (AMREX_SPACEDIM == 3)
    call amrex_hoextraptocc_3d(q,qlo(1),qlo(2),qlo(3),qhi(1),qhi(2),qhi(3),domlo,domhi,dx,xlo)
#elif (AMREX_SPACEDIM == 2)
    call amrex_hoextraptocc_2d(q,qlo(1),qlo(2),qhi(1),qhi(2),domlo,domhi,dx,xlo)
#endif
  end subroutine amrex_hoextraptocc


#if (AMREX_SPACEDIM == 3)
subroutine amrex_hoextraptocc_3d(q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,domlo,domhi,dx,xlo)

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module

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

end subroutine amrex_hoextraptocc_3d
#endif

#if (AMREX_SPACEDIM == 2)
subroutine amrex_hoextraptocc_2d(q,q_l1,q_l2,q_h1,q_h2,domlo,domhi,dx,xlo)

  use amrex_fort_module
  use amrex_constants_module

  implicit none

  integer    q_l1, q_l2, q_h1, q_h2
  integer    domlo(2), domhi(2)
  real(amrex_real)     xlo(2), dx(2)
  real(amrex_real)     q(q_l1:q_h1,q_l2:q_h2)

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

end subroutine amrex_hoextraptocc_2d
#endif

end module amrex_filcc_module
