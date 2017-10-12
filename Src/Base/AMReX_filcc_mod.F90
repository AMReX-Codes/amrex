
module amrex_filcc_module

  use amrex_fort_module, only : amrex_real, amrex_spacedim, get_loop_bounds
  use filcc_module, only: filccn
#ifdef AMREX_USE_CUDA
  use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

  implicit none

  interface amrex_filcc
     module procedure amrex_filcc_1
     module procedure amrex_filcc_n
  end interface amrex_filcc

  private
  public :: amrex_filcc, amrex_fab_filcc

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

  AMREX_LAUNCH_SUBROUTINE subroutine amrex_fab_filcc (q, qlo, qhi, nq, domlo, domhi, dx, xlo, bc) &
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

end module amrex_filcc_module
