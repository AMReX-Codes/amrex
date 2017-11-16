
module init_prob_module
  use amrex_base_module
  implicit none
  private
  public :: init_prob_poisson, init_prob_abeclaplacian

contains

  subroutine init_prob_poisson (geom, solution, rhs, exact_solution)
    type(amrex_geometry), intent(in   ) :: geom(0:)
    type(amrex_multifab), intent(inout) :: solution(0:)
    type(amrex_multifab), intent(inout) :: rhs(0:)
    type(amrex_multifab), intent(inout) :: exact_solution(0:)

    integer :: ilev
    integer :: rlo(4), rhi(4), elo(4), ehi(4)
    type(amrex_box) :: bx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: prhs, pexact

    do ilev = 0, size(rhs)-1
       !$omp parallel private(rlo,rhi,elo,ehi,bx,mfi,prhs,pexact)
       call amrex_mfiter_build(mfi, rhs(ilev), tiling=.true.)
       
       do while (mfi%next())
          bx = mfi%tilebox()
          prhs   =>            rhs(ilev) % dataptr(mfi)
          pexact => exact_solution(ilev) % dataptr(mfi)
          rlo = lbound(prhs)
          rhi = ubound(prhs)
          elo = lbound(pexact)
          ehi = ubound(pexact)
          call actual_init_poisson(bx%lo, bx%hi, prhs, rlo(1:3), rhi(1:3), pexact, elo(1:3), ehi(1:3), &
               amrex_problo, amrex_probhi, )
       end do

       call amrex_mfiter_destroy(mfi)
       !$omp end parallel
    end do

  end subroutine init_prob_poisson


  subroutine init_prob_abeclaplacian (geom, solution, rhs, exact_solution, acoef, bcoef)
    type(amrex_geometry), intent(in   ) :: geom(0:)
    type(amrex_multifab), intent(inout) :: solution(0:)
    type(amrex_multifab), intent(inout) :: rhs(0:)
    type(amrex_multifab), intent(inout) :: exact_solution(0:)
    type(amrex_multifab), intent(inout) :: acoef(0:)
    type(amrex_multifab), intent(inout) :: bcoef(0:)
  end subroutine init_prob_abeclaplacian


  subroutine actual_init_poisson (lo, hi, rhs, rlo, rhi, exact, elo, ehi)
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, elo, ehi
    real(amrex_real), intent(inout) :: rhs  (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(inout) :: exact(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))
  end subroutine actual_init_poisson

end module init_prob_module
