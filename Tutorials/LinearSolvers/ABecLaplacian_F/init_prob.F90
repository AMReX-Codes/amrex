
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
               amrex_problo, amrex_probhi, geom(ilev)%dx)
       end do

       call amrex_mfiter_destroy(mfi)
       !$omp end parallel

       call solution(ilev)%setVal(0.0_amrex_real)  ! This will be used to provide bc.
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


  subroutine actual_init_poisson (lo, hi, rhs, rlo, rhi, exact, elo, ehi, prob_lo, prob_hi, dx)
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, elo, ehi
    real(amrex_real), intent(inout) :: rhs  (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(inout) :: exact(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))
    real(amrex_real), dimension(3), intent(in) :: prob_lo, prob_hi, dx

    integer :: i,j,k
    real(amrex_real) :: x, y, z
    real(amrex_real), parameter :: tpi =  8.d0*atan(1.0)
    real(amrex_real), parameter :: fpi = 16.d0*atan(1.0)
    real(amrex_real), parameter :: fac = tpi*tpi*3.d0

    do k = lo(3), hi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = lo(2), hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = lo(1), hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
             
             exact(i,j,k) = 1.d0 * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z))  &
                  &      + .25d0 * (sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
                
             rhs(i,j,k) = -fac * (sin(tpi*x) * sin(tpi*y) * sin(tpi*z))  &
                  &       -fac * (sin(fpi*x) * sin(fpi*y) * sin(fpi*z))
          end do
       end do
    end do

  end subroutine actual_init_poisson

end module init_prob_module
