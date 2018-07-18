module rhs_helmholtz
  use amrex_base_module
  use amrex_amr_module
  use amr_data_module

  implicit none
  private
  public :: init_prob_abeclaplacian

  real(amrex_real), private, parameter :: a = 1.0d0
  real(amrex_real), private, parameter :: b = 2.0d0

contains

  subroutine init_prob_abeclaplacian

    integer :: ilev
    integer :: rlo(4), rhi(4), elo(4), ehi(4), alo(4), ahi(4), blo(4), bhi(4)
    type(amrex_box) :: bx, gbx
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: prhs, pexact, pa, pb
    integer :: nlevs
    integer :: nc,ng,idim
    logical :: nodal(3)

    ascalar = a
    bscalar = b

    nlevs = amrex_get_finest_level()

    do ilev = 0, nlevs

      call amrex_multifab_build(solution(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc=1, ng=1)
      call amrex_multifab_build(exact_solution(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc=1, ng=0)
      call amrex_multifab_build(rhs(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc=1, ng=0)
      call amrex_multifab_build(acoef(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc=1, ng=0)
      call amrex_multifab_build(bcoef(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, nc=1, ng=1)

      do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(beta(idim,ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, 1, 0, nodal)
       end do
   enddo

    do ilev = 0, nlevs

       !$omp parallel private(rlo,rhi,elo,ehi,alo,ahi,blo,bhi,bx,mfi,prhs,pexact,pa,pb)
       call amrex_mfiter_build(mfi, rhs(ilev), tiling=.true.)
       
       do while (mfi%next())
          bx = mfi%tilebox()
          gbx = mfi%growntilebox(1)
          prhs   =>            rhs(ilev) % dataptr(mfi)
          pexact => exact_solution(ilev) % dataptr(mfi)
          pa     =>          acoef(ilev) % dataptr(mfi)
          pb     =>          bcoef(ilev) % dataptr(mfi)
          rlo = lbound(prhs)
          rhi = ubound(prhs)
          elo = lbound(pexact)
          ehi = ubound(pexact)
          alo = lbound(pa)
          ahi = ubound(pa)
          blo = lbound(pb)
          bhi = ubound(pb)
          call actual_init_abeclapcian(bx%lo, bx%hi, gbx%lo, gbx%hi, &
               prhs, rlo(1:3), rhi(1:3), pexact, elo(1:3), ehi(1:3), &
               pa, alo(1:3), ahi(1:3), pb, blo(1:3), bhi(1:3), &
               amrex_problo, amrex_probhi, amrex_geom(ilev)%dx)

       end do

       call amrex_mfiter_destroy(mfi)
       !$omp end parallel

       ! This will be used to provide bc and initial guess for the solver.
       call solution(ilev)%setVal(0.0_amrex_real)
    end do

  end subroutine init_prob_abeclaplacian

  subroutine actual_init_abeclapcian (lo, hi, glo, ghi, rhs, rlo, rhi, exact, elo, ehi, &
       alpha, alo, ahi, beta, blo, bhi, prob_lo, prob_hi, dx)
    integer, dimension(3), intent(in) :: lo, hi, glo, ghi, rlo, rhi, elo, ehi, alo, ahi, blo, bhi
    real(amrex_real), intent(inout) :: rhs  (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(inout) :: exact(elo(1):ehi(1),elo(2):ehi(2),elo(3):ehi(3))
    real(amrex_real), intent(inout) :: alpha(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    real(amrex_real), intent(inout) :: beta (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))
    real(amrex_real), dimension(3), intent(in) :: prob_lo, prob_hi, dx

    integer :: i,j,k
    real(amrex_real) x, y, z, xc, yc, zc
    real(amrex_real) r, theta, dbdrfac
    real(amrex_real) pi, fpi, tpi, fac
    real(amrex_real), parameter :: w = 0.05d0
    real(amrex_real), parameter :: sigma = 10.d0

    pi = 4.d0 * atan(1.d0)
    tpi = 2.0d0 * pi
    fpi = 4.0d0 * pi
#if BL_SPACEDIM == 2
    fac = 8.d0 * pi**2
#endif
#if BL_SPACEDIM == 3
    fac = 12.d0 * pi**2
#endif

    xc = (prob_hi(1) + prob_lo(1))/2.d0
    yc = (prob_hi(2) + prob_lo(2))/2.d0
    zc = (prob_hi(3) + prob_lo(3))/2.d0

    theta = 0.5d0*log(3.d0) / (w + 1.d-50)
    
    do k = glo(3), ghi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = glo(2), ghi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = glo(1), ghi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
    
#if BL_SPACEDIM == 2
             r = sqrt((x-xc)**2 + (y-yc)**2)
#endif
#if BL_SPACEDIM == 3	
             r = sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
#endif	
             beta(i,j,k) = (sigma-1.d0)/2.d0*tanh(theta*(r-0.25d0)) + (sigma+1.d0)/2.d0
          end do
       end do
    end do
    
    do k = lo(3), hi(3)
       z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
       do j = lo(2), hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i = lo(1), hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
             
#if BL_SPACEDIM == 2
             r = sqrt((x-xc)**2 + (y-yc)**2)
#endif
#if BL_SPACEDIM == 3	
             r = sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
#endif	
            
             dbdrfac = (sigma-1.d0)/2.d0/(cosh(theta*(r-0.25d0)))**2 * theta/r
             dbdrfac = dbdrfac * b
             
             alpha(i,j,k) = 1.d0

#if BL_SPACEDIM == 2
             exact(i,j,k) = 1.d0 * cos(tpi*x) * cos(tpi*y)   &
                  &      + .25d0 * cos(fpi*x) * cos(fpi*y)


             rhs(i,j,k) = beta(i,j,k)*b*fac*(cos(tpi*x) * cos(tpi*y)    &
                  &                        + cos(fpi*x) * cos(fpi*y))   &
                  &   + dbdrfac*((x-xc)*(tpi*sin(tpi*x) * cos(tpi*y)    &
                  &                     + pi*sin(fpi*x) * cos(fpi*y))   &
                  &            + (y-yc)*(tpi*cos(tpi*x) * sin(tpi*y)    &
                  &                     + pi*cos(fpi*x) * sin(fpi*y)))  &
                  &                   + a * (cos(tpi*x) * cos(tpi*y)    &
                  &               + 0.25d0 * cos(fpi*x) * cos(fpi*y))
#endif



#if BL_SPACEDIM == 3
             exact(i,j,k) = 1.d0 * cos(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
                  &      + .25d0 * cos(fpi*x) * cos(fpi*y) * cos(fpi*z)


             rhs(i,j,k) = beta(i,j,k)*b*fac*(cos(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
                  &                        + cos(fpi*x) * cos(fpi*y) * cos(fpi*z))  &
                  &   + dbdrfac*((x-xc)*(tpi*sin(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
                  &                     + pi*sin(fpi*x) * cos(fpi*y) * cos(fpi*z))  &
                  &            + (y-yc)*(tpi*cos(tpi*x) * sin(tpi*y) * cos(tpi*z)   &
                  &                     + pi*cos(fpi*x) * sin(fpi*y) * cos(fpi*z))  &
                  &            + (z-zc)*(tpi*cos(tpi*x) * cos(tpi*y) * sin(tpi*z)   &
                  &                     + pi*cos(fpi*x) * cos(fpi*y) * sin(fpi*z))) &
                  &                   + a * (cos(tpi*x) * cos(tpi*y) * cos(tpi*z)   &
                  &               + 0.25d0 * cos(fpi*x) * cos(fpi*y) * cos(fpi*z))
#endif

          end do
       end do
    end do

  end subroutine actual_init_abeclapcian

end module rhs_helmholtz
