subroutine integrate_state_force(lo, hi, &
                                 state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                 diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                 dx, time, a, half_dt)
!
!   Calculates the sources to be added later on.
!
!   Parameters
!   ----------
!   lo : double array (3)
!       The low corner of the current box.
!   hi : double array (3)
!       The high corner of the current box.
!   state_* : double arrays
!       The state vars
!   diag_eos_* : double arrays
!       Temp and Ne
!   src_* : doubles arrays
!       The source terms to be added to state (iterative approx.)
!   double array (3)
!       The low corner of the entire domain
!   a : double
!       The current a
!   half_dt : double
!       time step size, in Mpc km^-1 s ~ 10^12 yr.
!
!   Returns
!   -------
!   state : double array (dims) @todo
!       The state vars
!
    use amrex_fort_module, only : rt => amrex_real

    use forcing_spect_module
    use eos_module, only: eos_assume_neutral
    use network, only: nspec, aion, zion
    use atomic_rates_module, only: XHYDROGEN
    use probdata_module, only: prob_lo, prob_hi, alpha, rho0, temp0
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   NDIAG, TEMP_COMP, NE_COMP, small_pres, small_temp, gamma_minus_1
    use bl_constants_module, only : TWO, ONE, HALF, ZERO, M_PI, M_SQRT_2
    use fundamental_constants_module
 
    implicit none

    ! get the mass of a nucleon from Avogadro's number.
    real(rt), parameter :: m_nucleon = 1.d0/n_A

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in)    :: dx(3), time, a, half_dt

    integer :: i, j, k
    integer :: m, mi, mj, mk, n
    integer :: alloc, neg_e_count
    real(rt) :: m_nucleon_over_kB, mu, sum_y

    real(rt) :: rho, rho_e_orig, rho_K_res, T_orig, ne
    real(rt) :: delta, delta_re, eint0, press, small_eint

    integer :: num_phases(3)
    real(rt) :: xn_eos(nspec), ymass(nspec) 
    real(rt) :: delta_phase(3), phase_lo(3)
    real(rt) :: accel(3), buf(num_modes_ext) 
    real(rt) :: phasefct_init_even(num_modes_ext), phasefct_init_odd(num_modes_ext)
    real(rt) :: phasefct_mult_even(num_modes_ext,3), phasefct_mult_odd(num_modes_ext,3)
    real(rt) :: phasefct_yz(num_modes_ext,2)
    real(rt), allocatable :: phasefct_even_x(:), phasefct_even_y(:), phasefct_even_z(:)
    real(rt), allocatable :: phasefct_odd_x(:), phasefct_odd_y(:), phasefct_odd_z(:)

    ! Compute mu and small_eint to avoid EOS calls that prevent loop vectorization
    m_nucleon_over_kB = m_nucleon / k_B

    xn_eos(1) = XHYDROGEN
    xn_eos(2) = (1.d0 - XHYDROGEN)

    sum_y  = 0.d0

    if (eos_assume_neutral) then
       ! assume completely neutral atoms
       do n = 1, nspec
          ymass(n) = xn_eos(n)/aion(n)
          sum_y = sum_y + ymass(n)
       enddo
    else
       ! assume completely ionized species
       do n = 1, nspec
          ymass(n) = xn_eos(n)*(1.d0 + zion(n))/aion(n)
          sum_y = sum_y + ymass(n)
       enddo
    endif

    mu = 1.d0/sum_y

    small_eint = small_temp / (mu * m_nucleon_over_kB * gamma_minus_1)

    neg_e_count = 0

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT
    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                ! Original values
                rho        = state(i,j,k,URHO)
                rho_e_orig = state(i,j,k,UEINT)
		rho_K_res  = state(i,j,k,UEDEN) - state(i,j,k,UEINT)
                T_orig     = diag_eos(i,j,k,TEMP_COMP)
                !ne         = diag_eos(i,j,k,  NE_COMP)

                if (rho_e_orig .lt. 0.d0) neg_e_count = neg_e_count + 1

                ! Compute temperature increment and ensure that new temperature is positive
		delta = half_dt * alpha * (temp0 - T_orig) / a
		diag_eos(i,j,k,TEMP_COMP) = max(T_orig + delta, small_temp)

                ! Call EOS to get internal energy for constant equilibrium temperature
                eint0 = temp0 / (mu * m_nucleon_over_kB * gamma_minus_1)
		delta_re = half_dt * alpha * (rho*eint0 - rho_e_orig) / a

                ! Call EOS to get the internal energy floor

                ! Update cell quantities
 		state(i,j,k,UEINT) = max(rho_e_orig + delta_re, rho*small_eint)
                state(i,j,k,UEDEN) = state(i,j,k,UEINT) + rho_K_res
            end do
        end do
    end do

    if (neg_e_count > 0) call bl_abort('bad rho e in integrate_state_force_3d')

    delta_phase(:) = TWO*M_PI * dx(:) / (prob_hi(:) - prob_lo(:)) ! phase increment per cell
    phase_lo(:) = (dble(lo(:)) + HALF) * delta_phase(:)           ! phase of low corner

    ! compute initial phase factors and multiplying factors
    ! (sin and cos are expensive, so we do that only for low corner and cell width)
    do m = 1, num_modes_ext
       i = wavevectors(1,m)
       j = wavevectors(2,m)
       k = wavevectors(3,m)

       phasefct_init_even(m) = &
          (cos(i*phase_lo(1)) * cos(j*phase_lo(2)) - &
           sin(i*phase_lo(1)) * sin(j*phase_lo(2))) * cos(k*phase_lo(3)) - &
          (cos(i*phase_lo(1)) * sin(j*phase_lo(2)) + &
           sin(i*phase_lo(1)) * cos(j*phase_lo(2))) * sin(k*phase_lo(3))

       phasefct_init_odd(m) = &
          (cos(i*phase_lo(1)) * cos(j*phase_lo(2)) - &
           sin(i*phase_lo(1)) * sin(j*phase_lo(2))) * sin(k*phase_lo(3)) + &
          (cos(i*phase_lo(1)) * sin(j*phase_lo(2)) + &
           sin(i*phase_lo(1)) * cos(j*phase_lo(2))) * cos(k*phase_lo(3))

       phasefct_mult_even(m,1) = cos(i*delta_phase(1));
       phasefct_mult_odd (m,1) = sin(i*delta_phase(1));

       phasefct_mult_even(m,2) = cos(j*delta_phase(2));
       phasefct_mult_odd (m,2) = sin(j*delta_phase(2));

       phasefct_mult_even(m,3) = cos(k*delta_phase(3));
       phasefct_mult_odd (m,3) = sin(k*delta_phase(3));
    end do

    num_phases(:) = (hi(:)-lo(:)+1)*num_modes_ext

    allocate(phasefct_even_x(num_phases(1)), phasefct_even_y(num_phases(2)), phasefct_even_z(num_phases(3)), &
             phasefct_odd_x(num_phases(1)),  phasefct_odd_y(num_phases(2)),  phasefct_odd_z(num_phases(3)), &
             STAT=alloc)

    if (alloc > 0) call bl_abort('failed to allocate arrays for phase factors')      
 
    ! initialize phase factors for each coordinate axis: 
    ! since phase factors for inverse FT are given by 
    ! exp(i*(k1*x + k2*y + k3*z)) = exp(i*k1*x) * exp(i*k2*y)*...,
    ! we iteratively multiply with exp(i*k1*delta_x), etc.
    do m = 1, num_modes_ext
       phasefct_even_x(m) = ONE 
       phasefct_odd_x(m)  = ZERO
    end do
    do i = lo(1)+1,hi(1)
       mi = (i-lo(1))*num_modes_ext + 1
       do m = 1, num_modes_ext
            buf(m) = phasefct_even_x(mi-num_modes_ext);
            phasefct_even_x(mi) = phasefct_mult_even(m,1) * phasefct_even_x(mi-num_modes_ext) - &
                                  phasefct_mult_odd (m,1) * phasefct_odd_x(mi-num_modes_ext)
            phasefct_odd_x(mi)  = phasefct_mult_even(m,1) * phasefct_odd_x(mi-num_modes_ext) + &
                                  phasefct_mult_odd (m,1) * buf(m)
            mi = mi + 1
       end do
    end do         

    do m = 1, num_modes_ext
       phasefct_even_y(m) = ONE 
       phasefct_odd_y(m)  = ZERO 
    end do
    do j = lo(2)+1,hi(2)
       mj = (j-lo(2))*num_modes_ext + 1
       do m = 1, num_modes_ext
            buf(m) = phasefct_even_y(mj-num_modes_ext);
            phasefct_even_y(mj) = phasefct_mult_even(m,2) * phasefct_even_y(mj-num_modes_ext) - &
                                  phasefct_mult_odd (m,2) * phasefct_odd_y(mj-num_modes_ext)
            phasefct_odd_y(mj)  = phasefct_mult_even(m,2) * phasefct_odd_y(mj-num_modes_ext) + &
                                  phasefct_mult_odd (m,2) * buf(m)
            mj = mj + 1
       end do
    end do         

    do m = 1, num_modes_ext
       phasefct_even_z(m) = phasefct_init_even(m)  
       phasefct_odd_z(m)  = phasefct_init_odd(m)
    end do
    do k = lo(3)+1, hi(3)
       mk = (k-lo(3))*num_modes_ext + 1
       do m = 1, num_modes_ext
            buf(m) = phasefct_even_z(mk-num_modes_ext);
            phasefct_even_z(mk) = phasefct_mult_even(m,3) * phasefct_even_z(mk-num_modes_ext) - &
                                  phasefct_mult_odd (m,3) * phasefct_odd_z(mk-num_modes_ext)
            phasefct_odd_z(mk)  = phasefct_mult_even(m,3) * phasefct_odd_z(mk-num_modes_ext) + &
                                  phasefct_mult_odd (m,3) * buf(m)
            mk = mk + 1
       end do
    end do

    ! apply forcing in physical space
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          mj = (j-lo(2))*num_modes_ext + 1 ! offset in y-direction
          mk = (k-lo(3))*num_modes_ext + 1 ! offset in z-direction

          ! pre-compute products of phase factors depending on y- and z-coordinates 
          do m = 1, num_modes_ext
             phasefct_yz(m,1) = phasefct_even_y(mj) * phasefct_even_z(mk) - phasefct_odd_y(mj)  * phasefct_odd_z(mk)
             phasefct_yz(m,2) = phasefct_odd_y(mj)  * phasefct_even_z(mk) + phasefct_even_y(mj) * phasefct_odd_z(mk)
             mj = mj + 1
             mk = mk + 1
          end do

          do i = lo(1),hi(1)

             accel(:) = (/ ZERO, ZERO, ZERO /)

             ! compute components of acceleration via inverse FT 
             do n = 1, 3
                mi = (i-lo(1))*num_modes_ext + 1 ! offset in x-direction
  
                !dir$ vector
                do m = 1, num_modes_ext
                   ! sum up even modes
                   accel(n) = accel(n) + (phasefct_even_x(mi) * phasefct_yz(m,1) - &
                                          phasefct_odd_x(mi)  * phasefct_yz(m,2)) * modes_even(m,n)
                   ! sum up odd modes
                   accel(n) = accel(n) - (phasefct_even_x(mi) * phasefct_yz(m,2) + &
                                          phasefct_odd_x(mi)  * phasefct_yz(m,1)) * modes_odd(m,n)
                   mi = mi + 1
                 end do
             end do

             accel(:) = M_SQRT_2 * accel(:)

             ! add forcing to state		
             state(i,j,k,UMX) = state(i,j,k,UMX) + half_dt * state(i,j,k,URHO)*accel(1) / a
             state(i,j,k,UMY) = state(i,j,k,UMY) + half_dt * state(i,j,k,URHO)*accel(2) / a
             state(i,j,k,UMZ) = state(i,j,k,UMZ) + half_dt * state(i,j,k,URHO)*accel(3) / a
          end do
       end do
    end do

    ! update total energy
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5d0*(state(i,j,k,UMX)*state(i,j,k,UMX) + &
                                                              state(i,j,k,UMY)*state(i,j,k,UMY) + &
                                                              state(i,j,k,UMZ)*state(i,j,k,UMZ))/state(i,j,k,URHO)
          end do
       end do
    end do

end subroutine integrate_state_force

