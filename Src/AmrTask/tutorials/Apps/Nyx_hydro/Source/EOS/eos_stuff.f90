!
! We compute the mean molecular weight, mu, based on the mass fractions
! of the different species.
!
!   NOTE: in the helmholtz EOS, we use Abar, which is the weighted
!   ion mass.  Abar does not include any free electron contribution.
!
! There are 2 ways to compute mu, depending on whether the gas is
! completely ionized or not.
!
! For a neutral gas, the mean molecular weight is:
!
!   1/mu = sum_k { X_k / A_k }
!
! For a completely ionized gas, the mean molecular weight is:
!
!   1/mu = sum_k { (1 + Z_k) X_k / A_k }
!
! At the moment, we will select between these 2 ionization extremes
! (completely neutral vs. completely ionized) by a hard-coded
! parameter, eos_assume_neutral
!
! The ratio of specific heats (gamma) is allowed to vary.  NOTE: the
! expression for entropy is only valid for an ideal MONATOMIC gas
! (gamma = 5/3).  

module eos_module

  use bl_constants_module, only: M_PI, ONE
  use network, only: nspec, aion, zion
  use atomic_rates_module, only: XHYDROGEN

  implicit none

  private

  integer, parameter, private :: eos_input_rt = 1   ! density, temperature are inputs
  integer, parameter, private :: eos_input_re = 5   ! density, internal energy are inputs

  logical, parameter, public :: eos_assume_neutral = .false.

  private nspec, aion, zion

  public eos_init_small_pres, nyx_eos_T_given_Re, nyx_eos_T_given_Re_vec, nyx_eos_S_given_Re, &
         nyx_eos_soundspeed, nyx_eos_given_RT, nyx_eos_given_RT_vec, eos

contains

  !---------------------------------------------------------------------------
  ! Nyx interfaces 
  !---------------------------------------------------------------------------

  subroutine eos_init_small_pres(R, T, Ne, P, comoving_a)

     use amrex_fort_module, only : rt => amrex_real

     ! In/out variables
     real(rt), intent(  out) :: P
     real(rt), intent(in   ) :: R, T, Ne
     real(rt), intent(in   ) :: comoving_a

     ! Local variables
     logical :: do_diag

     real(rt) :: xn_eos(nspec)
     real(rt) :: temp_eos
     real(rt) :: den_eos
     real(rt) :: e_eos
     real(rt) :: p_eos
     real(rt) :: cv_eos
     real(rt) :: dpdt_eos
     real(rt) :: dpdr_eos
     real(rt) :: dedt_eos
     real(rt) ::    s_eos
     real(rt) :: comoving_a_cubed

     do_diag = .false.

     comoving_a_cubed = (comoving_a*comoving_a*comoving_a)

     ! Density is the only variable we convert from comoving to proper coordinates
     den_eos = R / comoving_a_cubed

     temp_eos = T

     xn_eos(1) = XHYDROGEN
     xn_eos(2) = (1.d0 - XHYDROGEN)

     call eos(eos_input_rt, den_eos, temp_eos, &
              xn_eos, &
              p_eos, e_eos, cv_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, &
              s_eos, &
              do_diag)

    ! Pressure must be converted from proper to comoving coordinates
    P  = p_eos * comoving_a_cubed

  end subroutine eos_init_small_pres

  subroutine nyx_eos_soundspeed(c, R, e)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module, only: gamma_const, gamma_minus_1

     ! In/out variables
     real(rt), intent(in   ) :: R, e
     real(rt), intent(  out) :: c

     ! Pressure
     real(rt) :: P

     P = R * e * gamma_minus_1

     ! sound speed
     c = sqrt(gamma_const * P / R)

  end subroutine nyx_eos_soundspeed

  subroutine nyx_eos_T_given_Re(JH, JHe, T, Ne, R, e, comoving_a)

     use amrex_fort_module, only : rt => amrex_real

     ! In/out variables
     integer, intent(in) :: JH, JHe ! stubs here
     real(rt),           intent(inout) :: T, Ne
     real(rt),           intent(in   ) :: R, e
     real(rt),           intent(in   ) :: comoving_a

     ! Local variables
     logical :: do_diag

     real(rt) :: xn_eos(nspec)
     real(rt) :: temp_eos
     real(rt) :: den_eos
     real(rt) :: e_eos
     real(rt) :: p_eos
     real(rt) :: cv_eos
     real(rt) :: dpdt_eos
     real(rt) :: dpdr_eos
     real(rt) :: dedt_eos
     real(rt) ::    s_eos
     real(rt) :: comoving_a_cubed

     do_diag = .false.

     comoving_a_cubed = comoving_a**3

     ! Density is the only variable we convert from comoving to proper coordinates
     den_eos = R / comoving_a_cubed

     temp_eos = T 
     e_eos = e 

     xn_eos(1) = XHYDROGEN
     xn_eos(2) = (1.d0 - XHYDROGEN)

     call eos(eos_input_re, den_eos, temp_eos, &
              xn_eos, &
              p_eos, e_eos, cv_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, &
              s_eos, &
              do_diag)

    T  = temp_eos

    Ne = 1.d0

  end subroutine nyx_eos_T_given_Re

  subroutine nyx_eos_S_given_Re(S, R, e, T, Ne, comoving_a)

     use amrex_fort_module, only : rt => amrex_real
     implicit none

     ! In/out variables
     real(rt),           intent(  out) :: S
     real(rt),           intent(in   ) :: R, e, T, Ne
     real(rt),           intent(in   ) :: comoving_a

     ! Local variables
     logical :: do_diag

     real(rt) :: xn_eos(nspec)
     real(rt) :: temp_eos
     real(rt) :: den_eos
     real(rt) :: e_eos
     real(rt) :: p_eos
     real(rt) :: cv_eos
     real(rt) :: dpdt_eos
     real(rt) :: dpdr_eos
     real(rt) :: dedt_eos
     real(rt) ::    s_eos

     do_diag = .false.

     ! Density is the only variable we convert from comoving to proper coordinates
     den_eos = R / (comoving_a*comoving_a*comoving_a)

     temp_eos = T
     e_eos = e

     xn_eos(1) = XHYDROGEN
     xn_eos(2) = (1.d0 - XHYDROGEN)

     call eos(eos_input_re, den_eos, temp_eos, &
              xn_eos, &
              p_eos, e_eos, cv_eos,  &
              dpdt_eos, dpdr_eos, dedt_eos, &
              s_eos, &
              do_diag)

    ! Should this be weighted by comoving a??
    S  = s_eos 

  end subroutine nyx_eos_S_given_Re

  subroutine nyx_eos_given_RT(e, P, R, T, Ne, comoving_a)

     use amrex_fort_module, only : rt => amrex_real

     ! In/out variables
     real(rt),           intent(  out) :: e, P
     real(rt),           intent(in   ) :: R, T, Ne
     real(rt),           intent(in   ) :: comoving_a


     ! Local variables
     logical :: do_diag
     
     real(rt) :: xn_eos(nspec)
     real(rt) :: temp_eos
     real(rt) :: den_eos
     real(rt) :: e_eos
     real(rt) :: p_eos
     real(rt) :: cv_eos
     real(rt) :: dpdt_eos
     real(rt) :: dpdr_eos
     real(rt) :: dedt_eos
     real(rt) ::    s_eos
     real(rt) :: comoving_a_cubed

     do_diag = .false.

     comoving_a_cubed = (comoving_a*comoving_a*comoving_a)

     ! Density is the only variable we convert from comoving to proper coordinates
     den_eos = R / comoving_a_cubed

     temp_eos = T 

     xn_eos(1) = XHYDROGEN
     xn_eos(2) = (1.d0 - XHYDROGEN)

     call eos(eos_input_rt, den_eos, temp_eos, &
              xn_eos, &
              p_eos, e_eos, cv_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, &
              s_eos, &
              do_diag)

    ! Pressure must be converted from proper to comoving coordinates
    P  = p_eos * comoving_a_cubed

    e  = e_eos

  end subroutine nyx_eos_given_RT

  !---------------------------------------------------------------------------
  ! The main interface -- this is used directly by MAESTRO
  !---------------------------------------------------------------------------
  subroutine eos(input, dens, temp, &
                 xmass, &
                 pres, eint, c_v, &
                 dPdT, dPdR, dEdT, &
                 entropy, &
                 do_eos_diag)

     use amrex_fort_module, only : rt => amrex_real
     use bl_error_module
     use fundamental_constants_module, only: k_B, n_A, hbar
     use meth_params_module, only: gamma_minus_1

! dens     -- mass density (g/cc)
! temp     -- temperature (K)
! xmass    -- the mass fractions of the individual isotopes
! pres     -- the pressure (dyn/cm**2)
! eint     -- the internal energy (erg/g)
! c_v      -- specific heat at constant volume
! ne       -- number density of electrons + positrons
! dPdT     -- d pressure/ d temperature
! dPdR     -- d pressure/ d density
! dEdT     -- d energy/ d temperature
! entropy  -- entropy (erg/g/K)  NOTE: presently the entropy expression is 
!             valid only for an ideal MONATOMIC gas (gamma = 5/3).
!
! input = 1 means temp, dens    , and xmass are inputs, return eint    , etc
!       = 5 means dens, eint    , and xmass are inputs, return temp    , etc

    implicit none

    logical do_eos_diag
    integer, intent(in) :: input

    real(rt) :: dens, temp
    real(rt) :: xmass(nspec)
    real(rt) :: pres, eint
    real(rt) :: c_v
    real(rt) :: dPdT, dPdR, dedT
    real(rt) :: entropy

    ! local variables
    real(rt) :: ymass(nspec)    
    real(rt) :: mu
    real(rt) :: sum_y
    real(rt) :: m_nucleon_over_kB
    real(rt) :: t1,t2,t3

    ! get the mass of a nucleon from Avogadro's number.
    real(rt), parameter :: m_nucleon = 1.d0/n_A

    integer :: n

    !-------------------------------------------------------------------------
    ! compute mu -- the mean molecular weight
    !-------------------------------------------------------------------------

    sum_y  = 0.d0

    if (eos_assume_neutral) then
       ! assume completely neutral atoms
       do n = 1, nspec
          ymass(n) = xmass(n)/aion(n)
          sum_y = sum_y + ymass(n)
       enddo
    else
       ! assume completely ionized species
       do n = 1, nspec
          ymass(n) = xmass(n)*(1.d0 + zion(n))/aion(n)
          sum_y = sum_y + ymass(n)
       enddo
    endif

    mu = 1.d0/sum_y

    !-------------------------------------------------------------------------
    ! for all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs
    !-------------------------------------------------------------------------

    ! These are both very large numbers so we pre-divide
    m_nucleon_over_kB = m_nucleon / k_B

    if (input .eq. eos_input_rt) then

       ! Compute e = k T / [(mu m_nucleon)*(gamma-1)] from T
       eint = temp / (mu * m_nucleon_over_kB * gamma_minus_1)

    else ! We are given eint instead of temp

       ! Solve e = k T / [(mu m_nucleon)*(gamma-1)] for T
       temp = eint * (mu * m_nucleon_over_kB * gamma_minus_1)

    endif

    !-------------------------------------------------------------------------
    ! now we have the density and temperature (and mass fractions / mu)
    ! regardless of the inputs.
    !-------------------------------------------------------------------------

    ! compute the pressure simply from the ideal gas law
    pres = gamma_minus_1 * dens * eint

    ! entropy (per gram) of an ideal monoatomic gas (the Sactur-Tetrode equation)
    ! NOTE: this expression is only valid for gamma = 5/3.

    t1 = (mu*m_nucleon);     t1 = t1*t1*sqrt(t1)
    t2 = (k_B*temp);         t2 = t2*sqrt(t2)
    t3 = (2*M_PI*hbar*hbar); t3 = t3*sqrt(t3)

    entropy = (1.d0 / (mu*m_nucleon_over_kB)) * (2.5d0 + log(t1/dens*t2/t3))

    ! compute the thermodynamic derivatives and specific heats
    dPdT = pres/temp
    dPdR = pres/dens
    dedT = eint/temp

    c_v = dedT

  end subroutine eos


  subroutine nyx_eos_T_given_Re_vec(T, Ne, R_in, e_in, a, veclen)

    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only: amrex_abort
  
    ! In/out variables
    integer, intent(in) :: veclen
    real(rt), dimension(veclen), intent(inout) :: T, Ne
    real(rt), dimension(veclen), intent(in   ) :: R_in, e_in
    real(rt),                    intent(in   ) :: a
  
    call amrex_abort("nyx_eos_T_given_Re_vec supported only with USE_HEATCOOL=TRUE and USE_CVODE=TRUE")

  end subroutine nyx_eos_T_given_Re_vec


  subroutine nyx_eos_given_RT_vec(e, P, R, T, Ne, a, veclen)

    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only: amrex_abort
    implicit none
  
    integer, intent(in) :: veclen
    real(rt), dimension(veclen), intent(  out) :: e, P
    real(rt), dimension(veclen), intent(in   ) :: R, T, Ne
    real(rt),          intent(in   ) :: a
  
    call amrex_abort("nyx_eos_given_RT_vec supported only with USE_HEATCOOL=TRUE and USE_CVODE=TRUE")

  end subroutine nyx_eos_given_RT_vec

end module eos_module
