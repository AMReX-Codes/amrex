
program main

    use constants_module, only : rt => type_real, M_PI
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   NDIAG, TEMP_COMP, NE_COMP, ZHI_COMP, &
                                   gamma_minus_1
    use eos_params_module
!    use network
    use eos_module, only: nyx_eos_T_given_Re, nyx_eos_given_RT, fort_setup_eos_params
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB
    use comoving_nd_module, only: fort_integrate_comoving_a
    use atomic_rates_module, only: YHELIUM, fort_tabulate_rates, fort_interp_to_this_z
    use vode_aux_module    , only: JH_vode, JHe_vode, z_vode, i_vode, j_vode, k_vode, z_vode, fn_vode, NR_vode
    use reion_aux_module   , only: zhi_flash, zheii_flash, flash_h, flash_he, &
                                   T_zhi, T_zheii, inhomogeneous_on
    use, intrinsic :: iso_c_binding

  implicit none

    real(rt) :: a, half_dt
    integer :: i, j, k
    real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z
    real(rt) :: T_orig, ne_orig, e_orig
    real(rt) :: T_out , ne_out , e_out, mu, mean_rhob, T_H, T_He
    real(rt) :: species(5)

    integer(c_int) :: ierr       ! error flag from C functions
    real(c_double) :: tstart     ! initial time
    real(c_double) :: atol, rtol
    type(c_ptr) :: sunvec_y      ! sundials vector
    type(c_ptr) :: CVmem         ! CVODE memory
    integer(c_long), parameter :: neq = 1
    integer :: fn_out = 0
!  call amrex_init()
 
  real(c_double), pointer :: yvec(:)
  
    real(c_double) :: vode_atol_scaled_in, vode_rtol_in, xacc_in
  CHARACTER(LEN=80) :: FMT, arg
  CHARACTER(LEN=6)  :: string
  integer :: STRANG_COMP
!  integer :: i_loop, j_loop

    DO i = 1, command_argument_count()
       CALL getarg(i, arg)
       WRITE (*,*) arg
    END DO

    print*,"Created data"
    
    print*,"Read table"
    call fort_tabulate_rates()
    vode_atol_scaled_in = 1e-4
    vode_rtol_in = 1e-4
    xacc_in = 1e-6
    gamma_minus_1 = 2.d0/3.d0
    call fort_setup_eos_params(xacc_in, vode_rtol_in, vode_atol_scaled_in)
    
    print*,"Finished reading table"

    allocate(yvec(neq))
    
    fn_vode = 0
    NR_vode = 0
    print*,"Read parameters"

 FMT="(A6,I1,/,ES21.15,/,ES21.15E2,/,ES21.15,/,ES21.15,/,ES21.15,/,ES21.15,/,ES21.15)"

    open(1,FILE=arg)
    read(1,FMT) string, STRANG_COMP, a, half_dt, rho, T_orig, ne_orig, e_orig
    close(1)

    yvec(1) = e_orig

    print*,"Finished reading parameters:"
    print(FMT), string,STRANG_COMP, a, half_dt, rho, T_orig, ne_orig, e_orig

    z = 1.d0/a - 1.d0
    call fort_integrate_comoving_a(a, a_end, half_dt)
    z_end = 1.0d0/a_end - 1.0d0
    call fort_interp_to_this_z(z)
    z_vode = z

    mean_rhob = comoving_OmB * 3.d0*(comoving_h*100.d0)**2 / (8.d0*M_PI*Gconst)

    ! Flash reionization?
    if ((flash_h .eqv. .true.) .and. (z .gt. zhi_flash)) then
       JH_vode = 0
    else
       JH_vode = 1
    endif
    if ((flash_he .eqv. .true.) .and. (z .gt. zheii_flash)) then
       JHe_vode = 0
    else
       JHe_vode = 1
    endif

    if (flash_h ) H_reion_z  = zhi_flash
    if (flash_he) He_reion_z = zheii_flash

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT

!    sunvec_y = N_VMake_Serial(NEQ, yvec)
!    if (.not. c_associated(sunvec_y)) then
!        call amrex_abort('integrate_state_fcvode: sunvec = NULL')
!    end if

!    CVmem = FCVodeCreate(CV_BDF, CV_NEWTON)
!    if (.not. c_associated(CVmem)) then
!        call amrex_abort('integrate_state_fcvode: CVmem = NULL')
!    end if

    tstart = 0.0
    ! CVodeMalloc allocates variables and initialize the solver. We can initialize the solver with junk because once we enter the
    ! (i,j,k) loop we will immediately call fcvreinit which reuses the same memory allocated from CVodeMalloc but sets up new
    ! initial conditions.
!    ierr = FCVodeInit(CVmem, c_funloc(RhsFn), tstart, sunvec_y)
!    if (ierr /= 0) then
!       call amrex_abort('integrate_state_fcvode: FCVodeInit() failed')
!    end if

    ! Set dummy tolerances. These will be overwritten as soon as we enter the loop and reinitialize the solver.
    rtol = 1.0d-5
    atol = 1.0d-10
!    ierr = FCVodeSStolerances(CVmem, rtol, atol)
!    if (ierr /= 0) then
!      call amrex_abort('integrate_state_fcvode: FCVodeSStolerances() failed')
!    end if

!    ierr = FCVDense(CVmem, neq)
!    if (ierr /= 0) then
!       call amrex_abort('integrate_state_fcvode: FCVDense() failed')
!    end if
    
    !-----------------cut out do ijk loops        
               ! Original values
                rho     = rho !state(i,j,k,URHO)
                e_orig  = e_orig !state(i,j,k,UEINT) / rho
                T_orig  = T_orig!diag_eos(i,j,k,TEMP_COMP)
                ne_orig = ne_orig!diag_eos(i,j,k,  NE_COMP)
                
                if (inhomogeneous_on) then
                   H_reion_z = 1*H_reion_z!diag_eos(i,j,k,ZHI_COMP)
                   if (z .gt. H_reion_z) then
                      JH_vode = 0
                   else
                      JH_vode = 1
                   endif
                endif

                if (e_orig .lt. 0.d0) then
                    !$OMP CRITICAL
                    print *,'negative e entering strang integration ',z, i,j,k, rho/mean_rhob, e_orig
!                    call bl_abort('bad e in strang')
                    !$OMP END CRITICAL
                end if

                i_vode = i
                j_vode = j
                k_vode = k

                call vode_wrapper(half_dt,rho,T_orig,ne_orig,e_orig, &
                                              T_out ,ne_out ,e_out, fn_out)
                print*, "vw:e_out   = ",e_out
                print*, "vw:T_out   = ",T_out
                print*, "vw:fn_vode = ", fn_vode
                print*, "vw:NR_vode = ", NR_vode
                if (e_out .lt. 0.d0) then
                    !$OMP CRITICAL
                    print *,'negative e exiting strang integration ',z, i,j,k, rho/mean_rhob, e_out
                    call flush(6)
                    !$OMP END CRITICAL
                    T_out  = 10.0
                    ne_out = 0.0
                    mu     = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_out)
                    e_out  = T_out / (gamma_minus_1 * mp_over_kB * mu)
!                    call bl_abort('bad e out of strang')
                end if

                ! Update T and ne (do not use stuff computed in f_rhs, per vode manual)
                call nyx_eos_T_given_Re(JH_vode, JHe_vode, T_out, ne_out, rho, e_out, a, species)

                ! Instanteneous heating from reionization:
                T_H = 0.0d0
                if (inhomogeneous_on .or. flash_h) then
                   if ((H_reion_z  .lt. z) .and. (H_reion_z  .ge. z_end)) T_H  = (1.0d0 - species(2))*max((T_zhi-T_out), 0.0d0)
                endif

                T_He = 0.0d0
                if (flash_he) then
                   if ((He_reion_z .lt. z) .and. (He_reion_z .ge. z_end)) T_He = (1.0d0 - species(5))*max((T_zheii-T_out), 0.0d0)
                endif

                if ((T_H .gt. 0.0d0) .or. (T_He .gt. 0.0d0)) then
                   T_out = T_out + T_H + T_He                            ! For simplicity, we assume
                   ne_out = 1.0d0 + YHELIUM                              !    completely ionized medium at
                   if (T_He .gt. 0.0d0) ne_out = ne_out + YHELIUM        !    this point.  It's a very minor
                   mu = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_out)   !    detail compared to the overall approximation.
                   e_out  = T_out / (gamma_minus_1 * mp_over_kB * mu)
                   call nyx_eos_T_given_Re(JH_vode, JHe_vode, T_out, ne_out, rho, e_out, a, species)
                endif
    !-----------------cut out end do ijk loops        
    print*, "e_out   = ",e_out
    print*, "T_out   = ",T_out
    print*, "fn_vode = ", fn_vode
    print*, "NR_vode = ", NR_vode
!    call N_VDestroy_Serial(sunvec_y)
!    call FCVodeFree(cvmem)

    deallocate(yvec)

!  call amrex_finalize()

end program main

