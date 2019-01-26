subroutine integrate_state_fcvode_vec(lo, hi, &
                                  state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                  diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                  a, half_dt, min_iter, max_iter)
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
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   NDIAG, TEMP_COMP, NE_COMP, gamma_minus_1
    use bl_constants_module, only: M_PI
    use eos_params_module
    use network
    use eos_module, only: nyx_eos_T_given_Re, nyx_eos_T_given_Re_vec, nyx_eos_given_RT
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB
    use atomic_rates_module, only: YHELIUM
    use vode_aux_module    , only: z_vode, i_vode, j_vode, k_vode
    use cvode_interface
    use fnvector_serial
    use fcvode_extras
    use misc_params, only: simd_width
    use parallel, only : parallel_ioprocessor
    use, intrinsic :: iso_c_binding

    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in)    :: a, half_dt
    integer         , intent(inout) :: max_iter, min_iter

    integer :: i, j, k, ii
    real(rt) :: z
    real(rt), dimension(simd_width) :: rho
    real(rt), dimension(simd_width) :: T_orig, ne_orig, e_orig
    real(rt), dimension(simd_width) :: T_out , ne_out , e_out, mu
    real(rt) :: mean_rhob
    integer(c_int) :: ierr       ! error flag from C functions
    real(c_double) :: tstart     ! initial time
    real(c_double) :: rtol
    real(c_double), pointer, dimension(:) :: atol
    type(c_ptr) :: sunvec_y      ! sundials vector
    type(c_ptr) :: CVmem         ! CVODE memory
    type(c_ptr) :: sunvec_atol
    integer(c_long) :: neq
    real(c_double), pointer :: yvec(:)
    character(len=128) :: errmsg

    if (mod(hi(1)-lo(1)+1, simd_width) /= 0) then
      if (parallel_ioprocessor()) then
        !$omp single
        write(errmsg, *) "simd_width does not divide evenly to tile x-length! lo(1) = ", &
                         lo(1), " hi(1) = ", hi(1), " simd_width = ", simd_width
        call amrex_abort(errmsg)
        !$omp end single
      endif
    end if

    neq = int(simd_width, c_long)

    allocate(yvec(neq))
    allocate(atol(neq))

    z = 1.d0/a - 1.d0

    z_vode = z
    mean_rhob = comoving_OmB * 3.d0*(comoving_h*100.d0)**2 / (8.d0*M_PI*Gconst)

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT

    sunvec_y = N_VMake_Serial(NEQ, yvec)
    if (.not. c_associated(sunvec_y)) then
        call amrex_abort('integrate_state_fcvode_vec: sunvec_y = NULL')
    end if

    sunvec_atol = N_VMake_Serial(NEQ, atol)
    if (.not. c_associated(sunvec_atol)) then
        call amrex_abort('integrate_state_fcvode_vec: sunvec_atol = NULL')
    end if

    CVmem = FCVodeCreate(CV_BDF, CV_NEWTON)
    if (.not. c_associated(CVmem)) then
        call amrex_abort('integrate_state_fcvode_vec: CVmem = NULL')
    end if

    tstart = 0.0
    ! CVodeMalloc allocates variables and initialize the solver. We can
    ! initialize the solver with junk because once we enter the (i,j,k) loop we will
    ! immediately call fcvreinit which reuses the same memory allocated from
    ! CVodeMalloc but sets up new initial conditions.
    ierr = FCVodeInit(CVmem, c_funloc(RhsFn_vec), tstart, sunvec_y)
    if (ierr /= 0) then
       call amrex_abort('integrate_state_fcvode_vec: FCVodeInit() failed')
    end if

    ! Set dummy tolerances. These will be overwritten as soon as we enter the
    ! loop and reinitialize the solver.
    rtol = 1.0d-5
    atol(:) = 1.0d-10
    ierr = FCVodeSVtolerances(CVmem, rtol, sunvec_atol)
    if (ierr /= 0) then
      call amrex_abort('integrate_state_fcvode_vec: FCVodeSVtolerances() failed')
    end if

    ierr = FCVDiag(CVmem)
    if (ierr /= 0) then
       call amrex_abort('integrate_state_fcvode_vec: FCVDiag() failed')
    end if

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
          do i = lo(1),hi(1),simd_width

                ! Original values
                rho     = state(i:i+simd_width-1,j,k,URHO)
                e_orig  = state(i:i+simd_width-1,j,k,UEINT) / rho
                T_orig  = diag_eos(i:i+simd_width-1,j,k,TEMP_COMP)
                ne_orig = diag_eos(i:i+simd_width-1,j,k,  NE_COMP)

                do ii = 1, simd_width
                  if (e_orig(ii) .lt. 0.d0) then
                      !$OMP CRITICAL
                      print *,'negative e entering strang integration ',z, i+ii-1,j,k, rho(ii)/mean_rhob, e_orig(ii)
                      call bl_abort('bad e in strang')
                      !$OMP END CRITICAL
                  end if
                end do

                i_vode = i
                j_vode = j
                k_vode = k

                call fcvode_wrapper_vec(half_dt,rho,T_orig,ne_orig,e_orig,neq,CVmem,sunvec_y,yvec, &
                                              T_out ,ne_out ,e_out)

                do ii = 1, simd_width
                  if (e_out(ii) .lt. 0.d0) then
                      !$OMP CRITICAL
                      print *,'negative e exiting strang integration ',z, i,j,k, rho(ii)/mean_rhob, e_out(ii)
                      call flush(6)
                      !$OMP END CRITICAL
                      T_out(ii)  = 10.0
                      ne_out(ii) = 0.0
                      mu(ii)     = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_out(ii))
                      e_out(ii)  = T_out(ii) / (gamma_minus_1 * mp_over_kB * mu(ii))
  !                    call bl_abort('bad e out of strang')
                  end if
                end do

                ! Update (rho e) and (rho E)
                state(i:i+simd_width-1,j,k,UEINT) = state(i:i+simd_width-1,j,k,UEINT) + rho(1:simd_width) * (e_out(1:simd_width)-e_orig(1:simd_width))
                state(i:i+simd_width-1,j,k,UEDEN) = state(i:i+simd_width-1,j,k,UEDEN) + rho(1:simd_width) * (e_out(1:simd_width)-e_orig(1:simd_width))

                ! Update T and ne (do not use stuff computed in f_rhs, per vode manual)
                call nyx_eos_T_given_Re_vec(T_out(1:simd_width), ne_out(1:simd_width), rho(1:simd_width), e_out(1:simd_width), a, simd_width)
                diag_eos(i:i+simd_width-1,j,k,TEMP_COMP) = T_out(1:simd_width)
                diag_eos(i:i+simd_width-1,j,k,  NE_COMP) = ne_out(1:simd_width)

            end do ! i
        end do ! j
    end do ! k

    call N_VDestroy_Serial(sunvec_atol)
    call N_VDestroy_Serial(sunvec_y)
    call FCVodeFree(cvmem)

    deallocate(yvec)
    deallocate(atol)

end subroutine integrate_state_fcvode_vec
