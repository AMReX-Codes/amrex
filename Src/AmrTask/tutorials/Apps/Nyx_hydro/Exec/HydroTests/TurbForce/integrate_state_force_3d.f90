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

    use turbforce_module
    use probdata_module, only: prob_lo, prob_hi, alpha, rho0, temp0
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, &
                                   TEMP_COMP, NE_COMP, small_pres, small_temp, gamma_minus_1
    use bl_constants_module, only : TWO, HALF, ZERO, M_PI
    use eos_params_module
    use eos_module, only: nyx_eos_given_RT, nyx_eos_T_given_Re
    use fundamental_constants_module
 
    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    real(rt), intent(inout) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, 2)
    real(rt), intent(in)    :: dx(3), time, a, half_dt

    integer :: i, j, k
    integer :: kx,ky,kz
    integer :: xstep,ystep,zstep
    real(rt) :: scaled_time
    real(rt) :: xpos, ypos, zpos
    real(rt) :: Lx, Ly, Lz, freqx, freqy, freqz
    real(rt) :: cosx,cosy,cosz,sinx,siny,sinz
    real(rt) :: HLx,HLy,HLz
    real(rt) :: kxd,kyd,kzd
    real(rt) :: kappa,kappaMax,Lmin,xT
    real(rt) :: f1,f2,f3
    real(rt) :: twicePi
    real(rt) :: divf, totf
    real(rt) :: rho, rho_e_orig, rho_K_res, T_orig, ne
    real(rt) :: delta, delta_re, eint0, press, small_eint

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
                ne         = diag_eos(i,j,k,  NE_COMP)

                if (rho_e_orig .lt. 0.d0) then
                    print *,'(rho e) entering strang integration negative ',i,j,k, rho_e_orig
                    call bl_abort('bad rho e in strang')
                end if

                ! Compute temperature increment and ensure that new temperature is positive
		delta = half_dt * alpha * (temp0 - T_orig) / a
		diag_eos(i,j,k,TEMP_COMP) = max(T_orig + delta, small_temp)

                ! Call EOS to get internal energy for constant equilibrium temperature
                call nyx_eos_given_RT(eint0, press, rho, temp0, ne, a)
		delta_re = half_dt * alpha * (rho*eint0 - rho_e_orig) / a

                ! Call EOS to get the internal energy floor
                call nyx_eos_given_RT(small_eint, press, rho, small_temp, ne, a)
		
                ! Update cell quantities
 		state(i,j,k,UEINT) = max(rho_e_orig + delta_re, rho*small_eint)
                state(i,j,k,UEDEN) = state(i,j,k,UEINT) + rho_K_res

		!if ((i.eq.16).and.(j.eq.16)) then
                !   print *, "temp: ", k, ne, temp0, T_orig, diag_eos(i,j,k,TEMP_COMP), delta
                !   print *, "rhoe: ", k, rho, rho*eint0, rho_e_orig, state(i,j,k,UEINT), delta_re
                !endif
            end do ! i
        end do ! j
    end do ! k

    scaled_time =  time / a

    if (scaled_time.ge.stop_forcing*forcing_time_scale_max) return

    twicePi=two*M_PI

    Lx = prob_hi(1)-prob_lo(1)
    Ly = prob_hi(2)-prob_lo(2)
    Lz = prob_hi(3)-prob_lo(3)

    Lmin = min(Lx,Ly,Lz)
    kappaMax = dble(nmodes)/Lmin + 1.0d-8
    nxmodes = nmodes*int(HALF+Lx/Lmin)
    nymodes = nmodes*int(HALF+Ly/Lmin)
    nzmodes = nmodes*int(HALF+Lz/Lmin)

    xstep = int(Lx/Lmin+HALF)
    ystep = int(Ly/Lmin+HALF)
    zstep = int(Lz/Lmin+HALF)

    HLx = Lx
    HLy = Ly
    HLz = Lz

!   write(6,*)"setup", Lx,Ly,Lz,kappaMax,nxmodes,nymodes,nzmodes
!   write(6,*)dx(1),dx(2),dx(3)
!   write (6,*) "In add_turb_forcing"
    
    do k = lo(3),hi(3)
       zpos = (dble(k) + HALF) * dx(3)

       do j = lo(2),hi(2)
          ypos =  (dble(j) + HALF) * dx(2)

          do i = lo(1),hi(1)
             xpos = (dble(i) + HALF) * dx(1)

             f1 = ZERO
             f2 = ZERO
             f3 = ZERO

!            write(6,*)"i,j,k",i,j,k,state(i,j,k,URHO)
!            write(6,*)"inital",f1,f2,f3

             do kz = mode_start*zstep, nzmodes, zstep
                kzd = dble(kz)
                freqz = twicePi*kzd*HLz
                do ky = mode_start*ystep, nymodes, ystep
                   kyd = dble(ky)
                   freqy=twicePi*kyd/HLy
                   do kx = mode_start*xstep, nxmodes, xstep
                      kxd = dble(kx)
                      kappa = sqrt( (kxd*kxd)/(Lx*Lx) + (kyd*kyd)/(Ly*Ly) + (kzd*kzd)/(Lz*Lz) )
                      freqx = twicePi*kxd/HLx
                      if (kappa.le.kappaMax) then
                         xT = cos(FTX(kx,ky,kz)*scaled_time+TAT(kx,ky,kz))

                         f1 = f1 + xT * ( FAZ(kx,ky,kz)*freqy*sin(freqx*xpos+FPZX(kx,ky,kz)) * cos(freqy*ypos+FPZY(kx,ky,kz)) * &
                                                              sin(freqz*zpos+FPZZ(kx,ky,kz)) &
                              -           FAY(kx,ky,kz)*freqz*sin(freqx*xpos+FPYX(kx,ky,kz)) * sin(freqy*ypos+FPYY(kx,ky,kz)) * &
                                                              cos(freqz*zpos+FPYZ(kx,ky,kz)) )
                         f2 = f2 + xT * ( FAX(kx,ky,kz)*freqz*sin(freqx*xpos+FPXX(kx,ky,kz)) * sin(freqy*ypos+FPXY(kx,ky,kz)) * &
                                                              cos(freqz*zpos+FPXZ(kx,ky,kz)) &
                              -           FAZ(kx,ky,kz)*freqx*cos(freqx*xpos+FPZX(kx,ky,kz)) * sin(freqy*ypos+FPZY(kx,ky,kz)) * &
                                                              sin(freqz*zpos+FPZZ(kx,ky,kz)) )
                         f3 = f3 + xT * ( FAY(kx,ky,kz)*freqx*cos(freqx*xpos+FPYX(kx,ky,kz)) * sin(freqy*ypos+FPYY(kx,ky,kz)) * &
                                                              sin(freqz*zpos+FPYZ(kx,ky,kz)) &
                              -           FAX(kx,ky,kz)*freqy*sin(freqx*xpos+FPXX(kx,ky,kz)) * cos(freqy*ypos+FPXY(kx,ky,kz)) * &
                                                              sin(freqz*zpos+FPXZ(kx,ky,kz)) ) 
                      endif
                   enddo
                enddo
             enddo

             state(i,j,k,UMX) = state(i,j,k,UMX) + half_dt * state(i,j,k,URHO)*f1 / a
             state(i,j,k,UMY) = state(i,j,k,UMY) + half_dt * state(i,j,k,URHO)*f2 / a
             state(i,j,k,UMZ) = state(i,j,k,UMZ) + half_dt * state(i,j,k,URHO)*f3 / a
!            write(6,*)i,j,k,f1,f2,f3

             state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5d0*(state(i,j,k,UMX)*state(i,j,k,UMX) + &
                                                              state(i,j,k,UMY)*state(i,j,k,UMY) + &
                                                              state(i,j,k,UMZ)*state(i,j,k,UMZ))/state(i,j,k,URHO)
          enddo
       enddo
    enddo

end subroutine integrate_state_force

