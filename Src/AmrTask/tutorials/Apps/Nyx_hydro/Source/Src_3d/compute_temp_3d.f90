
      ! *************************************************************************************

      subroutine fort_compute_temp(lo,hi, &
                                   state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                   diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                   comoving_a, print_fortran_warnings) &
      bind(C, name = "fort_compute_temp")

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use atomic_rates_module, only: this_z
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UEDEN, &
                                     NDIAG, TEMP_COMP, NE_COMP, ZHI_COMP, &
                                     small_temp, heat_cool_type
      use reion_aux_module,    only: zhi_flash, zheii_flash, flash_h, flash_he, &
                                     inhomogeneous_on
      use  eos_params_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer         , intent(in   ) :: print_fortran_warnings
      real(rt), intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt), intent(in   ) :: comoving_a

      integer          :: i,j,k, JH, JHe
      real(rt) :: rhoInv,eint
      real(rt) :: ke,dummy_pres
      real(rt) :: z

      z = 1.d0/comoving_a - 1.d0

      ! Flash reionization?
      if ((flash_h .eqv. .true.) .and. (z .gt. zhi_flash)) then
         JH = 0
      else
         JH = 1
      endif
      if ((flash_he .eqv. .true.) .and. (z .gt. zheii_flash)) then
         JHe = 0
      else
         JHe = 1
      endif

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (state(i,j,k,URHO) <= 0.d0) then
                  print *,'   '
                  print *,'>>> Error: compute_temp ',i,j,k
                  print *,'>>> ... negative density ',state(i,j,k,URHO)
                  print *,'    '
                  call amrex_error("Error:: compute_temp_3d.f90 :: compute_temp")
               end if
            enddo
         enddo
      enddo

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rhoInv = 1.d0 / state(i,j,k,URHO)

               if (state(i,j,k,UEINT) > 0.d0) then

                   eint = state(i,j,k,UEINT) * rhoInv

                   JH = 1
                   if (inhomogeneous_on) then
                       if (z .gt. diag_eos(i,j,k,ZHI_COMP)) JH = 0
                   end if

                   call nyx_eos_T_given_Re(JH, JHe, diag_eos(i,j,k,TEMP_COMP), diag_eos(i,j,k,NE_COMP), &
                                           state(i,j,k,URHO), eint, comoving_a)

               else
                  if (print_fortran_warnings .gt. 0) then
                     print *,'   '
                     print *,'>>> Warning: (rho e) is negative in compute_temp: ',i,j,k
                  end if
                   ! Set temp to small_temp and compute corresponding internal energy
                   call nyx_eos_given_RT(eint, dummy_pres, state(i,j,k,URHO), small_temp, &
                                         diag_eos(i,j,k,NE_COMP), comoving_a)

                   ke = 0.5d0 * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) * rhoInv

                   diag_eos(i,j,k,TEMP_COMP) = small_temp
                   state(i,j,k,UEINT) = state(i,j,k,URHO) * eint
                   state(i,j,k,UEDEN) = state(i,j,k,UEINT) + ke

               end if

            enddo
         enddo
      enddo

      end subroutine fort_compute_temp

      subroutine fort_compute_temp_vec(lo,hi, &
                                   state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                   diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                   comoving_a, print_fortran_warnings) &
      bind(C, name = "fort_compute_temp_vec")

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use atomic_rates_module, only: this_z
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UEDEN, &
                                     NDIAG, TEMP_COMP, NE_COMP, small_temp, heat_cool_type
      use  eos_params_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer         , intent(in   ) :: print_fortran_warnings
      real(rt), intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt), intent(in   ) :: comoving_a

      integer          :: i,j,k
      real(rt) :: rhoInv,eint
      real(rt), dimension(hi(1)-lo(1)+1) :: ke,dummy_pres,small_temp_vec
      real(rt) :: z
      real(rt), dimension(hi(1)-lo(1)+1,4) :: eos_inputs_pos_ueint, eos_inputs_neg_ueint
      integer :: orig_indices(hi(1)-lo(1)+1,3)
      integer :: pos_eos_count, neg_eos_count

      z = 1.d0/comoving_a - 1.d0

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (state(i,j,k,URHO) <= 0.d0) then
                  print *,'   '
                  print *,'>>> Error: compute_temp ',i,j,k
                  print *,'>>> ... negative density ',state(i,j,k,URHO)
                  print *,'    '
                  call amrex_error("Error:: compute_temp_3d.f90 :: compute_temp")
               end if
            enddo
         enddo
      enddo

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)

            pos_eos_count = 0
            neg_eos_count = 0

            do i = lo(1),hi(1)
               rhoInv = 1.d0 / state(i,j,k,URHO)

               if (state(i,j,k,UEINT) > 0.d0) then

                   pos_eos_count = pos_eos_count + 1

                   eos_inputs_pos_ueint(pos_eos_count,1) = diag_eos(i,j,k,TEMP_COMP)
                   eos_inputs_pos_ueint(pos_eos_count,2) = diag_eos(i,j,k,NE_COMP)
                   eos_inputs_pos_ueint(pos_eos_count,3) = state(i,j,k,URHO)
                   eos_inputs_pos_ueint(pos_eos_count,4) = state(i,j,k,UEINT)*rhoInv

                   orig_indices(pos_eos_count,1) = i
                   orig_indices(pos_eos_count,2) = j
                   orig_indices(pos_eos_count,3) = k

               else

                   neg_eos_count = neg_eos_count + 1

                   eos_inputs_neg_ueint(neg_eos_count,1) = diag_eos(i,j,k,TEMP_COMP) ! DON'T NEED THIS; GET RID OF IT
                   eos_inputs_neg_ueint(neg_eos_count,2) = diag_eos(i,j,k,NE_COMP)
                   eos_inputs_neg_ueint(neg_eos_count,3) = state(i,j,k,URHO)
                   eos_inputs_neg_ueint(neg_eos_count,4) = state(i,j,k,UEINT)

                   orig_indices(neg_eos_count,1) = i
                   orig_indices(neg_eos_count,2) = j
                   orig_indices(neg_eos_count,3) = k

               end if
             end do

             ! For cells with positive E_int
             call nyx_eos_T_given_Re_vec(eos_inputs_pos_ueint(1:pos_eos_count,1), &
                                         eos_inputs_pos_ueint(1:pos_eos_count,2), &
                                         eos_inputs_pos_ueint(1:pos_eos_count,3), &
                                         eos_inputs_pos_ueint(1:pos_eos_count,4), &
                                         comoving_a, &
                                         pos_eos_count)
             diag_eos(orig_indices(1:pos_eos_count,1),j,k,TEMP_COMP) = eos_inputs_pos_ueint(1:pos_eos_count,1)
             diag_eos(orig_indices(1:pos_eos_count,1),j,k,NE_COMP)   = eos_inputs_pos_ueint(1:pos_eos_count,2)

             ! For cells with negative E_int
             call nyx_eos_given_RT_vec(eos_inputs_neg_ueint(1:neg_eos_count,4), &
                                   dummy_pres(1:neg_eos_count), &
                                   eos_inputs_neg_ueint(1:neg_eos_count,3), &
                                   small_temp_vec(1:neg_eos_count), &
                                   eos_inputs_neg_ueint(1:neg_eos_count,2), &
                                   comoving_a, &
                                   neg_eos_count)

             ke(1:neg_eos_count) = 0.5d0 * (state(orig_indices(1:neg_eos_count,1),j,k,UMX)*state(orig_indices(1:neg_eos_count,1),j,k,UMX) + &
                                   state(orig_indices(1:neg_eos_count,1),j,k,UMY)*state(orig_indices(1:neg_eos_count,1),j,k,UMY) + &
                                   state(orig_indices(1:neg_eos_count,1),j,k,UMZ)*state(orig_indices(1:neg_eos_count,1),j,k,UMZ)) * rhoInv

             diag_eos(orig_indices(1:neg_eos_count,1),j,k,TEMP_COMP) = small_temp_vec(1:neg_eos_count)
             state(orig_indices(1:neg_eos_count,1),j,k,UEINT) = eos_inputs_neg_ueint(1:neg_eos_count,3) * eos_inputs_neg_ueint(1:neg_eos_count,4)
             state(orig_indices(1:neg_eos_count,1),j,k,UEDEN) = eos_inputs_neg_ueint(1:neg_eos_count,4) + ke(1:neg_eos_count)

         enddo
      enddo

      end subroutine fort_compute_temp_vec

      subroutine fort_compute_rho_temp(lo,hi,dx, &
                                     state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                                  diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                                  rho_ave,rho_T_sum, &
                                  T_sum,Tinv_sum,T_meanrho_sum,rho_sum,vol_sum,vol_mn_sum) &
      bind(C, name = "fort_compute_rho_temp")

      use meth_params_module, only : NVAR, URHO, NDIAG, TEMP_COMP

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt), intent(in   ) :: dx(3)
      real(rt), intent(in   ) :: rho_ave
      real(rt), intent(in   ) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt), intent(inout) :: rho_T_sum, rho_sum, T_sum, Tinv_sum, T_meanrho_sum
      real(rt), intent(inout) :: vol_sum, vol_mn_sum

      integer          :: i,j,k
      real(rt) :: rho_hi, rho_lo, vol

      vol = dx(1)*dx(2)*dx(3)
      rho_hi = 1.1d0*rho_ave
      rho_lo = 0.9d0*rho_ave
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                   T_sum =     T_sum + vol*diag_eos(i,j,k,TEMP_COMP)
                Tinv_sum =  Tinv_sum + state(i,j,k,URHO)/diag_eos(i,j,k,TEMP_COMP)
               rho_T_sum = rho_T_sum + state(i,j,k,URHO)*diag_eos(i,j,k,TEMP_COMP)
                 rho_sum =   rho_sum + state(i,j,k,URHO)
                 if ( (state(i,j,k,URHO) .lt. rho_hi) .and. &
                      (state(i,j,k,URHO) .gt. rho_lo) .and. &
                      (diag_eos(i,j,k,TEMP_COMP) .le. 1.0e5) ) then
                         T_meanrho_sum = T_meanrho_sum + vol*dlog10(diag_eos(i,j,k,TEMP_COMP))
                         vol_mn_sum = vol_mn_sum + vol
                 endif
                 vol_sum = vol_sum + vol
            enddo
         enddo
      enddo

      end subroutine fort_compute_rho_temp

      subroutine fort_compute_gas_frac(lo,hi,dx, &
                                       state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                                       diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                                       rho_ave, T_cut, rho_cut, &
                                       whim_mass, whim_vol, hh_mass, &
                                       hh_vol, igm_mass, igm_vol, mass_sum, vol_sum) &
      bind(C, name = "fort_compute_gas_frac")

      use meth_params_module, only : NVAR, URHO, NDIAG, TEMP_COMP

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt), intent(in   ) :: dx(3)
      real(rt), intent(in   ) :: rho_ave, T_cut, rho_cut
      real(rt), intent(in   ) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt), intent(inout) :: whim_mass, whim_vol, hh_mass, hh_vol, igm_mass, igm_vol
      real(rt), intent(inout) :: mass_sum, vol_sum

      integer :: i,j,k
      real(rt) :: vol, T, R

      vol = dx(1)*dx(2)*dx(3)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                 T = diag_eos(i,j,k,TEMP_COMP)
                 R = state(i,j,k,URHO) / rho_ave
                 if ( (T .gt. T_cut) .and. (R .le. rho_cut) ) then
                     whim_mass = whim_mass + state(i,j,k,URHO)*vol
                     whim_vol  = whim_vol  + vol
                 else if ( (T .gt. T_cut) .and. (R .gt. rho_cut) ) then
                     hh_mass = hh_mass + state(i,j,k,URHO)*vol
                     hh_vol  = hh_vol  + vol
                 else if ( (T .le. T_cut) .and. (R .le. rho_cut) ) then
                     igm_mass = igm_mass + state(i,j,k,URHO)*vol
                     igm_vol  = igm_vol  + vol
                 endif
                 mass_sum = mass_sum + state(i,j,k,URHO)*vol
                 vol_sum  = vol_sum + vol
            enddo
         enddo
      enddo

      end subroutine fort_compute_gas_frac

      subroutine fort_compute_max_temp_loc(lo,hi, &
                                           state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                           diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                           max_temp, den_maxt, imax, jmax, kmax) &
      bind(C, name = "fort_compute_max_temp_loc")

      use meth_params_module, only : TEMP_COMP, NVAR, URHO, NDIAG

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      real(rt), intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      real(rt), intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt), intent(in   ) :: max_temp
      real(rt), intent(  out) :: den_maxt
      integer         , intent(inout) :: imax,jmax,kmax

      integer                         :: i,j,k
      real(rt)                :: one_minus_eps

      one_minus_eps = 1.d0 - 1.d-12

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (diag_eos(i,j,k,TEMP_COMP) .ge. one_minus_eps*max_temp) then
                  imax = i
                  jmax = j
                  kmax = k
                  den_maxt = state(i,j,k,URHO)
               end if
            enddo
         enddo
      enddo

      end subroutine fort_compute_max_temp_loc
