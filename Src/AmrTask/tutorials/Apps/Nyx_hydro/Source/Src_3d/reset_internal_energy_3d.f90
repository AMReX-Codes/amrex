      subroutine reset_internal_e(lo,hi, &
                                  u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3, &
                                  d,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                                  r,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3, &
                                  print_fortran_warnings, comoving_a) &
                                  bind(C, name="reset_internal_e")

      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, small_temp, &
                                     NDIAG, NE_COMP
      use  eos_params_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: print_fortran_warnings
      integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
      integer          :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer          :: r_l1,r_l2,r_l3,r_h1,r_h2,r_h3
      real(rt) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
      real(rt) :: d(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,NDIAG)
      real(rt) :: r(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3)

      real(rt), intent(in   ) :: comoving_a

      ! Local variables
      integer          :: i,j,k
      real(rt) :: Up, Vp, Wp, ke, rho_eint, eint_new
      real(rt) :: dummy_pres, rhoInv

      ! Reset internal energy if necessary
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

           rhoInv = 1.0d0 / u(i,j,k,URHO)
           Up     = u(i,j,k,UMX) * rhoInv
           Vp     = u(i,j,k,UMY) * rhoInv
           Wp     = u(i,j,k,UMZ) * rhoInv
           ke     = 0.5d0 * u(i,j,k,URHO) * (Up*Up + Vp*Vp + Wp*Wp)

           rho_eint = u(i,j,k,UEDEN) - ke

           ! Reset (e from e) if it's greater than 0.01% of big E.
           if (rho_eint .gt. 0.d0 .and. rho_eint / u(i,j,k,UEDEN) .gt. 1.d-6) then

              ! Create reset source so u(i,j,k,UEINT) = u(i,j,k,UEINT) + r(i,j,k) = rho_eint
               r(i,j,k) = rho_eint - u(i,j,k,UEINT)
               u(i,j,k,UEINT) = rho_eint

           ! If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
           else if (u(i,j,k,UEINT) .gt. 0.d0) then

              ! e is not updated, so reset source is zero
              r(i,j,k) = 0.d0
              u(i,j,k,UEDEN) = u(i,j,k,UEINT) + ke

           ! If not resetting and little e is negative ...
           else if (u(i,j,k,UEINT) .le. 0.d0) then

              call nyx_eos_given_RT(eint_new, dummy_pres, u(i,j,k,URHO), small_temp, &
                                    d(i,j,k,NE_COMP),comoving_a)

              if (print_fortran_warnings .gt. 0) then
                 print *,'   '
                 print *,'>>> Warning: Nyx_3d::reset_internal_energy  ',i,j,k
                 print *,'>>> ... Resetting neg. e from EOS using small_temp: ',small_temp,&
                         ' from ',u(i,j,k,UEINT)/u(i,j,k,URHO),' to ', eint_new
                 call flush(6)
              end if

              ! Create reset source so u(i,j,k,UEINT) = u(i,j,k,UEINT) + r(i,j,k) = u(i,j,k,URHO) * eint_new
              r(i,j,k) = u(i,j,k,URHO) *  eint_new - u(i,j,k,UEINT)
              u(i,j,k,UEINT) = u(i,j,k,URHO) *  eint_new

              u(i,j,k,UEDEN) = u(i,j,k,UEINT) + ke

           end if

           ! Scale reset source by 1/rho so src is in terms of e
!           r(i,j,k) = r(i,j,k) / u(i,j,k,URHO)

      enddo
      enddo
      enddo

      end subroutine reset_internal_e
