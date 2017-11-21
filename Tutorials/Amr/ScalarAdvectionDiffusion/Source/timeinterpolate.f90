


subroutine timeinterprk4_jbb(stage, lo, hi, &
     phi,  phi_lo, phi_hi, &
     old, old_lo, old_hi, &
     k1 , k1_lo, k1_hi, &
     k2 , k2_lo, k2_hi, &
     k3 , k3_lo, k3_hi, &
     k4 , k4_lo, k4_hi, &
     tf, tc_old, dt_c, dt_f &
     ) bind(C, name="timeinterprk4_jbb")

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: phi_lo(3), phi_hi(3)
  integer, intent(in) :: old_lo(3), old_hi(3)
  integer, intent(in) :: k1_lo(3), k1_hi(3)
  integer, intent(in) :: k2_lo(3), k2_hi(3)
  integer, intent(in) :: k3_lo(3), k3_hi(3)
  integer, intent(in) :: k4_lo(3), k4_hi(3), stage

  double precision, intent(in)    :: tf, tc_old, dt_c, dt_f
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
                                         phi_lo(2):phi_hi(2), &
                                         phi_lo(3):phi_hi(3))

  double precision, intent(in)    :: old(old_lo(1):old_hi(1), &
                                         old_lo(2):old_hi(2), &
                                         old_lo(3):old_hi(3))

  double precision, intent(in)    :: k1(k1_lo(1):k1_hi(1), &
                                        k1_lo(2):k1_hi(2), &
                                        k1_lo(3):k1_hi(3))
  double precision, intent(in)    :: k2(k2_lo(2):k2_hi(2), &
                                        k2_lo(2):k2_hi(2), &
                                        k2_lo(3):k2_hi(3))
  double precision, intent(in)    :: k3(k3_lo(1):k3_hi(1), &
                                        k3_lo(2):k3_hi(2), &
                                        k3_lo(3):k3_hi(3))
  double precision, intent(in)    :: k4(k4_lo(1):k4_hi(1), &
                                        k4_lo(2):k4_hi(2), &
                                        k4_lo(3):k4_hi(3))

  integer          :: i,j,k
  double precision :: k_1, k_2, k_3, k_4, squcoef, cubcoef, u0
  double precision  :: xi, fn, fnfprime, fnfdouble, dt2, dt3

  !this gets us to u0
  xi = (tf - tc_old)/dt_c

 !debug set to old method
  if(stage .eq. 0) then
     xi = (tf - tc_old)/dt_c
  else if(stage .eq. 1) then
     xi = (tf + 0.5*dt_f - tc_old)/dt_c
  else if(stage .eq. 2) then
     xi = (tf + 0.5*dt_f - tc_old)/dt_c !why, yes, this *is* the same as stage 1
  else
     xi = (tf + dt_f - tc_old)/dt_c;
  endif
  !end debug     

  dt2 = dt_f*dt_f
  dt3 = dt_f*dt_f*dt_f
  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           k_1 = k1(i,j,k)
           k_2 = k2(i,j,k)
           k_3 = k3(i,j,k)
           k_4 = k4(i,j,k)
           !straight outta mccorquodale
           squcoef = 0.5d0*(-3.0d0*k_1 + 2.0d0*k_2 + 2.0d0*k_3 - k_4)
           cubcoef = (2.0d0/3.0d0)*(k_1 - k_2 - k_3 + k_4)

           u0 = old(i,j,k) + xi*k_1 + xi*xi*squcoef + xi*xi*xi*cubcoef

           fn = k_1
           fnfprime  = (-3.0d0*k_1 + 2.0d0*k_2 + 2.0d0*k_3 - k_4)/dt_f
           fnfdouble = (-4.0d0*k_2 + k_3)/(dt_f*dt_f)
           if(stage.eq.0) then
              phi(i,j,k) = u0
           else if(stage.eq.1) then
              phi(i,j,k) = u0 + 0.5d0*dt_f*fn
           else if(stage.eq.2) then
              phi(i,j,k) = u0 + 0.5d0*dt_f*fn + 0.25d0*dt2*fnfprime + 0.0625d0*dt3*fn*fnfdouble
           else if(stage.eq.3) then
              phi(i,j,k) = u0 + 0.5d0*dt_f*fn + 0.5d0*dt2*fnfprime + 0.25d0*dt3*fn*fnfprime + 0.125d0*dt3*fn*fnfdouble
           else
              print*, "bogus stage"
              stop
           endif

!debug set to old method
           phi(i,j,k) = u0
!end debug              

        end do
     end do
  end do

end subroutine timeinterprk4_jbb

subroutine timeinterprk3_jbb(stage, lo, hi, &
     phi,  phi_lo, phi_hi, &
     old, old_lo, old_hi, &
     k1 , k1_lo, k1_hi, &
     k2 , k2_lo, k2_hi, &
     tf, tc_old, dt_c, dt_f  &
     ) bind(C, name="timeinterprk3_jbb")

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none
  double precision, intent(in)    :: tf, tc_old, dt_c, dt_f
  integer, intent(in) :: lo(3), hi(3),stage
  integer, intent(in) :: phi_lo(3), phi_hi(3)
  integer, intent(in) :: old_lo(3), old_hi(3)
  integer, intent(in) :: k1_lo(3), k1_hi(3)
  integer, intent(in) :: k2_lo(3), k2_hi(3)

  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
                                         phi_lo(2):phi_hi(2), &
                                         phi_lo(3):phi_hi(3))

  double precision, intent(in)    :: old(old_lo(1):old_hi(1), &
                                         old_lo(2):old_hi(2), &
                                         old_lo(3):old_hi(3))

  double precision, intent(in)    :: k1(k1_lo(1):k1_hi(1), &
                                        k1_lo(2):k1_hi(2), &
                                        k1_lo(3):k1_hi(3))
  double precision, intent(in)    :: k2(k2_lo(2):k2_hi(2), &
                                        k2_lo(2):k2_hi(2), &
                                        k2_lo(3):k2_hi(3))

  integer          :: i,j,k
  double precision :: k_1, k_2,  squcoef, u0
  double precision :: xi

  xi = (tf - tc_old)/dt_c

  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           k_1 = k1(i,j,k)
           k_2 = k2(i,j,k)
           !straight outta fok and rosales
           squcoef = 0.5d0*(k_2 - k_1)

           u0 = old(i,j,k) + xi*k_1 + xi*xi*squcoef 

           if(stage.eq.0) then
              phi(i,j,k) = u0
           else if(stage.eq.1) then
              phi(i,j,k) = u0 + k_1
           else if(stage.eq.2) then
              phi(i,j,k) = u0 + 0.25d0*dt_f*(k_1 + k_2)
           else
              print*, "bogus stage rk3"
              stop
           endif
        end do
     end do
  end do

end subroutine timeinterprk3_jbb


!      if(a_stage == 0)
!      {
!        xi = (tf - tc_old)/dt_c;
!      }
!      else if(a_stage == 1)
!      {
!        xi = (tf + 0.5*dt_f - tc_old)/dt_c;
!      }
!      else if(a_stage == 2)
!      {
!        xi = (tf + 0.5*dt_f - tc_old)/dt_c; //why, yes, this *is* the same as stage 1
!      }
!      else
!      {
!        xi = (tf + dt_f - tc_old)/dt_c;
!      }
subroutine timeinterprk4_mccorq(xi, lo, hi, &
     phi,  phi_lo, phi_hi, &
     old, old_lo, old_hi, &
     k1 , k1_lo, k1_hi, &
     k2 , k2_lo, k2_hi, &
     k3 , k3_lo, k3_hi, &
     k4 , k4_lo, k4_hi &
     ) bind(C, name="timeinterprk4_mccorq")

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: phi_lo(3), phi_hi(3)
  integer, intent(in) :: old_lo(3), old_hi(3)
  integer, intent(in) :: k1_lo(3), k1_hi(3)
  integer, intent(in) :: k2_lo(3), k2_hi(3)
  integer, intent(in) :: k3_lo(3), k3_hi(3)
  integer, intent(in) :: k4_lo(3), k4_hi(3)

  double precision, intent(in) :: xi
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
                                         phi_lo(2):phi_hi(2), &
                                         phi_lo(3):phi_hi(3))

  double precision, intent(in)    :: old(old_lo(1):old_hi(1), &
                                         old_lo(2):old_hi(2), &
                                         old_lo(3):old_hi(3))

  double precision, intent(in)    :: k1(k1_lo(1):k1_hi(1), &
                                        k1_lo(2):k1_hi(2), &
                                        k1_lo(3):k1_hi(3))
  double precision, intent(in)    :: k2(k2_lo(2):k2_hi(2), &
                                        k2_lo(2):k2_hi(2), &
                                        k2_lo(3):k2_hi(3))
  double precision, intent(in)    :: k3(k3_lo(1):k3_hi(1), &
                                        k3_lo(2):k3_hi(2), &
                                        k3_lo(3):k3_hi(3))
  double precision, intent(in)    :: k4(k4_lo(1):k4_hi(1), &
                                        k4_lo(2):k4_hi(2), &
                                        k4_lo(3):k4_hi(3))

  integer          :: i,j,k
  double precision :: k_1, k_2, k_3, k_4, squcoef, cubcoef, phival
  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           k_1 = k1(i,j,k)
           k_2 = k2(i,j,k)
           k_3 = k3(i,j,k)
           k_4 = k4(i,j,k)
           !straight outta mccorquodale
           squcoef = 0.5d0*(-3.0d0*k_1 + 2.0d0*k_2 + 2.0d0*k_3 - k_4)
           cubcoef = (2.0d0/3.0d0)*(k_1 - k_2 - k_3 + k_4)

           phival = old(i,j,k) + xi*k_1 + xi*xi*squcoef + xi*xi*xi*cubcoef
           phi(i,j,k) = phival

        end do
     end do
  end do

end subroutine timeinterprk4_mccorq


!      if(a_stage == 0)
!      {
!        xi = (tf - tc_old)/dt_c;
!      }
!      else if(a_stage == 1)
!      {
!        xi = (tf + dt_f - tc_old)/dt_c;
!      }
!      else if(a_stage == 2)
!      {
!        xi = (tf + 0.5*dt_f - tc_old)/dt_c; 
!      }
subroutine timeinterprk3_mccorq(xi, lo, hi, &
     phi,  phi_lo, phi_hi, &
     old, old_lo, old_hi, &
     k1 , k1_lo, k1_hi, &
     k2 , k2_lo, k2_hi &
     ) bind(C, name="timeinterpk3_mccorq")

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: phi_lo(3), phi_hi(3)
  integer, intent(in) :: old_lo(3), old_hi(3)
  integer, intent(in) :: k1_lo(3), k1_hi(3)
  integer, intent(in) :: k2_lo(3), k2_hi(3)

  double precision, intent(in) :: xi
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
                                         phi_lo(2):phi_hi(2), &
                                         phi_lo(3):phi_hi(3))

  double precision, intent(in)    :: old(old_lo(1):old_hi(1), &
                                         old_lo(2):old_hi(2), &
                                         old_lo(3):old_hi(3))

  double precision, intent(in)    :: k1(k1_lo(1):k1_hi(1), &
                                        k1_lo(2):k1_hi(2), &
                                        k1_lo(3):k1_hi(3))
  double precision, intent(in)    :: k2(k2_lo(2):k2_hi(2), &
                                        k2_lo(2):k2_hi(2), &
                                        k2_lo(3):k2_hi(3))

  integer          :: i,j,k
  double precision :: k_1, k_2,  squcoef,  phival
  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           k_1 = k1(i,j,k)
           k_2 = k2(i,j,k)
           !straight outta fok and rosales
           squcoef = 0.5d0*(k_2 - k_1)

           phival = old(i,j,k) + xi*k_1 + xi*xi*squcoef 
           phi(i,j,k) = phival

        end do
     end do
  end do

end subroutine timeinterprk3_mccorq
