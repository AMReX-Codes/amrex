



! this does the jbb algorithm for rk4
subroutine timeinterprk4_jbb(stage, lo, hi, &
     phi,  phi_lo, phi_hi, &
     old, old_lo, old_hi, &
     k1 , k1_lo, k1_hi, &
     k2 , k2_lo, k2_hi, &
     k3 , k3_lo, k3_hi, &
     k4 , k4_lo, k4_hi, &
     tf, tc_old, dt_c, dt_f, iter, nref &
     ) bind(C, name="timeinterprk4_jbb")

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none
  integer, intent(in) :: lo(3), hi(3), iter, nref
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
  double precision, intent(in)    :: k2(k2_lo(1):k2_hi(1), &
                                        k2_lo(2):k2_hi(2), &
                                        k2_lo(3):k2_hi(3))
  double precision, intent(in)    :: k3(k3_lo(1):k3_hi(1), &
                                        k3_lo(2):k3_hi(2), &
                                        k3_lo(3):k3_hi(3))
  double precision, intent(in)    :: k4(k4_lo(1):k4_hi(1), &
                                        k4_lo(2):k4_hi(2), &
                                        k4_lo(3):k4_hi(3))

  integer          :: i,j,k
  double precision :: k_1, k_2, k_3, k_4, squcoef, cubcoef, u0,  utemp(0:3)
  double precision  :: xi, fn, fnfprime, fnsqfdouble, dtc2, dtc3, dtf2, dtf3, fprimsqfn

  !this gets us to un
  xi = (tf - tc_old)/dt_c

  !print*, "*************IN JBB VERSION****************"
  dtf2 = dt_f*dt_f
  dtf3 = dt_f*dt_f*dt_f
  dtc2 = dt_c*dt_c
  dtc3 = dt_c*dt_c*dt_c
  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           !the ks get multiplied by dtc in the code
           k_1 = k1(i,j,k)
           k_2 = k2(i,j,k)
           k_3 = k3(i,j,k)
           k_4 = k4(i,j,k)
           !straight outta mccorquodale
           squcoef = 0.5d0*(-3.0d0*k_1 + 2.0d0*k_2 + 2.0d0*k_3 - k_4)
           cubcoef = (2.0d0/3.0d0)*(k_1 - k_2 - k_3 + k_4)

           u0 = old(i,j,k) + xi*k_1 + xi*xi*squcoef + xi*xi*xi*cubcoef

           fn = k_1/dt_c
           fnfprime    = (1.0d0/dtc2)*(4.0d0*k_3 - k_4 - 3.0d0*k_1)
           fnsqfdouble = (1.0d0/dtc3)*(4.0d0*k_4 + 4.0d0*k_1 + 8.0d0*k_2 - 16.0d0*k_3)
           fprimsqfn   = (4.0d0/dtc3)*(k_3 - k_2)

           utemp(0) =  u0
           utemp(1) =  u0 + 0.5d0*dt_f*fn
           utemp(2) =  u0 + 0.5d0*dt_f*fn + 0.25d0*dtf2*fnfprime + 0.0625d0*dtf3*fnsqfdouble
           utemp(3) =  u0 +       dt_f*fn +  0.5d0*dtf2*fnfprime + 0.1250d0*dtf3*fnsqfdouble + 0.25d0*dtf3*fprimsqfn


           phi(i,j,k) = utemp(stage)

!debug set to old method
!           phi(i,j,k) = u0
!end debug              

        end do
     end do
  end do

end subroutine timeinterprk4_jbb



! this does the jbb algorithm for rk3
subroutine timeinterprk3_jbb(stage, lo, hi, &
     phi,  phi_lo, phi_hi, &
     old, old_lo, old_hi, &
     k1 , k1_lo, k1_hi, &
     k2 , k2_lo, k2_hi, &
     tf, tc_old, dt_c, dt_f, iter, nref  &
     ) bind(C, name="timeinterprk3_jbb")

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none
  double precision, intent(in)    :: tf, tc_old, dt_c, dt_f
  integer, intent(in) :: lo(3), hi(3),stage, iter, nref
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
  double precision, intent(in)    :: k2(k2_lo(1):k2_hi(1), &
                                        k2_lo(2):k2_hi(2), &
                                        k2_lo(3):k2_hi(3))

  integer          :: i,j,k
  double precision :: k_1, k_2,  squcoef, u0, fn, fnfprime, utemp(0:2)
  double precision :: xi


  xi = (tf              - tc_old)/dt_c

  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           k_1 = k1(i,j,k)
           k_2 = k2(i,j,k)

           !straight outta fok and rosales
           squcoef = 0.5d0*(k_2 - k_1)

           u0 = old(i,j,k) + xi*k_1 + xi*xi*squcoef 

           fn       = k_1/dt_c
           fnfprime = (k_2 - k_1)/(dt_c*dt_c)

           utemp(0) = u0
           utemp(1) = u0 +       dt_f*fn
           utemp(2) = u0 + 0.5d0*dt_f*fn + 0.25d0*dt_f*dt_f*fnfprime

           phi(i,j,k) = utemp(stage)

        end do
     end do
  end do

end subroutine timeinterprk3_jbb


!this is the simply polynomial algorithm introduced by mccorquodale for rk4
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
subroutine timeinterprk4_simplepoly(stage, lo, hi, &
     phi,  phi_lo, phi_hi, &
     old, old_lo, old_hi, &
     k1 , k1_lo, k1_hi, &
     k2 , k2_lo, k2_hi, &
     k3 , k3_lo, k3_hi, &
     k4 , k4_lo, k4_hi, &
     tf, tc_old, dt_c, dt_f, iter, nref  &
     ) bind(C, name="timeinterprk4_simplepoly")

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none
  integer, intent(in) :: lo(3), hi(3), stage, iter, nref
  integer, intent(in) :: phi_lo(3), phi_hi(3)
  integer, intent(in) :: old_lo(3), old_hi(3)
  integer, intent(in) :: k1_lo(3), k1_hi(3)
  integer, intent(in) :: k2_lo(3), k2_hi(3)
  integer, intent(in) :: k3_lo(3), k3_hi(3)
  integer, intent(in) :: k4_lo(3), k4_hi(3)

  double precision, intent(in) ::  tf, tc_old, dt_f, dt_c
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
                                         phi_lo(2):phi_hi(2), &
                                         phi_lo(3):phi_hi(3))

  double precision, intent(in)    :: old(old_lo(1):old_hi(1), &
                                         old_lo(2):old_hi(2), &
                                         old_lo(3):old_hi(3))

  double precision, intent(in)    :: k1(k1_lo(1):k1_hi(1), &
                                        k1_lo(2):k1_hi(2), &
                                        k1_lo(3):k1_hi(3))
  double precision, intent(in)    :: k2(k2_lo(1):k2_hi(1), &
                                        k2_lo(2):k2_hi(2), &
                                        k2_lo(3):k2_hi(3))
  double precision, intent(in)    :: k3(k3_lo(1):k3_hi(1), &
                                        k3_lo(2):k3_hi(2), &
                                        k3_lo(3):k3_hi(3))
  double precision, intent(in)    :: k4(k4_lo(1):k4_hi(1), &
                                        k4_lo(2):k4_hi(2), &
                                        k4_lo(3):k4_hi(3))

  integer          :: i,j,k
  double precision :: k_1, k_2, k_3, k_4, squcoef, cubcoef, phival,xi

  if(stage .eq. 0) then
     xi = (tf - tc_old)/dt_c
  else if(stage .eq. 1) then
     xi = (tf + 0.5*dt_f - tc_old)/dt_c
  else if(stage .eq. 2) then
     xi = (tf + 0.5*dt_f - tc_old)/dt_c !why, yes, this *is* the same as stage 1
  else
     xi = (tf + dt_f - tc_old)/dt_c;
  endif

  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           k_1 = k1(i,j,k)
           k_2 = k2(i,j,k)
           k_3 = k3(i,j,k)
           k_4 = k4(i,j,k)
           !straight outta  mccorquodale
           squcoef = 0.5d0*(-3.0d0*k_1 + 2.0d0*k_2 + 2.0d0*k_3 - k_4)
           cubcoef = (2.0d0/3.0d0)*(k_1 - k_2 - k_3 + k_4)

           phival = old(i,j,k) + xi*k_1 + xi*xi*squcoef + xi*xi*xi*cubcoef
           phi(i,j,k) = phival

        end do
     end do
  end do

end subroutine timeinterprk4_simplepoly


!this is the simply polynomial algorithm introduced by mccorquodale for rk3
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
subroutine timeinterprk3_simplepoly(stage, lo, hi, &
     phi,  phi_lo, phi_hi, &
     old, old_lo, old_hi, &
     k1 , k1_lo, k1_hi, &
     k2 , k2_lo, k2_hi, &
     tf, tc_old, dt_c, dt_f, iter, nref  &
     ) bind(C, name="timeinterprk3_simplepoly")

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  implicit none
  integer, intent(in) :: lo(3), hi(3), stage, iter, nref
  integer, intent(in) :: phi_lo(3), phi_hi(3)
  integer, intent(in) :: old_lo(3), old_hi(3)
  integer, intent(in) :: k1_lo(3), k1_hi(3)
  integer, intent(in) :: k2_lo(3), k2_hi(3)

  double precision, intent(in) :: tf, tc_old, dt_c, dt_f
  double precision, intent(inout) :: phi(phi_lo(1):phi_hi(1), &
                                         phi_lo(2):phi_hi(2), &
                                         phi_lo(3):phi_hi(3))

  double precision, intent(in)    :: old(old_lo(1):old_hi(1), &
                                         old_lo(2):old_hi(2), &
                                         old_lo(3):old_hi(3))

  double precision, intent(in)    :: k1(k1_lo(1):k1_hi(1), &
                                        k1_lo(2):k1_hi(2), &
                                        k1_lo(3):k1_hi(3))
  double precision, intent(in)    :: k2(k2_lo(1):k2_hi(1), &
                                        k2_lo(2):k2_hi(2), &
                                        k2_lo(3):k2_hi(3))

  integer          :: i,j,k
  double precision :: k_1, k_2,  squcoef,  phival, xi

  if(stage .eq. 0) then
     xi = (tf              - tc_old)/dt_c
  else if(stage .eq. 1) then
     xi = (tf +       dt_f - tc_old)/dt_c
  else
     xi = (tf + 0.5d0*dt_f - tc_old)/dt_c;
  endif

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

end subroutine timeinterprk3_simplepoly
