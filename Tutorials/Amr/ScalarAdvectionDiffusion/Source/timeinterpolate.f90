


subroutine timeinterpolaterk4(xi, lo, hi, &
     phi,  phi_lo, phi_hi, &
     old, old_lo, old_hi, &
     k1 , k1_lo, k1_hi, &
     k2 , k2_lo, k2_hi, &
     k3 , k3_lo, k3_hi, &
     k4 , k4_lo, k4_hi &
     ) bind(C, name="timeinterpolaterk4")

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

end subroutine timeinterpolaterk4


subroutine timeinterpolaterk3(xi, lo, hi, &
     phi,  phi_lo, phi_hi, &
     old, old_lo, old_hi, &
     k1 , k1_lo, k1_hi, &
     k2 , k2_lo, k2_hi &
     ) bind(C, name="timeinterpolaterk3")

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

end subroutine timeinterpolaterk3
