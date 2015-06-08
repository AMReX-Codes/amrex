module derivative_stencil_module

  implicit none

  public

  integer, parameter :: stencil_ng = 4

  ! for 8th-order first derivatives
  double precision,dimension(4),parameter :: D8 = (/ 0.8d0, -0.2d0, 4.d0/105.d0, -1.d0/280.d0 /)

  ! coefficients for 8th-order stencil of second derivatives
  ! d(a*du/dx)/dx = H_{i+1/2} - H_{i-1/2},
  ! where H = a.M.u
  !
  double precision, save, dimension(8,8) :: M8, M8T
  !
  ! optimized for more zeros
  !  double precision, private, parameter :: M8_47 = 683.d0/10080.d0, M8_48 = -1.d0/224.d0
  !
  ! optimized for leading order truncation error assuming equal weight for the error terms
  double precision, private, parameter :: M8_47 = 3557.d0/44100.d0, M8_48 = -2083.d0/117600.d0

contains
  
  subroutine stencil_init

    ! 8th-order
    M8(1,1) = 5.d0/336.d0 + M8_48
    M8(2,1) = -11.d0/560.d0 - 2.d0*M8_48
    M8(3,1) = -1.d0/280.d0
    M8(4,1) = 17.d0/1680.d0 + 2.d0*M8_48
    M8(5,1) = -M8_48
    M8(6,1) = 0.d0
    M8(7,1) = 0.d0
    M8(8,1) = 0.d0

    M8(1,2) = -83.d0/3600.d0 - M8_47/5.d0 - 14.d0*M8_48/5.d0
    M8(2,2) = -31.d0/360.d0 + M8_47 + 3.d0*M8_48
    M8(3,2) = 1097.d0/5040.d0 - 2.d0*M8_47 + 6.d0*M8_48
    M8(4,2) = -319.d0/2520.d0 + 2.d0*M8_47 - 8.d0*M8_48
    M8(5,2) = -M8_47
    M8(6,2) = -139.d0/25200.d0 + M8_47/5.d0 + 9.d0*M8_48/5.d0
    M8(7,2) = 0.d0
    M8(8,2) = 0.d0

    M8(1,3) = 299.d0/50400.d0 + 2.d0*M8_47/5.d0 + 13.d0*M8_48/5.d0
    M8(2,3) = 41.d0/200.d0 - 9.d0*M8_47/5.d0 + 4.d0*M8_48/5.d0
    M8(3,3) = -1349.d0/10080.d0 + 3.d0*M8_47 - 12.d0*M8_48
    M8(4,3) = -919.d0/5040.d0 - 2.d0*M8_47 + 6.d0*M8_48
    M8(5,3) = 65.d0/224.d0 + 7.d0*M8_48
    M8(6,3) = -467.d0/25200.d0 + 3.d0*M8_47/5.d0 - 18.d0*M8_48/5.d0
    M8(7,3) = 503.d0/50400.d0 - M8_47/5.d0 - 4.d0*M8_48/5.d0
    M8(8,3) = 0.d0

    M8(1,4) = 17.d0/12600.d0 - M8_47/5.d0 - 4.d0*M8_48/5.d0
    M8(2,4) = -5927.d0/50400.d0 + 4.d0*M8_47/5.d0 - 9.d0*M8_48/5.d0
    M8(3,4) = -887.d0/5040.d0 - M8_47 + 6.d0*M8_48
    M8(4,4) = -445.d0/2016.d0
    M8(5,4) = -583.d0/720.d0 + M8_47 - 6.d0*M8_48
    M8(6,4) = -3613.d0/50400.d0 - 4.d0*M8_47/5.d0 + 9.d0*M8_48/5.d0
    M8(7,4) = -17.d0/600.d0 + M8_47/5.d0 + 4.d0*M8_48/5.d0
    M8(8,4) = -1.d0/1120.d0

    M8(1,5) = -M8(8,4)
    M8(2,5) = -M8(7,4)
    M8(3,5) = -M8(6,4)
    M8(4,5) = -M8(5,4)
    M8(5,5) = -M8(4,4)
    M8(6,5) = -M8(3,4)
    M8(7,5) = -M8(2,4)
    M8(8,5) = -M8(1,4)

    M8(1,6) = 0.d0
    M8(2,6) = -M8(7,3)
    M8(3,6) = -M8(6,3)
    M8(4,6) = -M8(5,3)
    M8(5,6) = -M8(4,3)
    M8(6,6) = -M8(3,3)
    M8(7,6) = -M8(2,3)
    M8(8,6) = -M8(1,3)

    M8(1,7) = 0.d0
    M8(2,7) = 0.d0
    M8(3,7) = -M8(6,2)
    M8(4,7) = -M8(5,2)
    M8(5,7) = -M8(4,2)
    M8(6,7) = -M8(3,2)
    M8(7,7) = -M8(2,2)
    M8(8,7) = -M8(1,2)

    M8(1,8) = 0.d0
    M8(2,8) = 0.d0
    M8(3,8) = 0.d0
    M8(4,8) = -M8(5,1)
    M8(5,8) = -M8(4,1)
    M8(6,8) = -M8(3,1)
    M8(7,8) = -M8(2,1)
    M8(8,8) = -M8(1,1)

    M8T = transpose(M8)
    
  end subroutine stencil_init

end module derivative_stencil_module
