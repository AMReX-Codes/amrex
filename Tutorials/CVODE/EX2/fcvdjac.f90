SUBROUTINE FCVDJAC(N, T, Y, FY, JAC, H, IPAR, RPAR, &
                         V1, V2, V3, IER)
! Fortran routine for dense user-supplied Jacobian.
      IMPLICIT NONE
!
! The following declaration specification should match C type long int.
      INTEGER*8 N, IPAR(*)
      INTEGER IER
      DOUBLE PRECISION T, Y(*), FY(*), JAC(N,*), H, RPAR(*)
      DOUBLE PRECISION V1(*), V2(*), V3(*)
!
      DOUBLE PRECISION  Y1, Y2, Y3

      Y1 = Y(1)
      Y2 = Y(2)
      Y3 = Y(3)
      JAC(1,1) = -0.04D0
      JAC(1,2) = 1.0D4 * Y3
      JAC(1,3) = 1.0D4 * Y2
      JAC(2,1) =  0.04D0
      JAC(2,2) = -1.0D4 * Y3 - 6.0D7 * Y2
      JAC(2,3) = -1.0D4 * Y2
      JAC(3,2) = 6.0D7 * Y2
!
      IER = 0
!
      RETURN
END subroutine fcvdjac
