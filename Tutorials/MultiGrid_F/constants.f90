module constants

    integer, parameter :: V_CYCLE_TYPE = 1
    integer, parameter :: FMG_CYCLE_TYPE = 2
    integer, parameter :: CONST_INTERP = 1
    integer, parameter :: LIN_INTERP = 2
    integer, parameter :: SIN_RHS = 1
    integer, parameter :: RAND_RHS = 2
    integer, parameter :: DELTA_RHS = 3
    integer, parameter :: CONST_COEFFS = 1
    integer, parameter :: RAND_COEFFS = 2

    double precision, parameter :: L = 1.d0
    double precision, parameter :: pi = 4.d0*datan(1.d0)

end module constants
