
! This module stores the runtime parameters.  
! These parameter are initialized in set_method_params().

module meth_params_module

  implicit none

  double precision, save :: difmag        ! used only in consup to weight the divu contributin
  integer         , save :: iorder        ! used only in uslope and uflaten

  integer, parameter     :: NHYP    = 4
  integer, parameter     :: MAXADV  = 2

  integer         , save :: NTHERM, NVAR
  integer         , save :: URHO, UX, UY, UZ, UFA, UFS

  integer         , save :: QTHERM, QVAR
  integer         , save :: QRHO, QU, QV, QW
  integer         , save :: QFA, QFS
  integer         , save :: nadv
  integer         , save :: normalize_species

end module meth_params_module
