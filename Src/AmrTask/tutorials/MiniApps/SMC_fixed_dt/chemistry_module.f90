module chemistry_module

  implicit none

  integer, parameter :: nspecies = 9   ! number of species

  logical, save :: chemistry_initialized = .false.

  integer, private, parameter :: L_spec_name = 16 ! Each species name has at most 8 characters
  character*(L_spec_name), allocatable, save :: spec_names(:)

  double precision, save :: Ru, Ruc, Patm
  double precision, save :: molecular_weight(nspecies), inv_mwt(nspecies)

contains

  subroutine chemistry_init() bind(c, name='chemistry_init')
    integer :: iwrk
    double precision :: rwrk

    allocate(spec_names(nspecies))
    spec_names(1) = "H2"
    spec_names(2) = "O2"
    spec_names(3) = "H2O"
    spec_names(4) = "H"
    spec_names(5) = "O"
    spec_names(6) = "OH"
    spec_names(7) = "HO2"
    spec_names(8) = "H2O2"
    spec_names(9) = "N2"

    call ckrp(iwrk, rwrk, Ru, Ruc, Patm)

    call ckwt(iwrk, rwrk, molecular_weight)
    inv_mwt = 1.d0 / molecular_weight

    chemistry_initialized = .true.
  end subroutine chemistry_init

  function get_num_species () result(r) bind(c,name='get_num_species')
    integer :: r
    r = nspecies
  end function get_num_species

  subroutine chemistry_close() bind(c,name='chemistry_close')
    return
  end subroutine chemistry_close

end module chemistry_module
