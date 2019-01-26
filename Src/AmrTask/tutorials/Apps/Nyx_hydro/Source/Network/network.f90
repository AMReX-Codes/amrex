! The simple network is used for a gamma-law gas with H and He. 
!
! The network module provides the information about the species we are
! advecting:
!
! nspec      -- the number of species
!
! aion       -- atomic number
! zion       -- proton number
!
! spec_names       -- the           name of the isotope
! short_spec_names -- the shortened name of the isotope
!
!
! This module contains two routines:
!
!  network_init()        -- initialize the isotope properties
!
!  network_species_index -- return the index of the species given its name
!

module network

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: nspec = 2
  integer, parameter :: naux  = 0

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(kind=rt), save :: aion(nspec), zion(nspec), ebin(nspec)

  logical, save :: network_initialized = .false.

contains
  
  subroutine network_init()

    spec_names(1) = "H"
    spec_names(2) = "He"    

    short_spec_names(1) = "H"
    short_spec_names(2) = "He"

    aion(1) = 1.0_rt
    aion(2) = 4.0_rt
    
    zion(1) = 1.0_rt
    zion(2) = 2.0_rt

    network_initialized = .true.

  end subroutine network_init
  
  function network_species_index(name)

    character(len=*) :: name
    integer :: network_species_index, n

    network_species_index = -1

    do n = 1, nspec
       if (name == spec_names(n)) then
          network_species_index = n
          exit
       endif
    enddo
    
    return
  end function network_species_index

end module network
