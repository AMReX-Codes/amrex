! the network module provides the information about the species we are
! advecting:
!
! nspec      -- the number of species
! spec_names -- the name of the isotope
!
!
! This module contains two routines:
!
!  network_init()        -- initialize the isotope properties
!  network_species_index -- return the index of the species given its name
!

module network

  use bl_types

  implicit none

  integer, parameter :: nspec = 3

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)

  logical, save :: network_initialized = .false.

contains
  
  subroutine network_init()

    spec_names(1) = "carbon-12"
    spec_names(2) = "oxygen-16"
    spec_names(3) = "magnesium-24"

    short_spec_names(1) = "C12"
    short_spec_names(2) = "O16"
    short_spec_names(3) = "Mg24"

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
