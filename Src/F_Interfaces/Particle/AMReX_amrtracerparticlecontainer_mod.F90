module amrex_amrtracerparticlecontainer_module

  use iso_c_binding
  use amrex_base_module
 
  implicit none

  private

  ! public routines
  public :: amrex_amrtracerparticlecontainer_init, &
       amrex_amrtracerparticlecontainer_finalize

  interface
     subroutine amrex_fi_new_amrtracerparticlecontainer (tracerpc,amrcore) bind(c)
       import
       implicit none
       type(c_ptr) :: tracerpc, amrcore
     end subroutine amrex_fi_new_amrtracerparticlecontainer

     subroutine amrex_fi_delete_amrtracerparticlecontainer (tracerpc) bind(c)
       import
       implicit none
       type(c_ptr), value :: tracerpc
     end subroutine amrex_fi_delete_amrtracerparticlecontainer
  end interface

  type(c_ptr) :: amrtracerparticlecontainer = c_null_ptr
  
contains

  subroutine amrex_amrtracerparticlecontainer_init (amrcore)
    type(c_ptr), intent(in) :: amrcore
    call amrex_fi_new_amrtracerparticlecontainer(amrtracerparticlecontainer, amrcore)
  end subroutine amrex_amrtracerparticlecontainer_init

  subroutine amrex_amrtracerparticlecontainer_finalize ()
    call amrex_fi_delete_amrtracerparticlecontainer(amrtracerparticlecontainer)
    amrtracerparticlecontainer = c_null_ptr
  end subroutine amrex_amrtracerparticlecontainer_finalize

end module amrex_amrtracerparticlecontainer_module

