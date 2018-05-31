module amrex_amrtracerparticlecontainer_module

  use iso_c_binding
  use amrex_base_module
  use amrex_string_module
  use amrex_fort_module, only: amrex_particle_real
  
  implicit none

  private

  ! public routines
  public :: amrex_amrtracerparticlecontainer_init, amrex_amrtracerparticlecontainer_finalize
  public :: amrex_init_particles_one_per_cell, amrex_write_particles, amrex_particle_redistribute
  public :: amrex_get_particles, amrex_num_particles
  
  type, bind(C), public :: amrex_tracerparticle
     real(amrex_particle_real)    :: pos(AMREX_SPACEDIM) !< Position
     real(amrex_particle_real)    :: vel(AMREX_SPACEDIM) !< Particle velocity
     integer(c_int)               :: id
     integer(c_int)               :: cpu
  end type amrex_tracerparticle

  interface
     subroutine amrex_fi_new_amrtracerparticlecontainer (tracerpc,amrcore) bind(c)
       import
       implicit none
       type(c_ptr) :: tracerpc
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_new_amrtracerparticlecontainer

     subroutine amrex_fi_delete_amrtracerparticlecontainer (tracerpc) bind(c)
       import
       implicit none
       type(c_ptr), value :: tracerpc
     end subroutine amrex_fi_delete_amrtracerparticlecontainer

     subroutine amrex_fi_init_particles_one_per_cell (tracerpc) bind(c)
       import
       implicit none
       type(c_ptr), value :: tracerpc
     end subroutine amrex_fi_init_particles_one_per_cell

     subroutine amrex_fi_write_particles (tracerpc, dirname, pname, is_checkpoint) bind(c)
       import
       implicit none
       type(c_ptr), value :: tracerpc
       character(kind=c_char), intent(in) :: dirname(*)
       character(kind=c_char), intent(in) :: pname(*)
       logical, intent(in)                :: is_checkpoint
     end subroutine amrex_fi_write_particles

     subroutine amrex_fi_particle_redistribute (tracerpc,lev_min,lev_max,ng) bind(c)
       import
       implicit none
       type(c_ptr), value :: tracerpc
       integer(c_int), value :: lev_min, lev_max, ng
     end subroutine amrex_fi_particle_redistribute

     subroutine amrex_fi_get_particles(tracerpc, lev, mfi, dp, np) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr),    value :: tracerpc, mfi
       type(c_ptr)           :: dp
       integer(c_long)       :: np
     end subroutine amrex_fi_get_particles

     subroutine amrex_fi_num_particles(tracerpc, lev, mfi, np) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr),    value :: tracerpc, mfi
       integer(c_long)       :: np
     end subroutine amrex_fi_num_particles
     
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

  subroutine amrex_init_particles_one_per_cell ()
    call amrex_fi_init_particles_one_per_cell(amrtracerparticlecontainer)
  end subroutine amrex_init_particles_one_per_cell

  subroutine amrex_write_particles (dirname, pname, is_checkpoint)
    character(len=*), intent(in) :: dirname
    character(len=*), intent(in) :: pname
    logical, intent(in)          :: is_checkpoint    
    call amrex_fi_write_particles(amrtracerparticlecontainer, &
         amrex_string_f_to_c(dirname), amrex_string_f_to_c(pname), is_checkpoint)
  end subroutine amrex_write_particles

  subroutine amrex_particle_redistribute (lev_min,lev_max,nghost)
    integer, optional, intent(in) :: lev_min, lev_max, nghost
    integer(c_int) :: min, max, ng
    min = 0
    max = -1
    ng  = 0
    if(present(lev_min)) min = lev_min
    if(present(lev_max)) max = lev_max
    if(present(nghost )) ng = nghost
    call amrex_fi_particle_redistribute(amrtracerparticlecontainer, min, max, ng)
  end subroutine amrex_particle_redistribute

  function amrex_get_particles(lev, mfi) result(particles)
    integer(c_int),     intent(in) :: lev
    type(amrex_mfiter), intent(in) :: mfi
    type(amrex_tracerparticle), pointer :: particles(:)
    type(c_ptr), target            :: data
    integer(c_long)                :: np
    call amrex_fi_get_particles(amrtracerparticlecontainer, lev, mfi%p, data, np)
    call c_f_pointer(data, particles, [np])
  end function amrex_get_particles

  function amrex_num_particles(lev, mfi) result(np)
    integer(c_int), intent(in)     :: lev
    type(amrex_mfiter), intent(in) :: mfi
    integer(c_long) :: np
    call amrex_fi_num_particles(amrtracerparticlecontainer, lev, mfi%p, np)
  end function amrex_num_particles
  
end module amrex_amrtracerparticlecontainer_module

