module amrex_particlecontainer_module

  use iso_c_binding
  use amrex_base_module
  use amrex_string_module
  use amrex_fort_module, only: amrex_particle_real
  
  implicit none

  private

  ! public routines
  public :: amrex_particlecontainer_init, amrex_particlecontainer_finalize
  public :: amrex_init_particles_one_per_cell, amrex_write_particles, amrex_particle_redistribute
  public :: amrex_get_particles, amrex_num_particles
  
  type, bind(C), public :: amrex_particle
     real(amrex_particle_real)    :: pos(AMREX_SPACEDIM) !< Position
     real(amrex_particle_real)    :: vel(AMREX_SPACEDIM) !< Particle velocity
     integer(c_int)               :: id
     integer(c_int)               :: cpu
  end type amrex_particle

  interface
     subroutine amrex_fi_new_particlecontainer (pc,amrcore) bind(c)
       import
       implicit none
       type(c_ptr) :: pc
       type(c_ptr), value :: amrcore
     end subroutine amrex_fi_new_particlecontainer

     subroutine amrex_fi_delete_particlecontainer (pc) bind(c)
       import
       implicit none
       type(c_ptr), value :: pc
     end subroutine amrex_fi_delete_particlecontainer

     subroutine amrex_fi_init_particles_one_per_cell (pc) bind(c)
       import
       implicit none
       type(c_ptr), value :: pc
     end subroutine amrex_fi_init_particles_one_per_cell

     subroutine amrex_fi_write_particles (pc, dirname, pname, is_checkpoint) bind(c)
       import
       implicit none
       type(c_ptr), value :: pc
       character(kind=c_char), intent(in) :: dirname(*)
       character(kind=c_char), intent(in) :: pname(*)
       logical, intent(in)                :: is_checkpoint
     end subroutine amrex_fi_write_particles

     subroutine amrex_fi_particle_redistribute (pc,lev_min,lev_max,ng) bind(c)
       import
       implicit none
       type(c_ptr), value :: pc
       integer(c_int), value :: lev_min, lev_max, ng
     end subroutine amrex_fi_particle_redistribute

     subroutine amrex_fi_get_particles(pc, lev, mfi, dp, np) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr),    value :: pc, mfi
       type(c_ptr)           :: dp
       integer(c_long)       :: np
     end subroutine amrex_fi_get_particles

     subroutine amrex_fi_num_particles(pc, lev, mfi, np) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr),    value :: pc, mfi
       integer(c_long)       :: np
     end subroutine amrex_fi_num_particles
     
  end interface

  type(c_ptr) :: particlecontainer = c_null_ptr
  
contains

  subroutine amrex_particlecontainer_init (amrcore)
    type(c_ptr), intent(in) :: amrcore
    call amrex_fi_new_particlecontainer(particlecontainer, amrcore)
  end subroutine amrex_particlecontainer_init

  subroutine amrex_particlecontainer_finalize ()
    call amrex_fi_delete_particlecontainer(particlecontainer)
    particlecontainer = c_null_ptr
  end subroutine amrex_particlecontainer_finalize

  subroutine amrex_init_particles_one_per_cell ()
    call amrex_fi_init_particles_one_per_cell(particlecontainer)
  end subroutine amrex_init_particles_one_per_cell

  subroutine amrex_write_particles (dirname, pname, is_checkpoint)
    character(len=*), intent(in) :: dirname
    character(len=*), intent(in) :: pname
    logical, intent(in)          :: is_checkpoint    
    call amrex_fi_write_particles(particlecontainer, &
         amrex_string_f_to_c(dirname), amrex_string_f_to_c(pname), is_checkpoint)
  end subroutine amrex_write_particles

  subroutine amrex_particle_redistribute (lev_min,lev_max,nghost)
    integer, optional, intent(in) :: lev_min, lev_max, nghost
    integer(c_int) :: default_min, default_max, default_ng
    default_min = 0
    default_max = -1
    default_ng  = 0
    if(present(lev_min)) default_min = lev_min
    if(present(lev_max)) default_max = lev_max
    if(present(nghost )) default_ng = nghost
    call amrex_fi_particle_redistribute(particlecontainer, &
         default_min, default_max, default_ng)
  end subroutine amrex_particle_redistribute

  function amrex_get_particles(lev, mfi) result(particles)
    integer(c_int),     intent(in) :: lev
    type(amrex_mfiter), intent(in) :: mfi
    type(amrex_particle),  pointer :: particles(:)
    type(c_ptr), target            :: data
    integer(c_long)                :: np
    call amrex_fi_get_particles(particlecontainer, lev, mfi%p, data, np)
    call c_f_pointer(data, particles, [np])
  end function amrex_get_particles

  function amrex_num_particles(lev, mfi) result(np)
    integer(c_int), intent(in)     :: lev
    type(amrex_mfiter), intent(in) :: mfi
    integer(c_long) :: np
    call amrex_fi_num_particles(particlecontainer, lev, mfi%p, np)
  end function amrex_num_particles
  
end module amrex_particlecontainer_module

