module amrex_particlecontainer_module

  use iso_c_binding
  use amrex_base_module
  use amrex_string_module
  use amrex_fort_module, only: amrex_particle_real, amrex_long

  implicit none

  private

  ! public routines
  public :: amrex_particlecontainer_build, amrex_particlecontainer_destroy
  public :: amrex_get_next_particle_id, amrex_get_cpu
  public :: amrex_get_particle_id, amrex_set_particle_id
  public :: amrex_get_particle_cpu, amrex_set_particle_cpu

  type, bind(C), public :: amrex_particle
     real(amrex_particle_real)    :: pos(AMREX_SPACEDIM) !< Position
     real(amrex_particle_real)    :: vel(AMREX_SPACEDIM) !< Particle velocity
     integer(c_int)               :: id
     integer(c_int)               :: cpu
  end type amrex_particle

  type, public :: amrex_particlecontainer
     type(c_ptr) :: p = c_null_ptr
   contains
     procedure :: write                       => amrex_write_particles
     procedure :: redistribute                => amrex_particle_redistribute
     generic :: get_particles                 => amrex_get_particles_mfi, amrex_get_particles_i
     generic :: num_particles                 => amrex_num_particles_mfi, amrex_num_particles_i
     generic :: add_particle                  => amrex_add_particle_mfi,  amrex_add_particle_i
     procedure, private :: amrex_get_particles_mfi
     procedure, private :: amrex_num_particles_mfi
     procedure, private :: amrex_add_particle_mfi
     procedure, private :: amrex_get_particles_i
     procedure, private :: amrex_num_particles_i
     procedure, private :: amrex_add_particle_i
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_particlecontainer_destroy
#endif
  end type amrex_particlecontainer

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

     subroutine amrex_fi_get_next_particle_id (id) bind(c)
       import
       implicit none
       integer(amrex_long) id
     end subroutine amrex_fi_get_next_particle_id

     subroutine amrex_fi_get_cpu (cpu) bind(c)
       import
       implicit none
       integer(c_int) cpu
     end subroutine amrex_fi_get_cpu

     subroutine amrex_fi_get_particle_id (id, p) bind(c)
       import
       implicit none
       integer(amrex_long) id
       type(c_ptr),    value :: p
     end subroutine amrex_fi_get_particle_id

     subroutine amrex_fi_set_particle_id (id, p) bind(c)
       import
       implicit none
       integer(amrex_long) id
       type(c_ptr),    value :: p
     end subroutine amrex_fi_set_particle_id

     subroutine amrex_fi_get_particle_cpu (cpu, p) bind(c)
       import
       implicit none
       integer(c_int) cpu
       type(c_ptr),    value :: p
     end subroutine amrex_fi_get_particle_cpu

     subroutine amrex_fi_set_particle_cpu (cpu, p) bind(c)
       import
       implicit none
       integer(c_int) cpu
       type(c_ptr),    value :: p
     end subroutine amrex_fi_set_particle_cpu

     subroutine amrex_fi_write_particles (pc, dirname, pname, is_checkpoint) bind(c)
       import
       implicit none
       type(c_ptr), value :: pc
       character(kind=c_char), intent(in) :: dirname(*)
       character(kind=c_char), intent(in) :: pname(*)
       integer(c_int), value              :: is_checkpoint
     end subroutine amrex_fi_write_particles

     subroutine amrex_fi_particle_redistribute (pc,lev_min,lev_max,ng) bind(c)
       import
       implicit none
       type(c_ptr), value :: pc
       integer(c_int), value :: lev_min, lev_max, ng
     end subroutine amrex_fi_particle_redistribute

     subroutine amrex_fi_get_particles_mfi(pc, lev, mfi, dp, np) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr),    value :: pc, mfi
       type(c_ptr)           :: dp
       integer(amrex_long)   :: np
     end subroutine amrex_fi_get_particles_mfi

     subroutine amrex_fi_add_particle_mfi(pc, lev, mfi, p) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr),    value :: pc, mfi
       type(c_ptr),    value :: p
     end subroutine amrex_fi_add_particle_mfi

     subroutine amrex_fi_num_particles_mfi(pc, lev, mfi, np) bind(c)
       import
       implicit none
       integer(c_int), value :: lev
       type(c_ptr),    value :: pc, mfi
       integer(amrex_long)   :: np
     end subroutine amrex_fi_num_particles_mfi

     subroutine amrex_fi_get_particles_i(pc, lev, grid, tile, dp, np) bind(c)
       import
       implicit none
       integer(c_int), value :: lev, grid, tile
       type(c_ptr),    value :: pc
       type(c_ptr)           :: dp
       integer(amrex_long)   :: np
     end subroutine amrex_fi_get_particles_i

     subroutine amrex_fi_add_particle_i(pc, lev, grid, tile, p) bind(c)
       import
       implicit none
       integer(c_int), value :: lev, grid, tile
       type(c_ptr),    value :: pc
       type(c_ptr),    value :: p
     end subroutine amrex_fi_add_particle_i

     subroutine amrex_fi_num_particles_i(pc, lev, grid, tile, np) bind(c)
       import
       implicit none
       integer(c_int), value :: lev, grid, tile
       type(c_ptr),    value :: pc
       integer(amrex_long)   :: np
     end subroutine amrex_fi_num_particles_i

  end interface

contains

  subroutine amrex_particlecontainer_build (pc, amrcore)
    type(amrex_particlecontainer), intent(inout) :: pc
    type(c_ptr),                   intent(in)    :: amrcore
    call amrex_fi_new_particlecontainer(pc%p, amrcore)
  end subroutine amrex_particlecontainer_build

  subroutine amrex_particlecontainer_destroy (this)
    type(amrex_particlecontainer), intent(inout) :: this
    call amrex_fi_delete_particlecontainer(this%p)
    this%p = c_null_ptr
  end subroutine amrex_particlecontainer_destroy

  function amrex_get_next_particle_id() result(id)
    integer(amrex_long) :: id
    call amrex_fi_get_next_particle_id(id)
  end function amrex_get_next_particle_id

  function amrex_get_cpu() result(cpu)
    integer(c_int) :: cpu
    call amrex_fi_get_cpu(cpu)
  end function amrex_get_cpu

  subroutine amrex_get_particle_id (id, particle)
    integer(amrex_long), intent(inout) :: id
    type(amrex_particle), intent(in), target :: particle
    type(amrex_particle), pointer :: ptr
    type(c_ptr) :: dp
    ptr => particle
    dp = c_loc(ptr)
    call amrex_fi_get_particle_id(id, dp)
  end subroutine amrex_get_particle_id

  subroutine amrex_set_particle_id (id, particle)
    integer(amrex_long), intent(in) :: id
    type(amrex_particle), intent(inout), target :: particle
    type(amrex_particle), pointer :: ptr
    type(c_ptr) :: dp
    ptr => particle
    dp = c_loc(ptr)
    call amrex_fi_set_particle_id(id, dp)
  end subroutine amrex_set_particle_id

  subroutine amrex_get_particle_cpu (cpu, particle)
    integer(c_int), intent(inout) :: cpu
    type(amrex_particle), intent(in), target :: particle
    type(amrex_particle), pointer :: ptr
    type(c_ptr) :: dp
    ptr => particle
    dp = c_loc(ptr)
    call amrex_fi_get_particle_cpu(cpu, dp)
  end subroutine amrex_get_particle_cpu

  subroutine amrex_set_particle_cpu (cpu, particle)
    integer(c_int), intent(in) :: cpu
    type(amrex_particle), intent(inout), target :: particle
    type(amrex_particle), pointer :: ptr
    type(c_ptr) :: dp
    ptr => particle
    dp = c_loc(ptr)
    call amrex_fi_set_particle_cpu(cpu, dp)
  end subroutine amrex_set_particle_cpu

  subroutine amrex_write_particles (this, dirname, pname, is_checkpoint)
    class(amrex_particlecontainer), intent(inout) :: this
    character(len=*), intent(in) :: dirname
    character(len=*), intent(in) :: pname
    logical, intent(in)          :: is_checkpoint

    integer(c_int) :: is_check_flag
    if (is_checkpoint) then
       is_check_flag = 1
    else
       is_check_flag = 0
    end if

    call amrex_fi_write_particles(this%p, amrex_string_f_to_c(dirname), &
         amrex_string_f_to_c(pname), is_check_flag)
  end subroutine amrex_write_particles

  subroutine amrex_particle_redistribute (this, lev_min,lev_max,nghost)
    class(amrex_particlecontainer), intent(inout) :: this
    integer, optional, intent(in) :: lev_min, lev_max, nghost
    integer(c_int) :: default_min, default_max, default_ng
    default_min = 0
    default_max = -1
    default_ng  = 0
    if(present(lev_min)) default_min = lev_min
    if(present(lev_max)) default_max = lev_max
    if(present(nghost )) default_ng = nghost
    call amrex_fi_particle_redistribute(this%p, &
         default_min, default_max, default_ng)
  end subroutine amrex_particle_redistribute

  subroutine amrex_add_particle_mfi(this, lev, mfi, particle)
    class(amrex_particlecontainer), intent(inout) :: this
    integer(c_int),       intent(in)           :: lev
    type(amrex_mfiter),   intent(in)           :: mfi
    type(amrex_particle), intent(in), target :: particle
    type(amrex_particle), pointer :: ptr
    type(c_ptr) :: dp
    ptr => particle
    dp = c_loc(ptr)
    call amrex_fi_add_particle_mfi(this%p, lev, mfi%p, dp)
  end subroutine amrex_add_particle_mfi

  function amrex_get_particles_mfi(this, lev, mfi) result(particles)
    class(amrex_particlecontainer), intent(inout) :: this
    integer(c_int),     intent(in) :: lev
    type(amrex_mfiter), intent(in) :: mfi
    type(amrex_particle), pointer  :: particles(:)
    type(c_ptr)                    :: data
    integer(amrex_long)            :: np
    call amrex_fi_get_particles_mfi(this%p, lev, mfi%p, data, np)
    call c_f_pointer(data, particles, [np])
  end function amrex_get_particles_mfi

  function amrex_num_particles_mfi(this, lev, mfi) result(np)
    class(amrex_particlecontainer), intent(inout) :: this
    integer(c_int), intent(in)     :: lev
    type(amrex_mfiter), intent(in) :: mfi
    integer(amrex_long) :: np
    call amrex_fi_num_particles_mfi(this%p, lev, mfi%p, np)
  end function amrex_num_particles_mfi

  subroutine amrex_add_particle_i(this, lev, grid, tile, particle)
    class(amrex_particlecontainer), intent(inout) :: this
    integer(c_int),       intent(in)           :: lev, grid, tile
    type(amrex_particle), intent(in), target :: particle
    type(amrex_particle), pointer :: ptr
    type(c_ptr) :: dp
    ptr => particle
    dp = c_loc(ptr)
    call amrex_fi_add_particle_i(this%p, lev, grid, tile, dp)
  end subroutine amrex_add_particle_i

  function amrex_get_particles_i(this, lev, grid, tile) result(particles)
    class(amrex_particlecontainer), intent(inout) :: this
    integer(c_int),     intent(in) :: lev, grid, tile
    type(amrex_particle), pointer  :: particles(:)
    type(c_ptr)                    :: data
    integer(amrex_long)            :: np
    call amrex_fi_get_particles_i(this%p, lev, grid, tile, data, np)
    call c_f_pointer(data, particles, [np])
  end function amrex_get_particles_i

  function amrex_num_particles_i(this, lev, grid, tile) result(np)
    class(amrex_particlecontainer), intent(inout) :: this
    integer(c_int), intent(in)     :: lev, grid, tile
    integer(amrex_long) :: np
    call amrex_fi_num_particles_i(this%p, lev, grid, tile, np)
  end function amrex_num_particles_i

end module amrex_particlecontainer_module

