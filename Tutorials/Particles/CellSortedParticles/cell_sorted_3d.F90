module cell_sorted_particle_module
  use amrex_fort_module, only: amrex_real, amrex_particle_real
  use iso_c_binding ,    only: c_int
  
  implicit none
  private
  
  public particle_t, remove_particle_from_cell
  
  type, bind(C) :: particle_t
     real(amrex_particle_real) :: pos(3)     !< Position
     real(amrex_particle_real) :: vel(3)     !< Position
     integer(c_int)            :: id         !< Particle id
     integer(c_int)            :: cpu        !< Particle cpu
     integer(c_int)            :: sorted     !< Particle is in the right cell
     integer(c_int)            :: i          !< Particle cell x
     integer(c_int)            :: j          !< Particle cell y
     integer(c_int)            :: k          !< Particle cell z
  end type particle_t

contains
  
  subroutine remove_particle_from_cell(cell_parts, cell_np, new_np, i)
    
    use iso_c_binding, only: c_int
    
    implicit none
    
    integer(c_int), intent(in   ) :: cell_np
    integer(c_int), intent(inout) :: cell_parts(cell_np)
    integer(c_int), intent(inout) :: new_np
    integer(c_int), intent(in   ) :: i 

    cell_parts(i) = cell_parts(new_np)
    new_np = new_np - 1
        
  end subroutine remove_particle_from_cell
  
end module cell_sorted_particle_module

subroutine move_particles(particles, np, lo, hi, &
     cell_part_ids, cell_part_cnt, clo, chi, plo, dx, dt) &
     bind(c,name="move_particles")
  
  use amrex_fort_module, only: amrex_real
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use cell_sorted_particle_module, only: particle_t, remove_particle_from_cell
  
  implicit none

  type(particle_t), intent(inout), target :: particles(np)
  integer(c_int), intent(in) :: np
  integer(c_int), intent(in) :: lo(3), hi(3)
  integer(c_int), intent(in) :: clo(3), chi(3)
  type(c_ptr), intent(inout) :: cell_part_ids(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  integer(c_int), intent(inout) :: cell_part_cnt(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  real(amrex_real), intent(in) :: plo(3)
  real(amrex_real), intent(in) :: dx(3)
  real(amrex_real), intent(in) :: dt
  
  integer :: i, j, k, p, cell_np, new_np
  integer :: cell(3)
  integer(c_int), pointer :: cell_parts(:)
  type(particle_t), pointer :: part
  real(amrex_real) inv_dx(3)

  inv_dx = 1.d0/dx
  
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           cell_np = cell_part_cnt(i,j,k)
           call c_f_pointer(cell_part_ids(i,j,k), cell_parts, [cell_np])

           new_np = cell_np
           p = 1
           do while (p <= new_np)
              part => particles(cell_parts(p))
              
              ! move the particle in a straight line
              part%pos = part%pos + dt*part%vel

              ! if it has changed cells, remove from vector.
              ! otherwise continue
              cell = floor((part%pos - plo)*inv_dx)              
              if ((cell(1) /= i) .or. (cell(2) /= j) .or. (cell(3) /= k)) then
                 part%sorted = 0
                 call remove_particle_from_cell(cell_parts, cell_np, new_np, p)  
              else
                 p = p + 1
              end if
           end do

           cell_part_cnt(i,j,k) = new_np
           
        end do
     end do
  end do
  
end subroutine move_particles
