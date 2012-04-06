module initialize_module

  use BoxLib
  use define_bc_module
  use bl_error_module

  implicit none

  private

  public :: initialize_bc

contains

  subroutine initialize_bc(dm,the_bc_tower,num_levs,pmask)

     use bc_module

     integer       , intent(in   ) :: dm
     type(bc_tower), intent(  out) :: the_bc_tower
     integer       , intent(in   ) :: num_levs
     logical       , intent(in   ) :: pmask(:)

     integer :: domain_phys_bc(dm,2)
     integer :: bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi

     bcx_lo = BC_PER
     bcy_lo = BC_PER
     bcz_lo = BC_PER
     bcx_hi = BC_PER
     bcy_hi = BC_PER
     bcz_hi = BC_PER

     ! Define the physical boundary conditions on the domain
     ! Put the bc values from the inputs file into domain_phys_bc
     domain_phys_bc(1,1) = bcx_lo
     domain_phys_bc(1,2) = bcx_hi
     if (pmask(1)) then
        domain_phys_bc(1,:) = BC_PER
        if (bcx_lo .ne. -1 .or. bcx_hi .ne. -1) &
             call bl_error('MUST HAVE BCX = -1 if PMASK = T')
     end if
     if (dm > 1) then
        domain_phys_bc(2,1) = bcy_lo
        domain_phys_bc(2,2) = bcy_hi
        if (pmask(2)) then
           domain_phys_bc(2,:) = BC_PER
           if (bcy_lo .ne. -1 .or. bcy_hi .ne. -1) &
                call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
        end if
     end if
     if (dm > 2) then
        domain_phys_bc(3,1) = bcz_lo
        domain_phys_bc(3,2) = bcz_hi
        if (pmask(3)) then
           domain_phys_bc(3,:) = BC_PER
           if (bcz_lo .ne. -1 .or. bcz_hi .ne. -1) &
                call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
        end if
     end if

     ! Initialize the_bc_tower object.

     call bc_tower_init(the_bc_tower,num_levs,dm,domain_phys_bc)

  end subroutine initialize_bc

end module initialize_module
