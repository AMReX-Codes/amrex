module define_bc_module

  implicit none
  
  type bc_tower   ! ROSE will insert "public", and that's bad!
     integer :: l
  end type bc_tower

  private

  public :: bc_tower

end module define_bc_module
