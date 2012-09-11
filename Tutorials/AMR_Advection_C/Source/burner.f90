module burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use network

contains

  subroutine burner(dens, temp, Xin, dt, time, Xout)

    implicit none
    
    real(kind=dp_t), intent(in) :: dens, temp, Xin(nspec), dt, time
    real(kind=dp_t), intent(out) :: Xout(nspec)
    
    integer :: n
    real(kind=dp_t) :: dX
    
    Xout(:) = Xin(:)
  
  end subroutine burner

end module burner_module
