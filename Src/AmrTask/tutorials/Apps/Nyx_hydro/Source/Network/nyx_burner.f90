module nyx_burner_module

  use amrex_fort_module, only : rt => amrex_real
  use eos_module
  use network

contains

  subroutine burner(dens, temp, Xin, ein, dt, time, Xout, eout)

    implicit none
    
    real(rt), intent(in) :: dens, temp, Xin(nspec), ein, dt, time
    real(rt), intent(out) :: Xout(nspec), eout
    
    integer  :: n
    real(rt) :: enuc, dX
    
    Xout(:) = Xin(:)
    
    enuc = 0.0_rt
    do n = 1, nspec
       dX = Xout(n)-Xin(n) 
       enuc = enuc - ebin(n) * dX
    enddo
  
    eout = ein + enuc
  
  end subroutine burner

end module nyx_burner_module
