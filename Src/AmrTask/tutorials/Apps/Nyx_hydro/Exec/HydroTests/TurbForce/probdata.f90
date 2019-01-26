module probdata_module

      use amrex_fort_module, only : rt => amrex_real

!     These determine the refinement criteria
      real(rt), save ::  denerr, dengrad
      integer , save ::  max_denerr_lev, max_dengrad_lev

      integer , save ::  radiative_cooling_type

      real(rt), save ::  alpha, rho0, temp0


!     Use these in add_turb_forcing.
      real(rt), save ::  prob_lo(3), prob_hi(3)
      
!because of harald... ;)
      real(rt), save :: center(3)   
end module probdata_module
