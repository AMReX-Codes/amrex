module probdata_module

      use amrex_fort_module, only : rt => amrex_real

!     These determine the refinement criteria
      real(rt), save ::    denerr,  dengrad
      real(rt), save ::    velerr,  velgrad
      real(rt), save ::  presserr,pressgrad
      real(rt), save ::   temperr, tempgrad
      real(rt), save ::    raderr,  radgrad
      integer , save ::  max_denerr_lev   ,max_dengrad_lev
      integer , save ::  max_velerr_lev   ,max_velgrad_lev
      integer , save ::  max_presserr_lev, max_pressgrad_lev
      integer , save ::  max_temperr_lev,  max_tempgrad_lev
      integer , save ::  max_raderr_lev,   max_radgrad_lev

      real(rt), save ::  center(3)

!     Sod variables
      real(rt), save ::  p_ambient, dens_ambient, exp_energy
      real(rt), save ::  r_init
      integer , save ::  nsub

!     These help specify which specific problem
      integer         , save :: probtype,idir

end module probdata_module
