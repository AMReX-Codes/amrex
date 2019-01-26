module probdata_module

      use amrex_fort_module, only : rt => amrex_real

!     These determine the refinement criteria
      real(rt), save :: denerr,  dengrad
      real(rt), save :: velerr,  velgrad
      real(rt), save :: presserr,pressgrad
      real(rt), save :: temperr,tempgrad
      real(rt), save :: raderr,radgrad
      integer         , save :: max_denerr_lev   ,max_dengrad_lev
      integer         , save :: max_velerr_lev   ,max_velgrad_lev
      integer         , save :: max_presserr_lev, max_pressgrad_lev
      integer         , save :: max_temperr_lev, max_tempgrad_lev
      integer         , save :: max_raderr_lev, max_radgrad_lev

!     Sod variables
      real(rt), save ::  p_l, u_l, rho_l, p_r, u_r, rho_r, rhoe_l, rhoe_r, frac
      real(rt), save :: center(3)

!     These help specify which specific problem
      integer        , save ::  probtype,idir

end module probdata_module
