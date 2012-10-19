module probdata_module

!     These determine the refinement criteria
      double precision, save ::   specerr, specgrad
      integer         , save ::  max_specerr_lev   ,max_specgrad_lev

!     Diffusion coefficient
      double precision, save :: diff_coeff

!     These help specify which specific problem
      integer         , save :: probtype

end module probdata_module
