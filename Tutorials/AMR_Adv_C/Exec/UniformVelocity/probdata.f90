module probdata_module

!     These determine the refinement criteria
      double precision, save ::    adverr, advgrad

!     Advection velocities
      double precision, save :: uadv, vadv, wadv

!     Diffusion coefficient
      double precision, save :: diff_coeff

!     These help specify which specific problem
      integer         , save :: probtype

end module probdata_module
