module stencil_types_module

  ! If this file is updated, don't forget to update C_to_F_MG/stencil_types.H.

   ! Cell-centered stencils
   integer, parameter :: CC_CROSS_STENCIL = 11
   integer, parameter :: HO_CROSS_STENCIL = 12
   integer, parameter :: HO_DENSE_STENCIL = 13

   ! Node-centered stencils
   integer, parameter :: ND_CROSS_STENCIL = 21
   integer, parameter :: ND_DENSE_STENCIL = 22
   integer, parameter :: ND_VATER_STENCIL = 23

end module stencil_types_module
