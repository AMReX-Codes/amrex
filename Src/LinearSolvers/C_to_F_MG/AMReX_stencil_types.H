#ifndef AMREX_STENCIL_TYPES_H_
#define AMREX_STENCIL_TYPES_H_

namespace amrex {

   //
   // The values here must be consistent with the values in F_MG/stencil_types.f90.
   //

   // Cell-centered stencils
   const int CC_CROSS_STENCIL = 11;
   const int HO_CROSS_STENCIL = 12;
   const int HO_DENSE_STENCIL = 13;

   // Node-centered stencils
   const int ND_CROSS_STENCIL = 21;
   const int ND_DENSE_STENCIL = 22;
   const int ND_VATER_STENCIL = 23;

}

#endif
