#ifndef _STENCIL_TYPES_H_
#define _STENCIL_TYPES_H_

// The values here should be consistent with the values in F_MG/stencil_tyeps.f90.

   // Cell-centered stencils
   const int CC_CROSS_STENCIL = 11;
   const int HO_CROSS_STENCIL = 12;
   const int HO_DENSE_STENCIL = 13;

   // Node-centered stencils
   const int ND_CROSS_STENCIL = 21;
   const int ND_DENSE_STENCIL = 22;

#endif
