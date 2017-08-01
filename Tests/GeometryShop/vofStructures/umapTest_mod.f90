#include <umapTest.H>

module umap_test_module

  use amrex_fort_module, only : dim=>bl_spacedim, amrex_real
  implicit none

  type, bind(c) :: facedata
     real(amrex_real) :: aperature
     real(amrex_real) :: centroid(0:dim-1)
  end type facedata

  type, bind(c) :: ebbndrydata
     real(amrex_real) :: normal(0:dim-1)
     real(amrex_real) :: bndry_centroid(0:dim-1)
     real(amrex_real) :: value
     real(amrex_real) :: bndry_area
     real(amrex_real) :: vol_frac
  end type ebbndrydata

end module umap_test_module

