#ifndef AMREX_HYPRE_ML_ABECLAP_1D_K_H_
#define AMREX_HYPRE_ML_ABECLAP_1D_K_H_

#include <AMReX_Array.H>
#include <AMReX_Orientation.H>

namespace amrex {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hypmlabeclap_f2c_set_values (IntVect const& cell, Real* values,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dx, Real sb,
                                  GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& b,
                                  GpuArray<Array4<int const>,AMREX_SPACEDIM*2> const& bmask,
                                  IntVect const& refratio, int not_covered)
{
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hypmlabeclap_c2f (int i, int j, int k,
                       Array4<GpuArray<Real,2*AMREX_SPACEDIM+1>> const& stencil,
                       GpuArray<HYPRE_Int,AMREX_SPACEDIM>* civ, HYPRE_Int* nentries,
                       int* entry_offset, Real* entry_values,
                       Array4<int const> const& offset_from,
                       Array4<int const> const& nentries_to,
                       Array4<int const> const& offset_to,
                       GpuArray<Real,AMREX_SPACEDIM> const& dx, Real sb,
                       Array4<int const> const& offset_bx,
                       Array4<int const> const& offset_by,
                       Real const* bx, Real const* by,
                       Array4<int const> const& fine_mask,
                       IntVect const& rr, Array4<int const> const& crse_mask)
{
}

}

#endif
