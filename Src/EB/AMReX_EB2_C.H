#ifndef AMREX_EB2_C_H_
#define AMREX_EB2_C_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EB2_Graph.H>

#if (AMREX_SPACEDIM == 2)
#include <AMReX_EB2_2D_C.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_EB2_3D_C.H>
#endif

namespace amrex { namespace EB2 {

#if (AMREX_SPACEDIM == 2)

void build_faces (Box const& bx, Array4<EBCellFlag> const& cell,
                  Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                  Array4<Real const> const& levset,
                  Array4<Real const> const& interx, Array4<Real const> const& intery,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  Array4<Real> const& fcx, Array4<Real> const& fcy,
                  GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  GpuArray<Real,AMREX_SPACEDIM> const& problo);

#elif (AMREX_SPACEDIM == 3)

#endif

}}

#endif
