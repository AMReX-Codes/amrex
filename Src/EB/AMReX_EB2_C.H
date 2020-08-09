#ifndef AMREX_EB2_C_H_
#define AMREX_EB2_C_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EB2_Graph.H>
#include <AMReX_Geometry.H>

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

void build_cells (Box const& bx, Array4<EBCellFlag> const& cell,
                  Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm, Real small_volfrac, 
                  Geometry const& geom, bool extend_domain_face);

#elif (AMREX_SPACEDIM == 3)

void build_faces (Box const& bx, Array4<EBCellFlag> const& cell,
                  Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                  Array4<Type_t> const& fz, Array4<Type_t> const& ex,
                  Array4<Type_t> const& ey, Array4<Type_t> const& ez,
                  Array4<Real const> const& levset, Array4<Real const> const& interx,
                  Array4<Real const> const& intery, Array4<Real const> const& interz,
                  Array4<Real> const& apx, Array4<Real> const& apy,
                  Array4<Real> const& apz, Array4<Real> const& fcx,
                  Array4<Real> const& fcy, Array4<Real> const& fcz,
                  Array4<Real> const& m2x, Array4<Real> const& m2y,
                  Array4<Real> const& m2z,
                  GpuArray<Real,AMREX_SPACEDIM> const& dx,
                  GpuArray<Real,AMREX_SPACEDIM> const& problo);

void build_cells (Box const& bx, Array4<EBCellFlag> const& cell,
                  Array4<Type_t> const& fx, Array4<Type_t> const& fy,
                  Array4<Type_t> const& fz, Array4<Real> const& apx,
                  Array4<Real> const& apy, Array4<Real> const& apz,
                  Array4<Real const> const& fcx, Array4<Real const> const& fcy,
                  Array4<Real const> const& fcz, Array4<Real const> const& m2x,
                  Array4<Real const> const& m2y, Array4<Real const> const& m2z,
                  Array4<Real> const& vfrac, Array4<Real> const& vcent,
                  Array4<Real> const& barea, Array4<Real> const& bcent,
                  Array4<Real> const& bnorm, Array4<EBCellFlag> const& ctmp,
                  Real small_volfrac, 
                  Geometry const& geom, bool extend_domain_face);

#endif

}}

#endif
