#ifndef AMREX_TRACERPARTICLE_MOD_K_H
#define AMREX_TRACERPARTICLE_MOD_K_H

#include <AMReX_Config.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_TracerParticles.H>
#include <cmath>

namespace amrex {

/**
 \brief Linearly interpolates the mesh data to the particle position from cell-centered data.
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate (const P& p,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                      const amrex::Array4<amrex::Real const>&  data_arr,
                      amrex::ParticleReal * val, int M = AMREX_SPACEDIM)
{
    int start_comp = 0;
    int ncomp_per_array = M;
    int num_arrays = 1;
    IntVect is_nodal = amrex::IntVect::TheZeroVector();
    linear_interpolate_to_particle (p, plo, dxi, &data_arr, val, &is_nodal, start_comp, ncomp_per_array, num_arrays);
}

/**
 \brief Linearly interpolates the mesh data to the particle position from node-centered data.
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate_nd (const P& p,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                         const amrex::Array4<amrex::Real const>&  data_arr,
                         amrex::ParticleReal * val, int M = AMREX_SPACEDIM)
{
    int start_comp = 0;
    int ncomp_per_array = M;
    int num_arrays = 1;
    IntVect is_nodal = amrex::IntVect::TheUnitVector();
    linear_interpolate_to_particle (p, plo, dxi, &data_arr, val, &is_nodal, start_comp, ncomp_per_array, num_arrays);
}

/**
 \brief Linearly interpolates the mesh data to the particle position from face-centered data.
        The nth component of the data_arr array is nodal in the nth direction, and cell-centered in the others.
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mac_interpolate (const P& p,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                      amrex::GpuArray<amrex::Array4<amrex::Real const>,AMREX_SPACEDIM> const& data_arr,
                      amrex::ParticleReal * val)
{
    int start_comp = 0;
    int ncomp_per_array = 1;
    int num_arrays = AMREX_SPACEDIM;
    IntVect is_nodal[AMREX_SPACEDIM];
    for (int d=0; d < AMREX_SPACEDIM; ++d) {
        is_nodal[d] = amrex::IntVect::TheZeroVector();
        is_nodal[d][d] = 1;
    }
    linear_interpolate_to_particle (p, plo, dxi, data_arr.data(), val, &is_nodal[0], start_comp, ncomp_per_array, num_arrays);
}



/**
 \brief Linearly interpolates the mesh data to the particle position from mesh data.
        This general form can handle an arbitrary number of Array4s, each with different staggerings.
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void linear_interpolate_to_particle (const P& p,
                                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                     const Array4<amrex::Real const>* data_arr,
                                     amrex::ParticleReal * val,
                                     const IntVect* is_nodal,
                                     int start_comp, int ncomp, int num_arrays)
{
    AMREX_ASSERT(val != nullptr);

    int ctr = 0;

    for (int d = 0; d < num_arrays; d++)
    {
        AMREX_D_TERM(amrex::Real lx = (Real(p.pos(0))-plo[0])*dxi[0] - static_cast<Real>(!is_nodal[d][0])*Real(0.5);,
                     amrex::Real ly = (Real(p.pos(1))-plo[1])*dxi[1] - static_cast<Real>(!is_nodal[d][1])*Real(0.5);,
                     amrex::Real lz = (Real(p.pos(2))-plo[2])*dxi[2] - static_cast<Real>(!is_nodal[d][2])*Real(0.5));

        // (i0,j0,k0) is the lower corner of the box needed for interpolation
        // i0 = (i-1) if particle is lower  than center of cell i
        // i0 = (i  ) if particle is higher than center of cell i
        AMREX_D_TERM(int const i0 = static_cast<int>(amrex::Math::floor(lx));,
                     int const j0 = static_cast<int>(amrex::Math::floor(ly));,
                     int const k0 = static_cast<int>(amrex::Math::floor(lz)));

        AMREX_D_TERM(amrex::Real const xint = lx - static_cast<Real>(i0);,
                     amrex::Real const yint = ly - static_cast<Real>(j0);,
                     amrex::Real const zint = lz - static_cast<Real>(k0));

        amrex::Real sx[] = {amrex::Real(1.0) - xint, xint};
#if (AMREX_SPACEDIM > 1)
        amrex::Real sy[] = {amrex::Real(1.0) - yint, yint};
#endif
#if (AMREX_SPACEDIM > 2)
        amrex::Real sz[] = {amrex::Real(1.0) - zint, zint};
#endif

        for (int comp = start_comp; comp < ncomp; ++comp) {
            val[ctr] = ParticleReal(0.0);
#if (AMREX_SPACEDIM > 2)
            for (int kk = 0; kk <=1; ++kk) {
#endif

#if (AMREX_SPACEDIM > 1)
                for (int jj = 0; jj <= 1; ++jj) {
#endif
                    for (int ii = 0; ii <= 1; ++ii) {
                        val[ctr] += static_cast<ParticleReal>((data_arr[d])(IntVect(AMREX_D_DECL(i0+ii, j0+jj, k0+kk)), comp) *
                                                                                    AMREX_D_TERM(sx[ii],*sy[jj],*sz[kk]));
            AMREX_D_TERM(},},});
            ctr++;
        } // ncomp
    } // d
}

}  // namespace amrex
#endif  // include guard
