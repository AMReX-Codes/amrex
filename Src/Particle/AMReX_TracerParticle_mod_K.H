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

//
// **********************************************************************
// Regular coordinates
// **********************************************************************
//

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
    cic_interpolate_cc(p, plo, dxi, data_arr, val, M);
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate_cc (const P& p,
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

//
// **********************************************************************
// Terrain-fitted coordinates
// **********************************************************************
//

/**
 \brief Linearly interpolates the mesh data to the particle position from cell-centered data
        on a terrain-fitted grid.
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate_mapped_z (const P& p,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                               const amrex::Array4<amrex::Real const>&  data_arr,
                               const amrex::Array4<amrex::Real const>&  height_arr,
                               amrex::ParticleReal * val, int M = AMREX_SPACEDIM)
{
    cic_interpolate_cc_mapped_z(p, plo, dxi, data_arr, height_arr, val, M);
}

/**
 \brief Linearly interpolates the mesh data to the particle position from cell-centered data
        on a terrain-fitted grid.
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate_cc_mapped_z (const P& p,
                                  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                  const amrex::Array4<amrex::Real const>&  data_arr,
                                  const amrex::Array4<amrex::Real const>&  height_arr,
                                  amrex::ParticleReal * val, int M = AMREX_SPACEDIM)
{
    int icomp = 0;
    int ncomp_per_array = M;
    int num_arrays = 1;
    IntVect is_nodal = amrex::IntVect::TheZeroVector();
    linear_interpolate_to_particle_z(p, plo, dxi, &data_arr, height_arr,
                                     val, &is_nodal, icomp, ncomp_per_array, num_arrays);
}

/**
 \brief Linearly interpolates the mesh data to the particle position from node-centered data.
        on a terrain-fitted grid.
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate_nd_mapped_z (const P& p,
                                  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                  const amrex::Array4<amrex::Real const>&  data_arr,
                                  const amrex::Array4<amrex::Real const>&  height_arr,
                                  amrex::ParticleReal * val, int M = AMREX_SPACEDIM)
{
    int icomp = 0;
    int ncomp_per_array = M;
    int num_arrays = 1;
    IntVect is_nodal = amrex::IntVect::TheUnitVector();
    linear_interpolate_to_particle_z(p, plo, dxi, &data_arr, height_arr,
                                     val, &is_nodal, icomp, ncomp_per_array, num_arrays);
}

/**
 \brief Linearly interpolates the mesh data to the particle position from face-centered data
        on a terrain-fitted grid.
        The nth component of the data_arr array is nodal in the nth direction, and cell-centered in the others.
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mac_interpolate_mapped_z (const P& p,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                               amrex::GpuArray<amrex::Array4<amrex::Real const>,AMREX_SPACEDIM> const& data_arr,
                               const amrex::Array4<amrex::Real const>&  height_arr,
                               amrex::ParticleReal * val)
{
    int icomp = 0;
    int ncomp_per_array = 1;
    int num_arrays = AMREX_SPACEDIM;
    IntVect is_nodal[AMREX_SPACEDIM];
    for (int d=0; d < AMREX_SPACEDIM; ++d) {
        is_nodal[d] = amrex::IntVect::TheZeroVector();
        is_nodal[d][d] = 1;
    }
    linear_interpolate_to_particle_z(p, plo, dxi, data_arr.data(), height_arr,
                                     val, &is_nodal[0], icomp, ncomp_per_array, num_arrays);
}

/**
 \brief Linearly interpolates the mesh data to the particle position from mesh data.
        This general form can handle an arbitrary number of Array4s, each with different staggerings
        on a terrain-fitted grid.
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void linear_interpolate_to_particle_z (const P& p,
                                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                       const Array4<amrex::Real const>* data_arr,
                                       const amrex::Array4<amrex::Real const>&  height_arr,
                                       amrex::ParticleReal * val,
                                       const IntVect* is_nodal,
                                       int start_comp, int ncomp, int num_arrays)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(p, plo, dxi, data_arr, height_arr, val, is_nodal, start_comp, ncomp, num_arrays);
    amrex::Abort(" Terrain fitted grid interpolation is not supported in 1D\n");
#else
    AMREX_ASSERT(val != nullptr);

    int ctr = 0;

    for (int d = 0; d < num_arrays; d++)
    {
        AMREX_D_TERM(amrex::Real lx = (Real(p.pos(0))-plo[0])*dxi[0] - static_cast<Real>(!is_nodal[d][0])*Real(0.5);,
                     amrex::Real ly = (Real(p.pos(1))-plo[1])*dxi[1] - static_cast<Real>(!is_nodal[d][1])*Real(0.5);,);

        int const i0 = static_cast<int>(amrex::Math::floor(lx));
        int k = 0;

#if (AMREX_SPACEDIM == 2)
        amrex::ignore_unused(ly);
        int const j  = p.idata(0);

        amrex::Real hlo_xlo = amrex::Real(0.25)
                                   * ( height_arr(i0                       , j                     , k)
                                   +   height_arr(i0 + (!(is_nodal[d][0])) , j                     , k)
                                   +   height_arr(i0                       , j + (!is_nodal[d][1]) , k)
                                   +   height_arr(i0 + (!(is_nodal[d][0])) , j + (!is_nodal[d][1]) , k) );

        amrex::Real hlo_xhi = amrex::Real(0.25)
                                   * ( height_arr(i0 + 1                       , j                    , k )
                                   +   height_arr(i0 + 1 + (!(is_nodal[d][0])) , j                    , k )
                                   +   height_arr(i0 + 1                       , j + (!is_nodal[d][1]), k )
                                   +   height_arr(i0 + 1 + (!(is_nodal[d][0])) , j + (!is_nodal[d][1]), k ) );


        amrex::Real const xint = lx - static_cast<Real>(i0);
        amrex::Real sx[] = { amrex::Real(1.) - xint, xint};
        amrex::Real height_at_px = sx[0] * hlo_xlo + sx[1] * hlo_xhi;

        int const j0 = (amrex::Real(p.pos(1)) >= height_at_px) ? j : j-1;

        int yctr = 0;
        amrex::Real ht[4];
        for (int ii=0; ii < 2; ++ii) {
            for (int jj=0; jj < 2; ++jj) {
                ht[yctr] = amrex::Real(0.25)
                                * ( height_arr(i0 + ii                       , j0 + jj , k                    )
                                +   height_arr(i0 + ii + (!(is_nodal[d][0])) , j0 + jj , k                    )
                                +   height_arr(i0 + ii                       , j0 + jj + (!is_nodal[d][1]), k )
                                +   height_arr(i0 + ii + (!(is_nodal[d][0])) , j0 + jj + (!is_nodal[d][1]), k ) );
                ++yctr;
            }
        }
        amrex::Real hint_ilo = (p.pos(1) - ht[0]) / (ht[1] - ht[0]);
        amrex::Real hint_ihi = (p.pos(1) - ht[2]) / (ht[3] - ht[2]);

        amrex::Real sy[] = { amrex::Real(1.) - hint_ilo, amrex::Real(1.) - hint_ihi,
                                               hint_ilo,                   hint_ihi};

#elif (AMREX_SPACEDIM == 3)

        int const j0 = static_cast<int>(amrex::Math::floor(ly));
        k = p.idata(0);
        amrex::Real const xint = lx - static_cast<Real>(i0);
        amrex::Real const yint = ly - static_cast<Real>(j0);
        amrex::Real sx[] = { amrex::Real(1.) - xint, xint};
        amrex::Real sy[] = { amrex::Real(1.) - yint, yint};

        amrex::Real hlo[4];
        int ilo = 0;
        amrex::Real height_at_pxy = 0.;
        for (int ii = 0; ii < 2; ++ii) {
            for (int jj = 0; jj < 2; ++jj) {
                hlo[ilo] = amrex::Real(0.125)
                * ( height_arr(i0 + ii                    , j0 + jj                    , k                    )
                +   height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj                    , k                    )
                +   height_arr(i0 + ii                    , j0 + jj + (!is_nodal[d][1]), k                    )
                +   height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj + (!is_nodal[d][1]), k                    )
                +   height_arr(i0 + ii                    , j0 + jj                    , k + (!is_nodal[d][2]))
                +   height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj                    , k + (!is_nodal[d][2]))
                +   height_arr(i0 + ii                    , j0 + jj + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
                +   height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
                );
                height_at_pxy += hlo[ilo] * sx[ii] * sy[jj];
                ++ilo;
            }
        }

        int const k0 = (amrex::Real(p.pos(2)) >= height_at_pxy ) ? k : k-1;

        int zctr = 0;
        amrex::Real ht[8];
        for (int ii = 0; ii < 2; ++ii) {
        for (int jj = 0; jj < 2; ++jj) {
        for (int kk = 0; kk < 2; ++kk) {
            ht[zctr] = amrex::Real(0.125) *
             ( height_arr(i0 + ii                    , j0 + jj                    , k0 + kk                    )
             + height_arr(i0 + ii                    , j0 + jj                    , k0 + kk + (!is_nodal[d][2]))
             + height_arr(i0 + ii                    , j0 + jj + (!is_nodal[d][1]), k0 + kk                    )
             + height_arr(i0 + ii                    , j0 + jj + (!is_nodal[d][1]), k0 + kk + (!is_nodal[d][2]))
             + height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj                    , k0 + kk                    )
             + height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj                    , k0 + kk + (!is_nodal[d][2]))
             + height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj + (!is_nodal[d][1]), k0 + kk                    )
             + height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj + (!is_nodal[d][1]), k0 + kk + (!is_nodal[d][2]))
             );
            ++zctr;
        }}}

        amrex::Real hint_ilojlo = ( p.pos(2) - ht[0] ) / (ht[1] - ht[0]);
        amrex::Real hint_ilojhi = ( p.pos(2) - ht[2] ) / (ht[3] - ht[2]);
        amrex::Real hint_ihijlo = ( p.pos(2) - ht[4] ) / (ht[5] - ht[4]);
        amrex::Real hint_ihijhi = ( p.pos(2) - ht[6] ) / (ht[7] - ht[6]);

        amrex::Real sz[] = { amrex::Real(1.) - hint_ilojlo,
                             amrex::Real(1.) - hint_ihijlo,
                             amrex::Real(1.) - hint_ilojhi,
                             amrex::Real(1.) - hint_ihijhi,
                                               hint_ilojlo,
                                               hint_ihijlo,
                                               hint_ilojhi,
                                               hint_ihijhi};
#endif
        for (int comp = start_comp; comp < ncomp; ++comp) {
            val[ctr] = amrex::ParticleReal(0.);
#if (AMREX_SPACEDIM == 2)
            int k0 = 0;
            int sy_ctr = 0;
            for (int jj = 0; jj <= 1; ++jj) {
                for (int ii = 0; ii <=1; ++ii) {
                    val[ctr] += static_cast<ParticleReal>( (data_arr[d])(i0+ii, j0+jj, k0 ,comp)*sx[ii]*sy[sy_ctr] );
                    ++sy_ctr;
                }
            }
#elif (AMREX_SPACEDIM == 3)
            int sz_ctr = 0;
            for (int kk = 0; kk <= 1; ++kk) {
                for (int jj = 0; jj <= 1; ++jj) {
                    for (int ii = 0; ii <= 1; ++ii) {
                        val[ctr] += static_cast<ParticleReal>(
                                    (data_arr[d])(i0+ii, j0+jj, k0 + kk, comp)*sx[ii]*sy[jj]*sz[sz_ctr]);
                        ++sz_ctr;
                    }
                }
            }
#endif
            ctr++;
        } // ncomp
    } // d
#endif
}


}  // namespace amrex
#endif  // include guard
