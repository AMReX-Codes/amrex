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


//
// **********************************************************************
// General mapped coordinates
// **********************************************************************
//

/**
 \brief Helper functions for 3D interpolation from values on nodes/vertices of arbitrary hexahedra to particle location
*/

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void particle_interp_decomp (Array2D<Real,1,7,1,7>& a, Array1D<int,1,7>& ip, int neq)
{
    amrex::Real amult,apm,apn;
    for (int n=1; n<neq+1; n++)
    {
        ip(n) = n;
    }

    for (int i=1; i<neq; i++)
    {
        int ip1=i+1;
        int k=i;
        int ipi = ip(i);
        apm  = std::abs(a(ipi,i));

        for (int j=ip1; j<neq+1; j++)
        {
            int ipj = ip(j);
            apn=std::abs(a(ipj,i));
            if(apm < apn){
                apm = apn;
                k=j;
            }
        }

        int j=ip(k);
        ip(k)=ip(i);
        ip(i)=j;
        for (int l=ip1; l<neq+1; l++)
        {
            int n=ip(l);
            amult=a(n,i)/a(j,i);
            a(n,i)=amult;
            for (int kk = ip1; kk < neq+1; kk++)
            {
                a(n,kk) -= amult*a(j,kk);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void particle_interp_solve (Array2D<Real,1,7,1,7>& a, Array1D<Real,1,7>& b, Array1D<int,1,7>& ip, int neq)
{
    amrex::Real scr;
    for (int l=2; l<neq+1; l++)
    {
        int n = ip(l);

        for (int k=1; k<l; k++){
             b(n) -= a(n,k) * b(ip(k));
        }
    }
    b(ip(neq))=b(ip(neq))/a(ip(neq),neq);

    for (int l=1; l<neq; l++)
    {
        int j=neq-l;
        int jp1=j+1;
        int n=ip(j);
        for (int k=jp1; k<neq+1; k++){
           b(n) -= a(n,k  ) * b(ip(k));
        }
        b(n)=b(n)/a(n,j);
    }
    for (int n=1; n<neq+1; n++)
    {
        while(ip(n) != n) {
            int j=ip(n);
            scr=b(j);
            ip(n)=ip(j);
            b(j)=b(ip(j));
            b(ip(j))=scr;
            ip(j)=j;
        }
    }
}


/**
 \brief Linearly interpolates the mesh data to the particle position from node-centered data on a general mapped grid.
*/
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate_nd_mapped (const P& p,
                                const amrex::Array4<amrex::Real const>& data_arr,
                                const amrex::Array4<amrex::Real const>&  loc_arr,
                                amrex::ParticleReal * val, int M = AMREX_SPACEDIM)
{
    int icomp = 0;
    int ncomp_per_array = M;
    int num_arrays = 1;
    IntVect is_nodal = amrex::IntVect::TheUnitVector();
    linear_interpolate_to_particle_mapped(p, &data_arr, loc_arr,
                                          val, &is_nodal, icomp, ncomp_per_array, num_arrays);
}

/**
 \brief Linearly interpolates the mesh data to the particle position on a general mapped grid.
        Note that currently only node-centered data interpolation is supported in 2D and 3D.
*/

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void linear_interpolate_to_particle_mapped (const P& p,
                                            const Array4<amrex::Real const>* data_arr,
                                            const amrex::Array4<amrex::Real const>&  loc_arr,
                                            amrex::ParticleReal * val,
                                            const IntVect* is_nodal,
                                            int start_comp, int ncomp, int num_arrays)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(p, data_arr, loc_arr, val, is_nodal, start_comp, ncomp, num_arrays);
    amrex::Abort(" General mapped interpolation is not supported in 1D\n");
#else

    AMREX_ASSERT(val != nullptr);
    AMREX_ASSERT(num_arrays == 1);
    for (int d = 0; d < num_arrays; d++)
    {
        // Assertion to ensure that only nodal data is supported
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE((is_nodal[d][0] == 1),"For general mapped coordinates, interpolation is supported only for node-centered data");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE((is_nodal[d][1] == 1),"For general mapped coordinates, interpolation is supported only for node-centered data");
#if (AMREX_SPACEDIM == 3)
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE((is_nodal[d][2] == 1),"For general mapped coordinates, interpolation is supported only for node-centered data");
#endif
    }
    int n_of_data = 0;

    int i = p.idata(0);
    int j = p.idata(1);
    for (int comp = start_comp; comp < ncomp; ++comp) {

#if (AMREX_SPACEDIM == 2)
        // Value of data at surrounding nodes
        amrex::Real f1 = (data_arr[n_of_data])(i  ,j  ,0,comp);
        amrex::Real f2 = (data_arr[n_of_data])(i  ,j+1,0,comp) - f1;
        amrex::Real f3 = (data_arr[n_of_data])(i+1,j  ,0,comp) - f1;
        amrex::Real f4 = (data_arr[n_of_data])(i+1,j+1,0,comp) - f1;

        // Node locations
        amrex::Real x1 =  loc_arr(i  ,j  ,0,0);
        amrex::Real x2 =  loc_arr(i  ,j+1,0,0) - x1;
        amrex::Real x3 =  loc_arr(i+1,j  ,0,0) - x1;
        amrex::Real x4 =  loc_arr(i+1,j+1,0,0) - x1;

        amrex::Real y1 =  loc_arr(i  ,j  ,0,1);
        amrex::Real y2 =  loc_arr(i  ,j+1,0,1) - y1;
        amrex::Real y3 =  loc_arr(i+1,j  ,0,1) - y1;
        amrex::Real y4 =  loc_arr(i+1,j+1,0,1) - y1;

        // Particle location
        amrex::Real x = p.pos(0) - x1;
        amrex::Real y = p.pos(1) - y1;

        amrex::Real det = x2*x4*y3*(y2 - y4) - x3*x4*y2*(y3 - y4) - x2*x3*(y2 - y3)*y4;

        amrex::Real b = (  f4*(x2*y2*y3 - x3*y2*y3) + (f2*(x3 - x4)*y3)*y4 + f3*(-(x2*y2*y4) + x4*y2*y4) ) / det;
        amrex::Real c = ( -f2*x3*x4*y3 + f4*x2*x3*(-y2 + y3) + f3*x2*x4*(y2 - y4) + f2*x3*x4*y4          ) / det;
        amrex::Real d = (  f2*x4*y3 + f4*(x3*y2 - x2*y3) - f2*x3*y4 + f3*(-(x4*y2) + x2*y4)              ) / det;

        amrex::Real f = b*x + c*y + d*x*y;
        val[comp-start_comp] = f1 + f;

#elif (AMREX_SPACEDIM == 3)

        int k = p.idata(2);

        //
        // Ordering of particles (lo,lo,lo) (hi,lo,lo) (lo,hi,lo) (hi,hi,lo) (lo,lo,hi) (hi,lo,hi) (lo,hi,hi) (hi,hi,hi)
        //    All locations are calculated as offsets from (lo,lo,lo), which is only need 7 values are needed
        //
        Array1D<Real,0,7> x_vals;
        Array1D<Real,0,7> y_vals;
        Array1D<Real,0,7> z_vals;
        Array1D<Real,1,7> b;

        amrex::Real x1, y1, z1, f1;

        x1        = loc_arr(i  ,j  ,k  ,0);
        x_vals(1) = loc_arr(i+1,j  ,k  ,0) - x1;
        x_vals(2) = loc_arr(i  ,j+1,k  ,0) - x1;
        x_vals(3) = loc_arr(i+1,j+1,k  ,0) - x1;
        x_vals(4) = loc_arr(i  ,j  ,k+1,0) - x1;
        x_vals(5) = loc_arr(i+1,j  ,k+1,0) - x1;
        x_vals(6) = loc_arr(i  ,j+1,k+1,0) - x1;
        x_vals(7) = loc_arr(i+1,j+1,k+1,0) - x1;

        y1        = loc_arr(i  ,j  ,k  ,1);
        y_vals(1) = loc_arr(i+1,j  ,k  ,1) - y1;
        y_vals(2) = loc_arr(i  ,j+1,k  ,1) - y1;
        y_vals(3) = loc_arr(i+1,j+1,k  ,1) - y1;
        y_vals(4) = loc_arr(i  ,j  ,k+1,1) - y1;
        y_vals(5) = loc_arr(i+1,j  ,k+1,1) - y1;
        y_vals(6) = loc_arr(i  ,j+1,k+1,1) - y1;
        y_vals(7) = loc_arr(i+1,j+1,k+1,1) - y1;

        z1        = loc_arr(i  ,j  ,k  ,2);
        z_vals(1) = loc_arr(i+1,j  ,k  ,2) - z1;
        z_vals(2) = loc_arr(i  ,j+1,k  ,2) - z1;
        z_vals(3) = loc_arr(i+1,j+1,k  ,2) - z1;
        z_vals(4) = loc_arr(i  ,j  ,k+1,2) - z1;
        z_vals(5) = loc_arr(i+1,j  ,k+1,2) - z1;
        z_vals(6) = loc_arr(i  ,j+1,k+1,2) - z1;
        z_vals(7) = loc_arr(i+1,j+1,k+1,2) - z1;

        f1   = (data_arr[n_of_data])(i  ,j  ,k  ,comp);
        b(1) = (data_arr[n_of_data])(i+1,j  ,k  ,comp) - f1;
        b(2) = (data_arr[n_of_data])(i  ,j+1,k  ,comp) - f1;
        b(3) = (data_arr[n_of_data])(i+1,j+1,k  ,comp) - f1;
        b(4) = (data_arr[n_of_data])(i  ,j  ,k+1,comp) - f1;
        b(5) = (data_arr[n_of_data])(i+1,j  ,k+1,comp) - f1;
        b(6) = (data_arr[n_of_data])(i  ,j+1,k+1,comp) - f1;
        b(7) = (data_arr[n_of_data])(i+1,j+1,k+1,comp) - f1;

        Array2D<Real,1,7,1,7> a;
        Array1D<int,1,7> ip;

        int neq = 7;

        for (int n=1; n<neq+1; n++) {
            a(n,1) = x_vals(n);
            a(n,2) = y_vals(n);
            a(n,3) = z_vals(n);
            a(n,4) = x_vals(n)*y_vals(n);
            a(n,5) = x_vals(n)*z_vals(n);
            a(n,6) = z_vals(n)*y_vals(n);
            a(n,7) = x_vals(n)*y_vals(n)*z_vals(n);
        }

        particle_interp_decomp(a,ip,neq);
        particle_interp_solve(a,b,ip,neq);

        amrex::Real px, py, pz;
        px = p.pos(0) - x1;
        py = p.pos(1) - y1;
        pz = p.pos(2) - z1;

        amrex::Real f = b(1)*px + b(2)*py + b(3)*pz + b(4)*px*py + b(5)*px*pz + b(6)*py*pz + b(7)*px*py*pz;
        val[comp-start_comp] = static_cast<ParticleReal>(f + f1);
#endif

    } // comp
#endif
}

}  // namespace amrex
#endif  // include guard
