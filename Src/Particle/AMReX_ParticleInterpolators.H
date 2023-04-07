#ifndef AMREX_PARTICLEINTERPOLATORS_H_
#define AMREX_PARTICLEINTERPOLATORS_H_
#include <AMReX_Config.H>

#include <AMReX_IntVect.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>

namespace amrex::ParticleInterpolator
{

/** \brief A base class for doing general particle/mesh interpolation operations.
 *
 *  Note that you don't call the base class version of this directly,
 *  you use one of the derived versions below that implement specific interpolations.
 */
template <class Derived, class WeightType>
struct Base
{
    int index[3];
    WeightType* w;

    /** \brief A general function for doing particle-to-mesh interpolation for one particle.
     *
     *  Note that you don't call the base class version of this, you call one of the derived
     *  versions below that implement specific interpolations.
     *
     * \tparam P the particle type
     * \tparam V the mesh data type (i.e. float, double, int)
     * \tparam F callable that generates the particle quantity to interpolate
     *
     * \param p the particle to interpolate
     * \param arr the Array4 to interpolate to
     * \param src_comp the particle component to start at
     * \param dst_comp the mesh component to start at
     * \param num_comps the number of components to interpolate
     * \param f function for computing the particle quantity to interpolate
     *
     * Usage:
     *
     *    Interpolate particle comp 0 to mesh comp 0 with no weighting:
     * \code{.cpp}
     *        interp.ParticleToMesh(p, rho, 0, 0, 1,
     *                [=] AMREX_GPU_DEVICE (const MyPC::ParticleType& part, int comp)
     *                {
     *                     return part.rdata(comp);  // no weighting
     *                 });
     * \endcode
     *
     * Usage:
     *
     *    Interpolate 3 particle components 1-3 to mesh components 1-3, weighting
     *    by particle comp 0:
     * \code{.cpp}
     *        interp.ParticleToMesh(p, rho, 1, 1, 3,
     *                [=] AMREX_GPU_DEVICE (const MyPC::ParticleType& part, int comp)
     *                {
     *                     return part.rdata(0)*part.rdata(comp);  // weight by comp0
     *                 });
     * \endcode
     */
    template <typename P, typename V, typename F>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void ParticleToMesh (const P& p,
                         amrex::Array4<V> const& arr,
                         int src_comp, int dst_comp, int num_comps, F&& f)
    {
        static constexpr int stencil_width = Derived::stencil_width;
        for (int ic=0; ic < num_comps; ++ic) {
            for (int kk = 0; kk <= Derived::nz; ++kk) {
                for (int jj = 0; jj <= Derived::ny; ++jj) {
                    for (int ii = 0; ii <= Derived::nx; ++ii) {
                        const auto pval = f(p, src_comp+ic);
                        const auto val = w[0*stencil_width+ii] *
                                         w[1*stencil_width+jj] *
                                         w[2*stencil_width+kk] * pval;
                        Gpu::Atomic::AddNoRet(&arr(index[0]+ii, index[1]+jj, index[2]+kk, ic+dst_comp), val);
                    }
                }
            }
        }
    }

    /** \brief A general function for doing mesh-to-particle interpolation for one particle.
     *
     *  Note that you don't call the base class version of this, you call one of the derived
     *  versions below that implement specific interpolations.
     *
     * \tparam P the particle type
     * \tparam V the mesh data type (i.e. float, double, int)
     * \tparam F callable that generates the mesh quantity to interpolate
     * \tparam G callable that updates the particle given the mesh value
     *
     * \param p the particle to interpolate
     * \param arr the Array4 to interpolate to
     * \param src_comp the particle component to start at
     * \param dst_comp the mesh component to start at
     * \param num_comps the number of components to interpolate
     * \param f function for computing the particle quantity to interpolate
     * \tparam g function that updates the particle given the mesh value
     *
     * Usage:
     *
     *    Interpolate mesh comps 0-2 to particle comps 4-6 with no weighting
     *    using addition:
     * \code{.cpp}
     *    interp.MeshToParticle(p, acc, 0, 4, 3,
     *            [=] AMREX_GPU_DEVICE (amrex::Array4<const amrex::Real> const& arr,
     *                                  int i, int j, int k, int comp)
     *            {
     *                return arr(i, j, k, comp);  // no weighting
     *            },
     *            [=] AMREX_GPU_DEVICE (MyParticleContainer::ParticleType& part,
     *                                  int comp, amrex::Real val)
     *            {
     *                part.rdata(comp) += val;
     *            });
     * \endcode
     *
     * Usage:
     *
     *    Interpolate mesh comp 0 to particle comp 0, simply setting the value for
     *    the result instead of adding:
     * \code{.cpp}
     *        interp.MeshToParticle(p, count, 0, 0, 1,
     *            [=] AMREX_GPU_DEVICE (amrex::Array4<const int> const& arr,
     *                                  int i, int j, int k, int comp)
     *            {
     *                return arr(i, j, k, comp);  // no weighting
     *            },
     *            [=] AMREX_GPU_DEVICE (MyParticleContainer::ParticleType& part,
     *                                  int comp, int val)
     *            {
     *                part.idata(comp) = val;
     *            });
     * \endcode
     */
    template <typename P, typename V, typename F, typename G>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void MeshToParticle (P& p,
                         amrex::Array4<const V> const& arr,
                         int src_comp, int dst_comp, int num_comps, F&& f, G&& g)
    {
        static constexpr int stencil_width = Derived::stencil_width;
        for (int ic=0; ic < num_comps; ++ic) {
            for (int kk = 0; kk <= Derived::nz; ++kk) {
                for (int jj = 0; jj <= Derived::ny; ++jj) {
                    for (int ii = 0; ii <= Derived::nx; ++ii) {
                        const auto mval = f(arr,index[0]+ii,index[1]+jj,index[2]+kk,src_comp+ic);
                        const auto val = w[0*stencil_width+ii] *
                                         w[1*stencil_width+jj] *
                                         w[2*stencil_width+kk] * mval;
                        g(p, ic + dst_comp, val);
                    }
                }
            }
        }
    }
};

/** \brief A class the implements nearest grid point particle/mesh interpolation.
 *
 *   Usage:
 *   \code{.cpp}
 *      ParticleInterpolator::Nearest interp(p, plo, dxi);
 *
 *        interp.MeshToParticle(p, count, 0, 0, 1,
 *                [=] AMREX_GPU_DEVICE (amrex::Array4<const int> const& arr,
 *                                      int i, int j, int k, int comp)
 *                {
 *                    return arr(i, j, k, comp);  // no weighting
 *                },
 *                [=] AMREX_GPU_DEVICE (MyParticleContainer::ParticleType& part,
 *                                      int comp, int val)
 *                {
 *                    part.idata(comp) = val;
 *                });
 *   \endcode
 */
struct Nearest : public Base<Nearest, int>
{
    static constexpr int stencil_width = 1;
    int weights[3*stencil_width];

    static constexpr int nx = (AMREX_SPACEDIM >= 1) ? stencil_width - 1 : 0;
    static constexpr int ny = (AMREX_SPACEDIM >= 2) ? stencil_width - 1 : 0;
    static constexpr int nz = (AMREX_SPACEDIM >= 3) ? stencil_width - 1 : 0;

    template <typename P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Nearest (const P& p,
             amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
             amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
    {
        w = &weights[0];
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            amrex::Real l = (p.pos(i) - plo[i]) * dxi[i] + 0.5;
            index[i] = static_cast<int>(amrex::Math::floor(l));
            w[i] = 1;
        }
        for (int i = AMREX_SPACEDIM; i < 3; ++i) {
            index[i] = 0;
            w[i] = 1;
        }
    }
};

/** \brief A class the implements linear (CIC) particle/mesh interpolation.
 *
 *   Usage:
 *   \code{.cpp}
 *        ParticleInterpolator::Linear interp(p, plo, dxi);
 *
 *        interp.ParticleToMesh(p, rho, 0, 0, 1,
 *                    [=] AMREX_GPU_DEVICE (const MyPC::ParticleType& part, int comp)
 *                    {
 *                        return part.rdata(comp);  // no weighting
 *                    });
 *   \endcode
 */
struct Linear : public Base<Linear, amrex::Real>
{
    static constexpr int stencil_width = 2;

    static constexpr int nx = (AMREX_SPACEDIM >= 1) ? stencil_width - 1 : 0;
    static constexpr int ny = (AMREX_SPACEDIM >= 2) ? stencil_width - 1 : 0;
    static constexpr int nz = (AMREX_SPACEDIM >= 3) ? stencil_width - 1 : 0;

    amrex::Real weights[3*stencil_width];

    template <typename P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Linear (const P& p,
            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
    {
        w = &weights[0];
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            amrex::Real l = (p.pos(i) - plo[i]) * dxi[i] + 0.5;
            index[i] = static_cast<int>(amrex::Math::floor(l)) - 1;
            amrex::Real lint = l - (index[i] + 1);
            w[stencil_width*i + 0] = 1.-lint;
            w[stencil_width*i + 1] = lint;
        }
        for (int i = AMREX_SPACEDIM; i < 3; ++i) {
            index[i] = 0;
            w[stencil_width*i + 0] = 1.;
            w[stencil_width*i + 1] = 0.;
        }
    }
};
}

#endif // include guard
