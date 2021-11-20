#ifndef AMREX_PARTICLEINTERPOLATORS_H_
#define AMREX_PARTICLEINTERPOLATORS_H_
#include <AMReX_Config.H>

#include <AMReX_IntVect.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>

namespace amrex
{
namespace ParticleInterpolator
{

template <class Derived, class WeightType>
struct Base
{
    int index[3];
    WeightType* w;

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
}

#endif // include guard
