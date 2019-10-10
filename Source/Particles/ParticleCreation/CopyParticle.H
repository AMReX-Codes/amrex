#ifndef COPYPARTICLE_H_
#define COPYPARTICLE_H_

#include "WarpXParticleContainer.H"

class copyParticle
{

public:

    int m_cpuid;
    bool m_do_boosted_product;
    amrex::GpuArray<amrex::ParticleReal*,3> m_runtime_uold_source;
    amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> m_attribs_source;
    amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> m_attribs_product;
    amrex::GpuArray<amrex::ParticleReal*,6> m_runtime_attribs_product;

    AMREX_GPU_HOST_DEVICE
    copyParticle(
        const int cpuid, const int do_boosted_product,
        const amrex::GpuArray<amrex::ParticleReal*,3> runtime_uold_source,
        const amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_source,
        const amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_product,
        const amrex::GpuArray<amrex::ParticleReal*,6> runtime_attribs_product)
        :
        m_cpuid(cpuid),
        m_do_boosted_product(do_boosted_product),
        m_runtime_uold_source(runtime_uold_source),
        m_attribs_source(attribs_source),
        m_attribs_product(attribs_product),
        m_runtime_attribs_product(runtime_attribs_product){}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    //void operator () (int is, int ip, int pid_product,
    //                  WarpXParticleContainer::ParticleType& p_source,
    //                  WarpXParticleContainer::ParticleType& p_product) const noexcept
    void operator () (int is, int ip, int pid_product,
                      WarpXParticleContainer::ParticleType& p_source,
                      WarpXParticleContainer::ParticleType& p_product)
    {
        // Copy particle from source to product: AoS
        p_product.id() = pid_product + ip;
        p_product.cpu() = m_cpuid;
        p_product.pos(0) = p_source.pos(0);
        p_product.pos(1) = p_source.pos(1);
#if (AMREX_SPACEDIM == 3)
        p_product.pos(2) = p_source.pos(2);
#endif
        // Copy particle from source to product: SoA
        for (int ia = 0; ia < PIdx::nattribs; ++ia) {
            m_attribs_product[ia][ip] = m_attribs_source[ia][is];
        }
        // Update xold etc. if boosted frame diagnostics required
        // for product species. Fill runtime attribs with a copy of
        // current properties (xold = x etc.).
        if (m_do_boosted_product) {
            m_runtime_attribs_product[0][ip] = p_source.pos(0);
            m_runtime_attribs_product[1][ip] = p_source.pos(1);
            m_runtime_attribs_product[2][ip] = p_source.pos(2);
            m_runtime_attribs_product[3][ip] = m_runtime_uold_source[0][ip];
            m_runtime_attribs_product[4][ip] = m_runtime_uold_source[1][ip];
            m_runtime_attribs_product[5][ip] = m_runtime_uold_source[2][ip];
        }
    }
};

#endif // COPYPARTICLE_H_
