#ifndef COPYPARTICLE_H_
#define COPYPARTICLE_H_

#include "WarpXParticleContainer.H"

/**
 * \brief Functor to copy one particle
 *
 * This is meant to be a very small class captured by value in kernel launches,
 * that can be initialized on the host and copied to the device at each
 * iteration. It should be general enough to be used by all processes.
 *
 * Pointers to SoA data are saved when constructor is called.
 * AoS data is passed as argument of operator ().
 */
class CopyParticle
{

public:

    // ID of MPI rank
    int m_cpuid;
    // If true, will copy old attribs for back-transforme diagnostics
    bool m_do_back_transformed_product;
    // Source old (runtime) attribs for back-transformed diagnostics
    amrex::GpuArray<amrex::ParticleReal*,3> m_runtime_uold_source;
    // Source attribs
    amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> m_attribs_source;
    // Product attribs
    amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> m_attribs_product;
    // Product runtime attribs
    amrex::GpuArray<amrex::ParticleReal*,6> m_runtime_attribs_product;

    // Simple constructor
    AMREX_GPU_HOST_DEVICE
    CopyParticle(
        const int cpuid, const int do_back_transformed_product,
        const amrex::GpuArray<amrex::ParticleReal*,3> runtime_uold_source,
        const amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_source,
        const amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_product,
        const amrex::GpuArray<amrex::ParticleReal*,6> runtime_attribs_product)
        :
        m_cpuid(cpuid),
        m_do_back_transformed_product(do_back_transformed_product),
        m_runtime_uold_source(runtime_uold_source),
        m_attribs_source(attribs_source),
        m_attribs_product(attribs_product),
        m_runtime_attribs_product(runtime_attribs_product){}

    /**
     * \brief Overload operator () to do the copy of 1 particle
     *
     * \param is index of source particle
     * \param ip index of product particle
     * \param pid_product ID of first product particle
     * \param p_source Struct with data for 1 source particle (position etc.)
     * \param p_product Struct with data for 1 product particle (position etc.)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (int is, int ip, int pid_product,
                      WarpXParticleContainer::ParticleType& p_source,
                      WarpXParticleContainer::ParticleType& p_product) const noexcept
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
        if (m_do_back_transformed_product) {
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
