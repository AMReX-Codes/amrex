#ifndef ELEMENTARYPROCESS_H_
#define ELEMENTARYPROCESS_H_

#include "WarpXParticleContainer.H"
#include "CopyParticle.H"
#include "TransformParticle.H"

/**
 * \brief Base class for particle creation processes
 *
 * Particles in a source particle container are copied to a product particle
 * container if they are flagged. Both source and product particles can be
 * modified.
 *
 * initialize_copy initializes a general copy functor
 * createParticles formats the data to prepare for the copy, then
 *      calls copyAndTransformParticles
 * copyAndTransformParticles loops over particles, performs the copy and
 * transform source and product particles if needed.
 *
 * The class is templated on the process type. This gives flexibility
 * on source and product particle transformation.
 */
template<elementaryProcessType ProcessT>
class elementaryProcess
{
public:

    /**
     * \brief initialize and return functor for particle copy
     *
     * \param cpuid id of MPI rank
     * \param do_boosted_product whether to copy old attribs
     * \param runtime_uold_source momentum of source particles
     * \param attribs_source attribs of source particles
     * \param attribs_product attribs of product particles
     * \param runtime_attribs_product runtime attribs for product, to store old attribs
     */
    copyParticle initialize_copy(
        int cpuid, int do_boosted_product,
        const amrex::GpuArray<amrex::ParticleReal*,3> runtime_uold_source,
        const amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_source,
        const amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_product,
        const amrex::GpuArray<amrex::ParticleReal*,6> runtime_attribs_product) const noexcept
    {
        return copyParticle (
            cpuid, do_boosted_product, runtime_uold_source,
            attribs_source, attribs_product, runtime_attribs_product);
    };

    /**
     * \brief create particles in product particle containers
     *
     * For particle i in mfi, if is_flagged[i]=1, copy particle
     * particle i from container pc_source into pc_product
     *
     * \param lev MR level
     * \param mfi MFIter where source particles are located
     * \param pc_source source particle container
     * \param pc_product product particle container
     * \param is_flagged particles for which is_flagged is 1 are copied
     * \param do_boosted_product whether to copy old attribs
     */
    void createParticles (
        int lev, const amrex::MFIter& mfi,
        std::unique_ptr< WarpXParticleContainer>& pc_source,
        std::unique_ptr< WarpXParticleContainer>& pc_product,
        amrex::Gpu::ManagedDeviceVector<int>& is_flagged,
        bool do_boosted_product)
    {
        BL_PROFILE("createIonizedParticles");

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();

        // Get source particle data
        auto& ptile_source = pc_source->GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        const int np_source = ptile_source.GetArrayOfStructs().size();
        if (np_source == 0) return;
        // --- source AoS particle data
        WarpXParticleContainer::ParticleType* particles_source = ptile_source.GetArrayOfStructs()().data();
        // --- source SoA particle data
        auto& soa_source = ptile_source.GetStructOfArrays();
        amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_source;
        for (int ia = 0; ia < PIdx::nattribs; ++ia) {
            attribs_source[ia] = soa_source.GetRealData(ia).data();
        }
        // --- source runtime attribs
        amrex::GpuArray<amrex::ParticleReal*,3> runtime_uold_source;
        // Prepare arrays for boosted frame diagnostics.
        runtime_uold_source[0] = soa_source.GetRealData(PIdx::ux).data();
        runtime_uold_source[1] = soa_source.GetRealData(PIdx::uy).data();
        runtime_uold_source[2] = soa_source.GetRealData(PIdx::uz).data();

        // --- source runtime integer attribs
        amrex::GpuArray<int*,1> runtime_iattribs_source;
        std::map<std::string, int> icomps_source = pc_source->getParticleiComps();
        runtime_iattribs_source[0] = soa_source.GetIntData(icomps_source["ionization_level"]).data();

        // Indices of product particle for each ionized source particle.
        // i_product[i]-1 is the location in product tile of product particle
        // from source particle i.
        amrex::Gpu::ManagedDeviceVector<int> i_product;
        i_product.resize(np_source);
        // 0<i<np_source
        // 0<i_product<np_ionized
        // Strictly speaking, i_product should be an exclusive_scan of
        // is_flagged. However, for indices where is_flagged is 1, the
        // inclusive scan gives the same result with an offset of 1.
        // The advantage of inclusive_scan is that the sum of is_flagged
        // is in the last element, so no other reduction is required to get
        // number of particles.
        // Gpu::inclusive_scan runs on the current GPU stream, and synchronizes
        // with the CPU, so that the next line (executed by the CPU) has the
        // updated values of i_product
        amrex::Gpu::inclusive_scan(is_flagged.begin(), is_flagged.end(), i_product.begin());
        int np_ionized = i_product[np_source-1];
        if (np_ionized == 0) return;

        // Get product particle data
        auto& ptile_product = pc_product->GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        // old and new (i.e., including ionized particles) number of particles
        // for product species
        const int np_product_old = ptile_product.GetArrayOfStructs().size();
        const int np_product_new = np_product_old + np_ionized;
        // Allocate extra space in product species for ionized particles.
        ptile_product.resize(np_product_new);
        // --- product AoS particle data
        // First element is the first newly-created product particle
        WarpXParticleContainer::ParticleType* particles_product = ptile_product.GetArrayOfStructs()().data() + np_product_old;
        // --- product SoA particle data
        auto& soa_product = ptile_product.GetStructOfArrays();
        amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_product;
        for (int ia = 0; ia < PIdx::nattribs; ++ia) {
            // First element is the first newly-created product particle
            attribs_product[ia] = soa_product.GetRealData(ia).data() + np_product_old;
        }
        // --- product runtime attribs
        amrex::GpuArray<amrex::ParticleReal*,6> runtime_attribs_product;
        if (do_boosted_product) {
            std::map<std::string, int> comps_product = pc_product->getParticleComps();
            runtime_attribs_product[0] = soa_product.GetRealData(comps_product[ "xold"]).data() + np_product_old;
            runtime_attribs_product[1] = soa_product.GetRealData(comps_product[ "yold"]).data() + np_product_old;
            runtime_attribs_product[2] = soa_product.GetRealData(comps_product[ "zold"]).data() + np_product_old;
            runtime_attribs_product[3] = soa_product.GetRealData(comps_product["uxold"]).data() + np_product_old;
            runtime_attribs_product[4] = soa_product.GetRealData(comps_product["uyold"]).data() + np_product_old;
            runtime_attribs_product[5] = soa_product.GetRealData(comps_product["uzold"]).data() + np_product_old;
        }

        // --- product runtime integer attribs
        amrex::GpuArray<int*,1> runtime_iattribs_product;
        std::map<std::string, int> icomps_product = pc_product->getParticleiComps();
        runtime_iattribs_product[0] = soa_product.GetIntData(icomps_product["ionization_level"]).data() + np_product_old;

        int pid_product;
#pragma omp critical (doFieldIonization_nextid)
        {
            // ID of first newly-created product particle
            pid_product = pc_product->NextID();
            // Update NextID to include particles created in this function
            pc_product->setNextID(pid_product+np_ionized);
        }
        const int cpuid = amrex::ParallelDescriptor::MyProc();

        // Create instance of copy functor
        copyParticle copy_functor = initialize_copy(
            cpuid, do_boosted_product,
            runtime_uold_source,
            attribs_source,
            attribs_product,
            runtime_attribs_product);

        // Loop over source particles and create product particles
        copyAndTransformParticles(is_flagged, i_product, np_source, pid_product,
                                  particles_product, particles_source, copy_functor,
                                  runtime_iattribs_source, runtime_iattribs_product);
    }

    /**
     * \brief Loop over source particles and create product particles
     *
     * \param is_flagged particles for which is_flagged is 1 are copied
     */
    void copyAndTransformParticles(
        amrex::Gpu::ManagedDeviceVector<int>& is_flagged,
        amrex::Gpu::ManagedDeviceVector<int>& i_product,
        int np_source, int pid_product,
        WarpXParticleContainer::ParticleType* particles_product,
        WarpXParticleContainer::ParticleType* particles_source,
        copyParticle copy_functor,
        amrex::GpuArray<int*,1> runtime_iattribs_source,
        amrex::GpuArray<int*,1> runtime_iattribs_product)
    {
        int const * const AMREX_RESTRICT p_is_flagged = is_flagged.dataPtr();
        int const * const AMREX_RESTRICT p_i_product = i_product.dataPtr();
        // Loop over all source particles. If is_flagged, copy particle data
        // to corresponding product particle.
        amrex::For(
            np_source, [=] AMREX_GPU_DEVICE (int is) noexcept
            {
                if(p_is_flagged[is]){
                    // offset of 1 due to inclusive scan
                    int ip = p_i_product[is]-1;
                    // is: index of ionized particle in source species
                    // ip: index of corresponding new particle in product species
                    WarpXParticleContainer::ParticleType& p_product = particles_product[ip];
                    WarpXParticleContainer::ParticleType& p_source  = particles_source[is];
                    copy_functor(is, ip, pid_product, p_source, p_product);
                    transformSourceParticle<ProcessT>(is, p_source, runtime_iattribs_source);
                    transformProductParticle<ProcessT>(ip, p_product, runtime_iattribs_product);
                }
            }
            );
    }
};

// Derived class for ionization process
class IonizationProcess: public elementaryProcess<elementaryProcessType::Ionization>
{};

#endif // ELEMENTARYPROCESS_H_
