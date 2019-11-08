#ifndef ELEMENTARYPROCESS_H_
#define ELEMENTARYPROCESS_H_

#include "WarpXParticleContainer.H"
#include "CopyParticle.H"
#include "TransformParticle.H"

/**
 * \brief Base class for particle creation processes
 *
 * Particles in a source particle container are copied to product particle
 * containers if they are flagged. Both source and product particles can be
 * modified.
 *
 * - initialize_copy initializes a general copy functor
 * - createParticles formats the data to prepare for the copy, then
 *       calls copyAndTransformParticles
 * - copyAndTransformParticles loops over particles, performs the copy and
 *       transform source and product particles if needed.
 *
 * The class is templated on the process type. This gives flexibility
 * on source and product particle transformations.
 */
template<ElementaryProcessType ProcessT>
class ElementaryProcess
{
public:

    /**
     * \brief initialize and return functor for particle copy
     *
     * \param cpuid id of MPI rank
     * \param do_back_transformed_product whether to copy old attribs
     * \param runtime_uold_source momentum of source particles
     * \param attribs_source attribs of source particles
     * \param attribs_product attribs of product particles
     * \param runtime_attribs_product runtime attribs for product, to store old attribs
     */
    CopyParticle initialize_copy(
        const int cpuid, const int do_back_transformed_product,
        const amrex::GpuArray<amrex::ParticleReal*,3> runtime_uold_source,
        const amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_source,
        const amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_product,
        const amrex::GpuArray<amrex::ParticleReal*,6> runtime_attribs_product) const noexcept
    {
        return CopyParticle(
            cpuid, do_back_transformed_product, runtime_uold_source,
            attribs_source, attribs_product, runtime_attribs_product
        );
    };

    /**
     * \brief create particles in product particle containers
     *
     * For particle i in mfi, if is_flagged[i]=1, copy particle
     * particle i from container pc_source into each container in v_pc_product
     *
     * \param lev MR level
     * \param mfi MFIter where source particles are located
     * \param pc_source source particle container
     * \param v_pc_product vector of product particle containers
     * \param is_flagged particles for which is_flagged is 1 are copied
     * \param v_do_back_transformed_product vector: whether to copy old attribs for
     *        each product container
     */
    void createParticles (
        int lev, const amrex::MFIter& mfi,
        std::unique_ptr< WarpXParticleContainer>& pc_source,
        amrex::Vector<WarpXParticleContainer*> v_pc_product,
        amrex::Gpu::ManagedDeviceVector<int>& is_flagged,
        amrex::Gpu::ManagedDeviceVector<int> v_do_back_transformed_product)
    {

        BL_PROFILE("createParticles");
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

        int nproducts = v_pc_product.size();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            v_do_back_transformed_product.size() == nproducts,
            "v_pc_product and v_do_back_transformed_product must have the same size.");

        // Indices of product particle for each source particle.
        // i_product[i]-1 is the location in product tile of product particle
        // from source particle i.
        amrex::Gpu::ManagedDeviceVector<int> i_product;
        i_product.resize(np_source);
        // 0<i<np_source
        // 0<i_product<np_flagged
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
        int np_flagged = i_product[np_source-1];
        if (np_flagged == 0) return;

        amrex::Gpu::ManagedDeviceVector<CopyParticle> v_copy_functor;
        v_copy_functor.reserve(nproducts);
        amrex::Gpu::ManagedDeviceVector<int> v_pid_product(nproducts);
        amrex::Gpu::ManagedDeviceVector<WarpXParticleContainer::ParticleType*> v_particles_product(nproducts);
        for (int iproduct=0; iproduct<nproducts; iproduct++){
            WarpXParticleContainer*& pc_product = v_pc_product[iproduct];
            // Get product particle data
            auto& ptile_product = pc_product->GetParticles(lev)[std::make_pair(grid_id,tile_id)];
            // old and new (i.e., including flagged particles) number of particles
            // for product species
            const int np_product_old = ptile_product.GetArrayOfStructs().size();
            const int np_product_new = np_product_old + np_flagged;
            // Allocate extra space in product species for flagged particles.
            ptile_product.resize(np_product_new);
            // --- product AoS particle data
            // WarpXParticleContainer::ParticleType* particles_product = ptile_product.GetArrayOfStructs()().data() + np_product_old;
            v_particles_product[iproduct] = ptile_product.GetArrayOfStructs()().data() + np_product_old;
            // First element is the first newly-created product particle
            // --- product SoA particle data
            auto& soa_product = ptile_product.GetStructOfArrays();
            amrex::GpuArray<amrex::ParticleReal*,PIdx::nattribs> attribs_product;
            for (int ia = 0; ia < PIdx::nattribs; ++ia) {
                // First element is the first newly-created product particle
                attribs_product[ia] = soa_product.GetRealData(ia).data() + np_product_old;
            }
            // --- product runtime attribs
            amrex::GpuArray<amrex::ParticleReal*,6> runtime_attribs_product;
            const int do_boost = v_do_back_transformed_product[iproduct];
            if (do_boost) {
                std::map<std::string, int> comps_product = pc_product->getParticleComps();
                runtime_attribs_product[0] = soa_product.GetRealData(comps_product[ "xold"]).data() + np_product_old;
                runtime_attribs_product[1] = soa_product.GetRealData(comps_product[ "yold"]).data() + np_product_old;
                runtime_attribs_product[2] = soa_product.GetRealData(comps_product[ "zold"]).data() + np_product_old;
                runtime_attribs_product[3] = soa_product.GetRealData(comps_product["uxold"]).data() + np_product_old;
                runtime_attribs_product[4] = soa_product.GetRealData(comps_product["uyold"]).data() + np_product_old;
                runtime_attribs_product[5] = soa_product.GetRealData(comps_product["uzold"]).data() + np_product_old;
            }

            // --- product runtime integer attribs
            int pid_product;
#pragma omp critical (doParticleCreation_nextid)
            {
                // ID of first newly-created product particle
                pid_product = pc_product->NextID();
                // Update NextID to include particles created in this function
                pc_product->setNextID(pid_product+np_flagged);
            }
            const int cpuid = amrex::ParallelDescriptor::MyProc();

            // Create instance of copy functor, and add it to the vector
            v_copy_functor.push_back( initialize_copy(
                cpuid, v_do_back_transformed_product[iproduct],
                runtime_uold_source,
                attribs_source,
                attribs_product,
                runtime_attribs_product
            ) );
            v_pid_product[iproduct] = pid_product;

        }
        // Loop over source particles and create product particles
        copyAndTransformParticles(is_flagged, i_product, np_source, v_pid_product,
                                  v_particles_product, particles_source, v_copy_functor,
                                  runtime_iattribs_source);
        // Synchronize to prevent the destruction of temporary arrays (at the
        // end of the function call) before the kernel executes.
        amrex::Gpu::streamSynchronize();
    }

    /**
     * \brief Loop over source particles and create product particles
     *
     * \param is_flagged particles for which is_flagged is 1 are copied
     * \param i_product relative indices of product particle. This is the same
     *        for all product particle containers.
     * \param np_source number of particles in source container
     * \param v_pid_product for each product particle container: ID of the
     *        first product particle. Others IDs are incremented from this one
     * \param v_particles_product for each product particle container: pointer
     *        to the AoS particle data
     * \param particles_source pointter to the AoS source particle data
     * \param v_copy_functor vector of copy functors, one for each product particle container
     * \param runtime_iattribs_source pointer to source runtime integer attribs
     */
    void copyAndTransformParticles(
        amrex::Gpu::ManagedDeviceVector<int>& is_flagged,
        amrex::Gpu::ManagedDeviceVector<int>& i_product,
        int np_source,
        amrex::Gpu::ManagedDeviceVector<int> v_pid_product,
        amrex::Gpu::ManagedDeviceVector<WarpXParticleContainer::ParticleType*> v_particles_product,
        WarpXParticleContainer::ParticleType* particles_source,
        const amrex::Gpu::ManagedDeviceVector<CopyParticle>& v_copy_functor,
        amrex::GpuArray<int*,1> runtime_iattribs_source)
    {
        int const * const AMREX_RESTRICT p_is_flagged = is_flagged.dataPtr();
        int const * const AMREX_RESTRICT p_i_product = i_product.dataPtr();
        CopyParticle const * const AMREX_RESTRICT p_copy_functor = v_copy_functor.dataPtr();
        WarpXParticleContainer::ParticleType * * p_particles_product = v_particles_product.data();
        int const * const p_pid_product = v_pid_product.data();

        // Loop over all source particles. If is_flagged, copy particle data
        // to corresponding product particle.
        int nproducts = v_particles_product.size();
        amrex::For(
            np_source, [=] AMREX_GPU_DEVICE (int is) noexcept
            {
                if(p_is_flagged[is]){
                    // offset of 1 due to inclusive scan
                    int ip = p_i_product[is]-1;
                    WarpXParticleContainer::ParticleType& p_source  = particles_source[is];
                    for (int iproduct=0; iproduct<nproducts; iproduct++){
                        // is: index of flagged particle in source species
                        // ip: index of corresponding new particle in product species
                        WarpXParticleContainer::ParticleType* particles_product = p_particles_product[iproduct];
                        WarpXParticleContainer::ParticleType& p_product = particles_product[ip];
                        p_copy_functor[iproduct](is, ip, p_pid_product[iproduct], p_source, p_product);
                    }
                    // Transform source particle
                    transformSourceParticle<ProcessT>(is, p_source, runtime_iattribs_source);
                    // Transform product particles. The loop over product
                    // species is done inside the function to allow for
                    // more flexibility.
                    transformProductParticle<ProcessT>(ip, p_particles_product);
                }
            }
        );
    }
};

// Derived class for ionization process
class IonizationProcess: public ElementaryProcess<ElementaryProcessType::Ionization>
{};

#endif // ELEMENTARYPROCESS_H_
