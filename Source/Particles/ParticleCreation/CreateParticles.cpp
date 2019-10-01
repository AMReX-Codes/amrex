#include "IonizationEvent.H"
#include "MultiParticleContainer.H"
#include "WarpX.H"

using namespace amrex;

// For particle i in mfi, if is_flagged[i]=1, copy particle
// particle i from container pc_source into pc_product
void MultiParticleContainer::createParticles (
    int lev, const MFIter& mfi,
    std::unique_ptr< WarpXParticleContainer>& pc_source,
    std::unique_ptr< WarpXParticleContainer>& pc_product,
    amrex::Gpu::ManagedDeviceVector<int>& is_flagged,
    elementaryProcess elementary_process)
{
    BL_PROFILE("createIonizedParticles");

    const int * const AMREX_RESTRICT p_is_flagged = is_flagged.dataPtr();

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
    GpuArray<ParticleReal*,PIdx::nattribs> attribs_source;
    for (int ia = 0; ia < PIdx::nattribs; ++ia) {
        attribs_source[ia] = soa_source.GetRealData(ia).data();
    }
    // --- source runtime attribs
    GpuArray<ParticleReal*,3> runtime_uold_source;
    // Prepare arrays for boosted frame diagnostics.
    runtime_uold_source[0] = soa_source.GetRealData(PIdx::ux).data();
    runtime_uold_source[1] = soa_source.GetRealData(PIdx::uy).data();
    runtime_uold_source[2] = soa_source.GetRealData(PIdx::uz).data();

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
    int* AMREX_RESTRICT p_i_product = i_product.dataPtr();

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
    GpuArray<ParticleReal*,PIdx::nattribs> attribs_product;
    for (int ia = 0; ia < PIdx::nattribs; ++ia) {
        // First element is the first newly-created product particle
        attribs_product[ia] = soa_product.GetRealData(ia).data() + np_product_old;
    }
    // --- product runtime attribs
    GpuArray<ParticleReal*,6> runtime_attribs_product;
    bool do_boosted_product = WarpX::do_boosted_frame_diagnostic
        && pc_product->DoBoostedFrameDiags();
    if (do_boosted_product) {
        std::map<std::string, int> comps_product = pc_product->getParticleComps();
        runtime_attribs_product[0] = soa_product.GetRealData(comps_product[ "xold"]).data() + np_product_old;
        runtime_attribs_product[1] = soa_product.GetRealData(comps_product[ "yold"]).data() + np_product_old;
        runtime_attribs_product[2] = soa_product.GetRealData(comps_product[ "zold"]).data() + np_product_old;
        runtime_attribs_product[3] = soa_product.GetRealData(comps_product["uxold"]).data() + np_product_old;
        runtime_attribs_product[4] = soa_product.GetRealData(comps_product["uyold"]).data() + np_product_old;
        runtime_attribs_product[5] = soa_product.GetRealData(comps_product["uzold"]).data() + np_product_old;
    }

    int pid_product;
#pragma omp critical (doFieldIonization_nextid)
    {
        // ID of first newly-created product particle
        pid_product = pc_product->NextID();
        // Update NextID to include particles created in this function
        pc_product->setNextID(pid_product+np_ionized);
    }
    const int cpuid = ParallelDescriptor::MyProc();

    copyAndTransformParticle copy_and_transform_functor = elementary_process.initialize_functor(
        cpuid, do_boosted_product,
        runtime_uold_source,
        attribs_source,
        attribs_product,
        runtime_attribs_product);
    
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

                copy_and_transform_functor(
                    is, ip, pid_product,
                    p_source,
                    p_product);
            }
        }
        );
}
