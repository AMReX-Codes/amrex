#include "ElementaryProcess.H"

using namespace amrex;

void
elementaryProcess::createParticles (
    int lev, const MFIter& mfi,
    std::unique_ptr< WarpXParticleContainer>& pc_source,
    std::unique_ptr< WarpXParticleContainer>& pc_product,
    Gpu::ManagedDeviceVector<int>& is_flagged,
    bool do_boosted_product)

void
elementaryProcess::copyAndTransformParticles(
    Gpu::ManagedDeviceVector<int>& is_flagged,
    Gpu::ManagedDeviceVector<int>& i_product,
    int np_source, int pid_product,
    WarpXParticleContainer::ParticleType* particles_product,
    WarpXParticleContainer::ParticleType* particles_source,
    copyParticle copy_functor)
