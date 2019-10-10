#ifndef ELEMENTARYPROCESS_H_
#define ELEMENTARYPROCESS_H_

#include "WarpXParticleContainer.H"
#include "CopyParticle.H"
#include "TransformParticle.H"

/*
enum struct elementaryProcessType { Ionization };
enum struct PType {source, product};

  namespace transform
  {
  template < elementaryProcessType ElemT, PType PartT >
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void transformParticle(
  int i,
  WarpXParticleContainer::ParticleType& particle)
  {
  amrex::Print()<<"primary template\n";
  }

  template <>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void transformParticle<elementaryProcessType::Ionization, PType::source>(
  int i,
  WarpXParticleContainer::ParticleType& particle)
  {
  amrex::Print()<<"ionization source template\n";
  }

  template <>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void transformParticle<elementaryProcessType::Ionization, PType::product>(
  int i,
  WarpXParticleContainer::ParticleType& particle)
  {
  amrex::Print()<<"ionization product template\n";
  }
  }
*/

class elementaryProcess
{

public:

    elementaryProcess() {};

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

    // For particle i in mfi, if is_flagged[i]=1, copy particle
    // particle i from container pc_source into pc_product
    void createParticles (
        int lev, const amrex::MFIter& mfi,
        std::unique_ptr< WarpXParticleContainer>& pc_source,
        std::unique_ptr< WarpXParticleContainer>& pc_product,
        amrex::Gpu::ManagedDeviceVector<int>& is_flagged,
        bool do_boosted_product);

    void copyAndTransformParticles(
        amrex::Gpu::ManagedDeviceVector<int>& is_flagged,
        amrex::Gpu::ManagedDeviceVector<int>& i_product,
        int np_source, int pid_product,
        WarpXParticleContainer::ParticleType* particles_product,
        WarpXParticleContainer::ParticleType* particles_source,
        copyParticle copy_functor);
};

#endif // ELEMENTARYPROCESS_H_
