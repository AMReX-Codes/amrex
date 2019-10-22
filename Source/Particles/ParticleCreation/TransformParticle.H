#ifndef TRANSFORMPARTICLE_H_
#define TRANSFORMPARTICLE_H_

#include "WarpXParticleContainer.H"

enum struct elementaryProcessType { Ionization };

/**
 * \brief small modifications on source particle
 * \param i index of particle
 * \param particle particle struct
 * \param runtime_iattribs integer attribs
 */
template < elementaryProcessType ProcessT >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformSourceParticle(
    int i,
    WarpXParticleContainer::ParticleType& particle,
    amrex::GpuArray<int*,1> runtime_iattribs){}

/**
 * \brief small modifications on product particle
 * \param i index of particle
 * \param particle particle struct
 */
template < elementaryProcessType ProcessT >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformProductParticle(
    int i,
    WarpXParticleContainer::ParticleType* * v_particle){}

// For ionization, increase ionization level of source particle
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformSourceParticle < elementaryProcessType::Ionization > (
    int i,
    WarpXParticleContainer::ParticleType& particle,
    amrex::GpuArray<int*,1> runtime_iattribs)
{
    // increment particle's ionization level
    runtime_iattribs[0][i] += 1;
}

#endif // TRANSFORMPARTICLE_H_
