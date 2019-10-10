#ifndef TRANSFORMPARTICLE_H_
#define TRANSFORMPARTICLE_H_

#include "WarpXParticleContainer.H"

enum struct PType {source, product};

template < PType PT >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformParticle(
    int i,
    WarpXParticleContainer::ParticleType& particle)
{
    amrex::Print()<<"primary template\n";
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformParticle < PType::source > (
    int i,
    WarpXParticleContainer::ParticleType& particle)
{
    amrex::Print()<<"specialized source\n";
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformParticle < PType::product > (
    int i,
    WarpXParticleContainer::ParticleType& particle)
{
    amrex::Print()<<"specialized product\n";
}

#endif // TRANSFORMPARTICLE_H_
