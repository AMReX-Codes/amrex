#ifndef TRANSFORMPARTICLE_H_
#define TRANSFORMPARTICLE_H_

#include "WarpXParticleContainer.H"

enum struct elementaryProcessType { Ionization };

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformSourceParticle(
    int i,
    WarpXParticleContainer::ParticleType& particle)
{
    amrex::Print()<<"primary template\n";
}

/*
template < elementaryProcessType ProcessT >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformProductParticle(
    int i,
    WarpXParticleContainer::ParticleType& particle)
{
    amrex::Print()<<"primary template\n";
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformSourceParticle < elementaryProcessType:Ionization > (
    int i,
    WarpXParticleContainer::ParticleType& particle)
{
    amrex::Print()<<"specialized source\n";
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformProductParticle < elementaryProcessType:Ionization> (
    int i,
    WarpXParticleContainer::ParticleType& particle)
{
    amrex::Print()<<"specialized source\n";
}
*/









/*
enum struct PType {source, product};

template < elementaryProcess ProcessT, PType PT >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformParticle(
    int i,
    WarpXParticleContainer::ParticleType& particle)
{
    amrex::Print()<<"primary template\n";
}

template <elementaryProcess ProcessT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformParticle < ProcessT, PType::source > (
    int i,
    WarpXParticleContainer::ParticleType& particle)
{
    amrex::Print()<<"specialized source\n";
}

template <elementaryProcess ProcessT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformParticle < ProcessT, PType::product > (
    int i,
    WarpXParticleContainer::ParticleType& particle)
{
    amrex::Print()<<"specialized product\n";
}
*/
#endif // TRANSFORMPARTICLE_H_
