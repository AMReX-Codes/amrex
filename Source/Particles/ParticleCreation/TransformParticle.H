#ifndef TRANSFORMPARTICLE_H_
#define TRANSFORMPARTICLE_H_

#include "WarpXParticleContainer.H"

enum struct elementaryProcessType { Ionization };

template < elementaryProcessType ProcessT >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformSourceParticle(
    int i,
    WarpXParticleContainer::ParticleType& particle,
    amrex::GpuArray<int*,1> runtime_iattribs_source)
{
    amrex::Print()<<"primary template\n";
}

template < elementaryProcessType ProcessT >
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformProductParticle(
    int i,
    WarpXParticleContainer::ParticleType& particle,
    amrex::GpuArray<int*,1> runtime_iattribs_product)
{
    amrex::Print()<<"primary template\n";
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformSourceParticle < elementaryProcessType::Ionization > (
    int i,
    WarpXParticleContainer::ParticleType& particle,
    amrex::GpuArray<int*,1> runtime_iattribs_source)
{
    // increment particle's ionization level
    runtime_iattribs_source[0][i] += 1; 
    amrex::Print()<<"specialized source\n";
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void transformProductParticle < elementaryProcessType::Ionization > (
    int i,
    WarpXParticleContainer::ParticleType& particle,
    amrex::GpuArray<int*,1> runtime_iattribs_product)
{
    amrex::Print()<<"specialized source\n";
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
