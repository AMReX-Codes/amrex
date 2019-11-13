#ifndef WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_
#define WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_

#include "WarpXParticleContainer.H"
//#include <WarpX.H>
#include <AMReX_REAL.H>
#include <AMReX_DenseBins.H>

class CollisionType
{
public:
    int m_species1;
    int m_species2;

    CollisionType(
        const std::vector<std::string>& species_names,
        std::string collision_name);

    static void ElasticCollisionPerez(
        amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const I1s,
        amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const I1e,
        amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const I2s,
        amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const I2e,
        amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type *I1,
        amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type *I2,
        amrex::Real* u1x, amrex::Real* u1y, amrex::Real* u1z,
        amrex::Real* u2x, amrex::Real* u2y, amrex::Real* u2z,
        const amrex::Real* w1, const amrex::Real* w2,
        const amrex::Real  q1, const amrex::Real  q2,
        const amrex::Real  m1, const amrex::Real  m2,
        const amrex::Real  T1, const amrex::Real  T2,
        const amrex::Real  dt, const amrex::Real   L, const amrex::Real V);

    static void ShuffleFisherYates(
        amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type *array,
        amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const is,
        amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const ie);

};

#endif // WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_
