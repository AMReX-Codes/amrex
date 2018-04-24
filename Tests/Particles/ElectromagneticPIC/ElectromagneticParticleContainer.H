#ifndef ELECTROMAGNETIC_PARTICLE_CONTAINER_H_
#define ELECTROMAGNETIC_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

struct PIdx
{
    enum { // Particle Attributes stored in amrex::ParticleContainer's struct of arrays
        w = 0,  // weight
        ux, uy, uz, Ex, Ey, Ez, Bx, By, Bz, ginv,
        nattribs
    };
};

class ElectromagneticParIter
    : public amrex::ParIter<0,0,PIdx::nattribs>
{
public:
    using amrex::ParIter<0,0,PIdx::nattribs>::ParIter;

    const std::array<amrex::Vector<amrex::Real>, PIdx::nattribs>& GetAttribs () const {
        return GetStructOfArrays().GetRealData();
    }
    
    std::array<amrex::Vector<amrex::Real>, PIdx::nattribs>& GetAttribs () {
        return GetStructOfArrays().GetRealData();
    }
    
    const amrex::Vector<amrex::Real>& GetAttribs (int comp) const {
        return GetStructOfArrays().GetRealData(comp);
    }
    
    amrex::Vector<amrex::Real>& GetAttribs (int comp) {
        return GetStructOfArrays().GetRealData(comp);
    }
};

class ElectromagneticParticleContainer
    : public amrex::ParticleContainer<0,0,PIdx::nattribs>
{
public:

    ElectromagneticParticleContainer (const amrex::Geometry            & a_geom,
                                      const amrex::DistributionMapping & a_dmap,
                                      const amrex::BoxArray            & a_ba, 
                                      const int                         a_species_id,
                                      const amrex::Real                 a_charge,
                                      const amrex::Real                 a_mass);

    void InitParticles(const amrex::IntVect& a_num_particles_per_cell,
                       const amrex::Real     a_thermal_momentum_std, 
                       const amrex::Real     a_thermal_momentum_mean,
                       const amrex::Real     a_density);

    void PushAndDeposeParticles(const amrex::MultiFab& Ex,
                                const amrex::MultiFab& Ey,
                                const amrex::MultiFab& Ez,
                                const amrex::MultiFab& Bx,
                                const amrex::MultiFab& By,
                                const amrex::MultiFab& Bz,
                                      amrex::MultiFab& jx, 
                                      amrex::MultiFab& jy, 
                                      amrex::MultiFab& jz,
                                      amrex::Real      dt);
    
    void PushParticlesOnly(amrex::Real dt);

    void EnforcePeriodicBCs();

protected:
    
    int m_species_id;
    
    amrex::Real m_charge;
    amrex::Real m_mass;
};

#endif
