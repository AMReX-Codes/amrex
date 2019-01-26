
#ifndef _AGNParticleContainer_H_
#define _AGNParticleContainer_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include "NyxParticleContainer.H"

class AGNParticleContainer
    : public NyxParticleContainer<3+BL_SPACEDIM, 0, 0, 0>
{
public:
    
    using MyParIter = amrex::ParIter<3+BL_SPACEDIM>;

    AGNParticleContainer (amrex::Amr* amr, int nghost)
        : NyxParticleContainer<3+BL_SPACEDIM>(amr, nghost)
    {
        real_comp_names.clear();
        real_comp_names.push_back("mass");
        real_comp_names.push_back("xvel");
        real_comp_names.push_back("yvel");
        real_comp_names.push_back("zvel");
        real_comp_names.push_back("energy");
        real_comp_names.push_back("mdot");
    }
    
    virtual ~AGNParticleContainer () {}

    const int NumberOfParticles(MyParIter& pti) { return pti.GetArrayOfStructs().size(); }

    virtual void moveKickDrift (amrex::MultiFab& acceleration,
                                int level,
                                amrex::Real timestep,
                                amrex::Real a_old = 1.0,
                                amrex::Real a_half = 1.0,
                                int where_width = 0);

    virtual void moveKick      (amrex::MultiFab& acceleration,
                                int level,
                                amrex::Real timestep,
                                amrex::Real a_new = 1.0,
                                amrex::Real a_half = 1.0);

    void AddOneParticle (int lev,
                         int grid,
                         int tile,
                         amrex::Real mass, 
                         amrex::Real x,
                         amrex::Real y,
                         amrex::Real z)
    {
        auto& particle_tile = this->GetParticles(lev)[std::make_pair(grid,tile)];
        AddOneParticle(particle_tile, mass, x, y, z);
    }

    void AddOneParticle (ParticleTileType& particle_tile,
                         amrex::Real mass, 
                         amrex::Real x,
                         amrex::Real y,
                         amrex::Real z)
    {
        ParticleType p;
        p.id()  = ParticleType::NextID();
        p.cpu() = amrex::ParallelDescriptor::MyProc();
        p.pos(0) = x;
        p.pos(1) = y;
        p.pos(2) = z;

        // Set mass 
        p.rdata(0) = mass;

        // Zero initial velocity
        p.rdata(1) = 0.;
        p.rdata(2) = 0.;
        p.rdata(3) = 0.;

        // Zero initial energy
        p.rdata(4) = 0.;

        // Zero initial mdot
        p.rdata(5) = 0.;

        particle_tile.push_back(p);
    }
    ///
    /// Invalidate particles in cells that are already occupied
    ///
    void ComputeOverlap(int lev);

    ///
    /// Invalidate particles that have been merged with other particles
    ///
    void Merge(int lev);

    ///
    /// Compute the momentum that has been removed from the gas in order to define the particle velocity
    ///
    void ComputeParticleVelocity(int lev,
                                 amrex::MultiFab& state_old,
                                 amrex::MultiFab& state_new,
                                 int add_energy);

    ///
    /// Accrete mass from the grid onto the existing AGN particles
    ///
    void AccreteMass(int lev,
                     amrex::MultiFab& state,
                     amrex::MultiFab& density_lost,
                     amrex::Real dt);

    ///
    ///  Release energy if it exceeds thermal feedback threshold.
    ///
    void ReleaseEnergy(int lev,
                       amrex::MultiFab& state,
                       amrex::MultiFab& D_new,
                       amrex::Real a);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel(int lev);

 protected:
    
    bool sub_cycle;
    amrex::Vector<std::string> real_comp_names;

};

#endif /* _AGNParticleContainer_H_ */
