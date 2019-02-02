#ifndef MD_PARTICLE_CONTAINER_H_
#define MD_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

struct PIdx
{
    enum {
        vx = 0,
        vy, vz, ax, ay, az,
        nattribs
    };
};

class MDParIter
    : public amrex::ParIter<0,0,PIdx::nattribs,0>
{
public:
    using amrex::ParIter<0,0,PIdx::nattribs,0>::ParIter;

    const std::array<RealVector, PIdx::nattribs>& GetAttribs () const { 
        return GetStructOfArrays().GetRealData(); 
    }

    std::array<RealVector, PIdx::nattribs>& GetAttribs () { 
        return GetStructOfArrays().GetRealData(); 
    }

    const RealVector& GetAttribs (int comp) const { 
        return GetStructOfArrays().GetRealData(comp);
    }

    RealVector& GetAttribs (int comp) { 
        return GetStructOfArrays().GetRealData(comp);
    }
};

class MDParticleContainer
    : public amrex::ParticleContainer<0, 0, PIdx::nattribs, 0>
{

public:

    MDParticleContainer (const amrex::Geometry            & a_geom,
                         const amrex::DistributionMapping & a_dmap,
                         const amrex::BoxArray            & a_ba);

    void InitParticles(const amrex::IntVect& a_num_particles_per_cell,
                       const amrex::Real     a_thermal_momentum_std,
                       const amrex::Real     a_thermal_momentum_mean);

    void RedistributeLocal()
    {
        const int lev_min = 0;
        const int lev_max = 0;
        const int nGrow = 0;
        const int local = 1;
        Redistribute(lev_min, lev_max, nGrow, local);
    }

    void BuildNeighborList();
};

#endif
