#ifndef _NEIGHBORLISTPARTICLECONTAINER_H_
#define _NEIGHBORLISTPARTICLECONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_NeighborParticles.H"

namespace amrex
{

struct Tag
{
    int lev;
    int grid;
    int tile;
    
    Tag (int a_lev, int a_grid, int a_tile)
        : lev(a_lev), grid(a_grid), tile(a_tile)
        {}
    
    bool operator==(const Tag& other) const
    {
        return ((lev == other.lev) and (grid == other.grid) and (tile == other.tile));
    }
    
    bool operator!=(const Tag& other) const
    {
        return !operator==(other);
    }
    
    friend std::ostream& operator<< (std::ostream& stream, const Tag& tag);        
};
    
class NeighborListParticleContainer
    : public NeighborParticleContainer<2*BL_SPACEDIM, 0>
{
public:
    
    ///
    /// This particle container fills a mask for quickly computing
    /// neighbor grids / tiles for a given particle
    ///
    NeighborListParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba,
                                  int                         nneighbor);

    NeighborListParticleContainer(const Vector<Geometry>            & geom, 
                                  const Vector<DistributionMapping> & dmap,
                                  const Vector<BoxArray>            & ba,
                                  const Vector<int>                 & rr,
                                  int                               nneighbor);
    
    ///
    /// Init one particle per cell with random velocities
    ///
    void InitParticles();

    ///
    /// Compute the short range forces on a tile's worth of particles.
    /// fillNeighbors must have already been called.
    ///
    void computeForces();

    ///
    /// Compute the short range forces on a tile's worth of particles using
    /// the neighbor list instead of the N^2 approach.
    /// fillNeighbors must have already been called.
    ///
    void computeForcesNL();

    ///
    /// Move the particles according to their forces, reflecting at domain boundaries
    ///
    void moveParticles(const Real dt);

    ///
    /// Save the particle data in an ASCII format
    ///
    void writeParticles(int n);

    void getNeighborParticles(const int nGrow);
    
private:

    Vector<Tag> getNeighborGrids(const ParticleType& p,
                                 const Tag src,
                                 const int nGrow);
    
    Vector<Tag> getNeighborGrids(const ParticleType& p,
                                 const Tag src,
                                 const IntVect& nGrow);
    
    struct
    {
        bool operator()(const ParticleType& p1, const ParticleType& p2) const
        {
            return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)) ,
                                   + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
                                   + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) )
                <= cutoff*cutoff;
        }
    } CheckPair;
    
    static constexpr Real cutoff = 1.e-2;
    static constexpr Real min_r  = 1.e-4;
};
    
}

#endif // _NEIGHBORLISTPARTICLECONTAINER_H_
