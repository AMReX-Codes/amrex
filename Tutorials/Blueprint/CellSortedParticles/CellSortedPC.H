#ifndef CELL_SORTED_PARTICLE_CONTAINER_H_
#define CELL_SORTED_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

struct RealData
{
    enum {
        vx = 0,
        vy,
        vz,
        ncomps
    };
};

struct IntData
{
    enum {
        sorted,
        i,
	j,
	k,
        ncomps
    };
};

class CellSortedParticleContainer
    : public amrex::ParticleContainer<RealData::ncomps, IntData::ncomps>
{
public:

    using MyParIter = amrex::ParIter<RealData::ncomps, IntData::ncomps>;
    
    CellSortedParticleContainer (const amrex::Geometry            & a_geom,
                                 const amrex::DistributionMapping & a_dmap,
                                 const amrex::BoxArray            & a_ba);

    void InitParticles(const amrex::IntVect& a_num_particles_per_cell);

    void ReBin();

    void MoveParticles();

    int SumCellVectors();

    int numUnsorted();

    int numWrongCell();

    void visitAllParticles();
    
protected:

    void UpdateCellVectors();
    
    void UpdateFortranStructures();

    // used to store vectors of particle indices on a cell-by-cell basis
    std::map<int, amrex::BaseFab<std::vector<int> > > m_cell_vectors;

    // primitive data for passing the above into Fortran
    std::map<int, amrex::BaseFab<int> > m_vector_size;
    std::map<int, amrex::BaseFab<int*> > m_vector_ptrs;

    // used to keep track of whether we need to rebuild the above or not
    bool m_vectors_initialized = false;
    amrex::BoxArray::RefID m_BARef;
    amrex::DistributionMapping::RefID m_DMRef;

private:
  
  void correctCellVectors(int old_index, int new_index, int grid, const ParticleType& p) override;

};

#endif
