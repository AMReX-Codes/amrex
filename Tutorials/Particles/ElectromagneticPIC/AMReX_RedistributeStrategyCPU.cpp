#include <AMReX_RedistributeStrategy.H>

void RedistributeStrategyCPU::Redistribute(std::map<int, Particles>&         a_particles,
                                           const amrex::BoxArray&            a_ba,
                                           const amrex::DistributionMapping& a_dm,
                                           const amrex::Geometry&            a_geom,
                                           const amrex::MultiFab*            a_mask_ptr)
{
}
    
void RedistributeStrategyCPU::OK (std::map<int, Particles>&         a_particles,
                                  const amrex::BoxArray&            a_ba,
                                  const amrex::DistributionMapping& a_dm,
                                  const amrex::Geometry&            a_geom,
                                  const amrex::MultiFab*            a_mask_ptr)
{
}

void RedistributeStrategyCPU::RedistributeMPI (std::map<int, amrex::Vector<char> >& not_ours,
                                               std::map<int, Particles>& a_particles)
{
}
