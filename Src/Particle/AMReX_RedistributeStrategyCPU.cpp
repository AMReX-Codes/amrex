#include <AMReX_RedistributeStrategy.H>

void RedistributeStrategyCPU::Redistribute(std::map<std::pair<int, int>, Particles>&         a_particles,
                                           const amrex::BoxArray&            a_ba,
                                           const amrex::DistributionMapping& a_dm,
                                           const amrex::Geometry&            a_geom,
                                           const amrex::MultiFab*            a_mask_ptr)
{
}
    
void RedistributeStrategyCPU::OK (std::map<std::pair<int, int>, Particles>&         a_particles,
                                  const amrex::BoxArray&            a_ba,
                                  const amrex::DistributionMapping& a_dm,
                                  const amrex::Geometry&            a_geom,
                                  const amrex::MultiFab*            a_mask_ptr)
{
}

void RedistributeStrategyCPU::RedistributeMPI (std::map<int, amrex::Vector<char> >& not_ours,
                                               std::map<std::pair<int, int>, Particles>& a_particles)
{
}
