#ifndef AMREX_REDISTRIBUTESTRATEGY_H_
#define AMREX_REDISTRIBUTESTRATEGY_H_

#include <cstdlib>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <Particles.H>

#ifdef AMREX_USE_CUDA
#include <thrust/device_vector.h>
#endif

class RedistributeStrategyBase {

public:

    RedistributeStrategyBase () { m_superparticle_size = PIdx::nattribs * sizeof(amrex::Real); }
    
    virtual ~RedistributeStrategyBase () {}

    virtual void Redistribute (std::map<std::pair<int, int>, Particles>&         a_particles,
                               const amrex::BoxArray&            a_ba,
                               const amrex::DistributionMapping& a_dm,
                               const amrex::Geometry&            a_geom,
                               const amrex::MultiFab*            a_mask_ptr) = 0;
    
    virtual void OK (std::map<std::pair<int, int>, Particles>&         a_particles,
                     const amrex::BoxArray&            a_ba,
                     const amrex::DistributionMapping& a_dm,
                     const amrex::Geometry&            a_geom,
                     const amrex::MultiFab*            a_mask_ptr) = 0;
    
protected:
    
    std::size_t m_superparticle_size;
};

class RedistributeStrategyCPU : public RedistributeStrategyBase {

public:
    
    using RedistributeStrategyBase::RedistributeStrategyBase;
    
    ~RedistributeStrategyCPU () final {}
    
    void Redistribute (std::map<std::pair<int, int>, Particles>&         a_particles,
                       const amrex::BoxArray&            a_ba,
                       const amrex::DistributionMapping& a_dm,
                       const amrex::Geometry&            a_geom,
                       const amrex::MultiFab*            a_mask_ptr) final;
    
    void OK (std::map<std::pair<int, int>, Particles>& a_particles,
             const amrex::BoxArray&            a_ba,
             const amrex::DistributionMapping& a_dm,
             const amrex::Geometry&            a_geom,
             const amrex::MultiFab*            a_mask_ptr) final;

private:

    void RedistributeMPI (std::map<int, amrex::Vector<char> >& not_ours,
                          std::map<std::pair<int, int>, Particles>& a_particles);
};

#ifdef AMREX_USE_CUDA
class RedistributeStrategyGPU : public RedistributeStrategyBase {
    
public:
    
    using RedistributeStrategyBase::RedistributeStrategyBase;
    
    ~RedistributeStrategyGPU () final {}
    
    void Redistribute (std::map<std::pair<int, int>, Particles>&         a_particles,
                       const amrex::BoxArray&            a_ba,
                       const amrex::DistributionMapping& a_dm,
                       const amrex::Geometry&            a_geom,
                       const amrex::MultiFab*            a_mask_ptr) final;
    
    void OK (std::map<std::pair<int, int>, Particles>& a_particles,
             const amrex::BoxArray&            a_ba,
             const amrex::DistributionMapping& a_dm,
             const amrex::Geometry&            a_geom,
             const amrex::MultiFab*            a_mask_ptr) final;
private:

    void RedistributeMPI (std::map<int, thrust::device_vector<char> >& not_ours,
                          std::map<std::pair<int, int>, Particles>& a_particles);
};
#endif // AMREX_USE_CUDA

#endif // AMREX_REDISTRIBUTESTRATEGY_H_
