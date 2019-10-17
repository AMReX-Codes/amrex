#ifndef _ELECTROSTATICPARTICLECONTAINER_H_
#define _ELECTROSTATICPARTICLECONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"

#include "PhysConst.H"

namespace amrex {

struct PIdx
{
    enum { // Particle Attributes stored in amrex::ParticleContainer's struct-of-arrays
        w = 0,
        vx, vy,
#if BL_SPACEDIM == 3
        vz,
#endif
        Ex, Ey, 
#if BL_SPACEDIM == 3
        Ez,
#endif
        nattribs
    };
};

class MyParIter
    : public amrex::ParIter<0, 0, PIdx::nattribs>
{
public:
    using amrex::ParIter<0, 0, PIdx::nattribs>::ParIter;
    using RealVector = ParIter::ContainerType::RealVector;
    
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

using VectorMeshData = Vector<std::array<std::unique_ptr<MultiFab>, BL_SPACEDIM> >;
using ScalarMeshData = Vector<std::unique_ptr<MultiFab> >;

class ElectrostaticParticleContainer
    : public ParticleContainer<0, 0, PIdx::nattribs>
{
public:

    ElectrostaticParticleContainer(const Vector<Geometry>            & geom, 
                                   const Vector<DistributionMapping> & dmap,
                                   const Vector<BoxArray>            & ba,
                                   const Vector<int>                 & rr)
        : ParticleContainer<0, 0, PIdx::nattribs> (geom, dmap, ba, rr) 
        {
            charge = -PhysConst::q_e;
            mass   =  PhysConst::m_e;
        }
    
    void InitParticles();

    void pushX(const Real& dt);

    void DepositCharge(ScalarMeshData& rho);

    void FieldGather(const VectorMeshData& E,
                     const Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks);

    void Evolve(const VectorMeshData& E, ScalarMeshData& rho, const Real& dt);

    void writeParticles(int n);

private:
    
    Real charge, mass;

};

}

#endif // _ELECTROSTATICPARTICLECONTAINER_H_
