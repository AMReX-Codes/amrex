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
        w = 0,      // Weight
        vx, vy, vz, // Velocities
        Ex, Ey, Ez, // Electric fields, interpolated to the particle positions
        nattribs
    };
};

class MyParIter
    : public amrex::ParIter<0, 0, PIdx::nattribs>
{
public:
    using amrex::ParIter<0, 0, PIdx::nattribs>::ParIter;

    const std::array<amrex::Array<amrex::Real>, PIdx::nattribs>& GetAttribs () const {
        return GetStructOfArrays().GetRealData();
    }

    std::array<amrex::Array<amrex::Real>, PIdx::nattribs>& GetAttribs () {
        return GetStructOfArrays().GetRealData();
    }

    const amrex::Array<amrex::Real>& GetAttribs (int comp) const {
        return GetStructOfArrays().GetRealData(comp);
    }
    
    amrex::Array<amrex::Real>& GetAttribs (int comp) {
        return GetStructOfArrays().GetRealData(comp);
    }

};

using VectorMeshData = Array<std::array<std::unique_ptr<MultiFab>, 3> >;
using ScalarMeshData = Array<std::unique_ptr<MultiFab> >;

class ElectrostaticParticleContainer
    : public ParticleContainer<0, 0, PIdx::nattribs>
{
public:

    ElectrostaticParticleContainer(const Array<Geometry>            & geom, 
                                   const Array<DistributionMapping> & dmap,
                                   const Array<BoxArray>            & ba,
                                   const Array<int>                 & rr)
        : ParticleContainer<0, 0, PIdx::nattribs> (geom, dmap, ba, rr) 
        {
            charge = -PhysConst::q_e;
            mass   =  PhysConst::m_e;
        }
    
    void InitParticles();

    void pushX(const Real& dt);

    void DepositCharge(ScalarMeshData& rho);

    void FieldGather(const VectorMeshData& E,
                     const Array<std::unique_ptr<FabArray<BaseFab<int> > > >& masks);

    void Evolve(const VectorMeshData& E, ScalarMeshData& rho, const Real& dt);

    void writeParticles(int n);

private:
    
    Real charge, mass;

};

}

#endif // _ELECTROSTATICPARTICLECONTAINER_H_
