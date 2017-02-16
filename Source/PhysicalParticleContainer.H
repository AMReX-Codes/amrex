#ifndef WARPX_PhysicalParticleContainer_H_
#define WARPX_PhysicalParticleContainer_H_

#include <map>

#include <WarpXParticleContainer.H>

class PhysicalParticleContainer
    : public WarpXParticleContainer
{
public:
    PhysicalParticleContainer (amrex::AmrCore* amr_core, int ispecies);
    virtual ~PhysicalParticleContainer () {}

    virtual void AllocData () override;

    virtual void InitData () override;

    virtual void Evolve (int lev,
			 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
			 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
			 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz, amrex::Real t, amrex::Real dt) override;

    virtual void PostRestart () override {}

private:
    // The outermost Array is for AMR level.
    // map's key is grid id and value is an Array. The size of map is the size of local grids.
    // The Array (i.e., the value of the map) has npartdata unique_ptr to another Array.
    // The size of the map value Array is npartdata.
    // The innermost Array holds a particular data (say Ex or gaminv) for particles on a grid.
    // The size of the innermost array is the number of particles.
    amrex::Array<std::map<int, amrex::Array<std::unique_ptr<amrex::Array<amrex::Real> > > > > m_partdata;
};

#endif
