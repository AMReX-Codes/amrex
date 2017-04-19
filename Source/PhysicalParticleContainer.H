#ifndef WARPX_PhysicalParticleContainer_H_
#define WARPX_PhysicalParticleContainer_H_

#include <map>

#include <PlasmaInjector.H>
#include <WarpXParticleContainer.H>

class PhysicalParticleContainer
    : public WarpXParticleContainer
{
public:
    PhysicalParticleContainer (amrex::AmrCore* amr_core,
                               int ispecies,
                               const std::string& name);
    virtual ~PhysicalParticleContainer () {}

    virtual void InitData () override;

    virtual void FieldGather(int lev,
                             const amrex::MultiFab& Ex,
                             const amrex::MultiFab& Ey,
                             const amrex::MultiFab& Ez,
                             const amrex::MultiFab& Bx,
                             const amrex::MultiFab& By,
                             const amrex::MultiFab& Bz) override;

    virtual void Evolve (int lev,
			 const amrex::MultiFab& Ex,
                         const amrex::MultiFab& Ey,
                         const amrex::MultiFab& Ez,
			 const amrex::MultiFab& Bx,
                         const amrex::MultiFab& By,
                         const amrex::MultiFab& Bz,
			 amrex::MultiFab& jx,
                         amrex::MultiFab& jy,
                         amrex::MultiFab& jz,
                         amrex::Real t,
                         amrex::Real dt) override;

    virtual void PostRestart () override {}

    // Inject particles in Box 'part_box'
    void AddParticles (int lev, amrex::Box part_box = amrex::Box());

protected:

    std::string species_name;
    std::unique_ptr<PlasmaInjector> plasma_injector;
};

#endif
