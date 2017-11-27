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

    virtual void InitData () final;

    virtual void FieldGatherES(const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                               const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks) override;

    virtual void FieldGather(int lev,
                             const amrex::MultiFab& Ex,
                             const amrex::MultiFab& Ey,
                             const amrex::MultiFab& Ez,
                             const amrex::MultiFab& Bx,
                             const amrex::MultiFab& By,
                             const amrex::MultiFab& Bz) final;

    virtual void EvolveES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                                 amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                           amrex::Real t, amrex::Real dt) override;

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
                         amrex::MultiFab* rho,
                         amrex::Real t,
                         amrex::Real dt) final;

    virtual void PostRestart () final {}

    // Inject particles in Box 'part_box'
    void AddParticles (int lev, amrex::Box part_box = amrex::Box());

    void AddGaussianBeam(amrex::Real x_m, amrex::Real y_m, amrex::Real z_m,
                         amrex::Real x_rms, amrex::Real y_rms, amrex::Real z_rms,
                         amrex::Real q_tot, long npart);

protected:

    std::string species_name;
    std::unique_ptr<PlasmaInjector> plasma_injector;
};

#endif
