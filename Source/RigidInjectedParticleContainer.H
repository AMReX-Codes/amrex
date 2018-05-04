#ifndef WARPX_RigidInjectedParticleContainer_H_
#define WARPX_RigidInjectedParticleContainer_H_

#include <PhysicalParticleContainer.H>
#include <AMReX_Vector.H>

class RigidInjectedParticleContainer
    : public PhysicalParticleContainer
{
public:
    RigidInjectedParticleContainer (amrex::AmrCore* amr_core,
                                    int ispecies,
                                    const std::string& name);
    virtual ~RigidInjectedParticleContainer () {}

    virtual void InitData() override;

    virtual void BoostandRemapParticles();

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
                         amrex::MultiFab* rho2,
                         amrex::Real t,
                         amrex::Real dt) override;

    virtual void PushPX(WarpXParIter& pti,
	                amrex::Vector<amrex::Real>& xp,
                        amrex::Vector<amrex::Real>& yp,
                        amrex::Vector<amrex::Real>& zp,
                        amrex::Vector<amrex::Real>& giv,
                        amrex::Real dt) override;

    virtual void PushP (int lev, amrex::Real dt,
                        const amrex::MultiFab& Ex,
                        const amrex::MultiFab& Ey,
                        const amrex::MultiFab& Ez,
                        const amrex::MultiFab& Bx,
                        const amrex::MultiFab& By,
                        const amrex::MultiFab& Bz) override;

private:

    // User input quantities
    amrex::Real zinject_plane = 0.;
    bool projected = false;
    bool focused = true;

    amrex::Real vzbeam_ave_boosted;

    bool done_injecting = false;

    // Temporary quantites
    amrex::Real zinject_plane_previous;
    amrex::Vector<int> done_injecting_temp;

};

#endif
