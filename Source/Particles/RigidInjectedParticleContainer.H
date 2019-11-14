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

    virtual void RemapParticles();
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
                         amrex::MultiFab* cjx,
                         amrex::MultiFab* cjy,
                         amrex::MultiFab* cjz,
                         amrex::MultiFab* rho,
                         amrex::MultiFab* crho,
                         const amrex::MultiFab* cEx,
                         const amrex::MultiFab* cEy,
                         const amrex::MultiFab* cEz,
                         const amrex::MultiFab* cBx,
                         const amrex::MultiFab* cBy,
                         const amrex::MultiFab* cBz,
                         amrex::Real t,
                         amrex::Real dt,
                         DtType a_dt_type=DtType::Full) override;

    virtual void PushPX(WarpXParIter& pti,
                        amrex::Gpu::ManagedDeviceVector<amrex::ParticleReal>& xp,
                        amrex::Gpu::ManagedDeviceVector<amrex::ParticleReal>& yp,
                        amrex::Gpu::ManagedDeviceVector<amrex::ParticleReal>& zp,
                        amrex::Real dt, DtType a_dt_type=DtType::Full) override;

    virtual void PushP (int lev, amrex::Real dt,
                        const amrex::MultiFab& Ex,
                        const amrex::MultiFab& Ey,
                        const amrex::MultiFab& Ez,
                        const amrex::MultiFab& Bx,
                        const amrex::MultiFab& By,
                        const amrex::MultiFab& Bz) override;

    virtual void ReadHeader (std::istream& is) override;

    virtual void WriteHeader (std::ostream& os) const override;

private:

    // User input quantities
    amrex::Real zinject_plane = 0.;
    bool projected = true; // When true, particle transverse positions are directly projected (without adjusment)
    bool focused = false; // When true, particle transverse positions are adjusted to account for distance between zinject and z=0
    bool rigid_advance = true; // When true, particles are advance with vzbar before injection

    amrex::Real vzbeam_ave_boosted;

    amrex::Vector<int> done_injecting;
    amrex::Vector<amrex::Real> zinject_plane_levels;

    // Temporary quantites
    amrex::Real zinject_plane_lev;
    amrex::Real zinject_plane_lev_previous;
    bool done_injecting_lev;

};

#endif
