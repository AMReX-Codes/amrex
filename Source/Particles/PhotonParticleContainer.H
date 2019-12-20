#ifndef WARPX_PhotonParticleContainer_H_
#define WARPX_PhotonParticleContainer_H_

#include <PhysicalParticleContainer.H>
#include <AMReX_Vector.H>

/**
 * Photon particles have no mass, they deposit no charge, and see specific QED
 * effects. For these reasons, they are stored in the separate particle
 * container PhotonParticleContainer, that inherts from
 * PhysicalParticleContainer. The particle pusher and current deposition, in
 * particular, are overriden in this container.
 */
class PhotonParticleContainer
    : public PhysicalParticleContainer
{
public:
    PhotonParticleContainer (amrex::AmrCore* amr_core,
                                    int ispecies,
                                    const std::string& name);
    virtual ~PhotonParticleContainer () {}

    virtual void InitData() override;

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

    // Do nothing
    virtual void PushP (int lev,
                        amrex::Real dt,
                        const amrex::MultiFab& Ex,
                        const amrex::MultiFab& Ey,
                        const amrex::MultiFab& Ez,
                        const amrex::MultiFab& Bx,
                        const amrex::MultiFab& By,
                        const amrex::MultiFab& Bz) override {};


    // DepositCurrent should do nothing for photons
    virtual void DepositCurrent(WarpXParIter& pti,
                                RealVector& wp,
                                RealVector& uxp,
                                RealVector& uyp,
                                RealVector& uzp,
                                const int * const ion_lev,
                                amrex::MultiFab* jx,
                                amrex::MultiFab* jy,
                                amrex::MultiFab* jz,
                                const long offset,
                                const long np_to_depose,
                                int thread_num,
                                int lev,
                                int depos_lev,
                                amrex::Real dt) override {};
    //Photons are not leptons
    virtual bool AmIALepton () override
    {
        return false;
    };

#ifdef WARPX_QED
    /**
     * This function evolves the optical depth of the photons if QED effects
     * are enabled.
     * @param[in,out] pti particle iterator (optical depth will be modified)
     * @param[in] dt temporal step
     */
    virtual void EvolveOpticalDepth(WarpXParIter& pti,
        amrex::Real dt) override;

#endif

};

#endif // #ifndef WARPX_PhotonParticleContainer_H_
