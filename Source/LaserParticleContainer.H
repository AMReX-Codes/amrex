#ifndef WARPX_LaserParticleContainer_H_
#define WARPX_LaserParticleContainer_H_

#include <limits>

#include <WarpXParticleContainer.H>

enum class laser_t { Null, Gaussian, Harris };

class LaserParticleContainer
    : public WarpXParticleContainer
{
public:
    LaserParticleContainer (amrex::AmrCore* amr_core, int ispecies);
    virtual ~LaserParticleContainer () {}

    virtual void InitData () final;

    virtual void EvolveES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                                 amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                           amrex::Real t, amrex::Real dt) { BL_ASSERT(false); }

    virtual void Evolve (int lev,
			 const amrex::MultiFab&, const amrex::MultiFab&, const amrex::MultiFab&,
			 const amrex::MultiFab&, const amrex::MultiFab&, const amrex::MultiFab&,
			 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz,
                         amrex::MultiFab* rho, amrex::Real t, amrex::Real dt) final;

    virtual void PostRestart () final;

private:

    // runtime paramters
    laser_t                   profile     = laser_t::Null;
    amrex::Vector<amrex::Real> position;
    amrex::Vector<amrex::Real> nvec;
    amrex::Vector<amrex::Real> p_X;
    long                      pusher_algo = -1;
    amrex::Real               e_max       = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real               wavelength  = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real               gamma_boost = 1;
    amrex::Real               beta_boost = 0;

    // computed using runtime parameters
    amrex::Vector<amrex::Real> p_Y;
    amrex::Vector<amrex::Real> u_X;
    amrex::Vector<amrex::Real> u_Y;
    amrex::Real weight   = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real mobility = std::numeric_limits<amrex::Real>::quiet_NaN();

    // Gaussian profile
    amrex::Real profile_waist          = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real profile_duration       = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real profile_t_peak         = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real profile_focal_distance = std::numeric_limits<amrex::Real>::quiet_NaN();

    // laser particle domain
    amrex::RealBox prob_domain;

    void ComputeSpacing (int lev, amrex::Real& Sx, amrex::Real& Sy) const;
    void ComputeWeightMobility (amrex::Real Sx, amrex::Real Sy);
    void InitData (int lev);
};

#endif
