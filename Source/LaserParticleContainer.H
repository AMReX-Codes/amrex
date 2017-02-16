#ifndef WARPX_LaserParticleContainer_H_
#define WARPX_LaserParticleContainer_H_

#include <limits>

#include <WarpXParticleContainer.H>

enum class laser_t { Null, Gaussian };

class LaserParticleContainer
    : public WarpXParticleContainer
{
public:
    LaserParticleContainer (AmrCore* amr_core, int ispecies);
    virtual ~LaserParticleContainer () {}

    virtual void AllocData () override;

    virtual void InitData () override;

    virtual void Evolve (int lev,
			 const MultiFab&, const MultiFab&, const MultiFab&,
			 const MultiFab&, const MultiFab&, const MultiFab&,
			 MultiFab& jx, MultiFab& jy, MultiFab& jz, Real t, Real dt) override;

    virtual void PostRestart () override;

private:

    // runtime paramters
    laser_t     profile     = laser_t::Null;
    Array<Real> position;
    Array<Real> nvec;
    Array<Real> p_X;
    long        pusher_algo = -1;
    Real        e_max       = std::numeric_limits<Real>::quiet_NaN();
    Real        wavelength  = std::numeric_limits<Real>::quiet_NaN();

    // computed using runtime parameters
    Array<Real> p_Y;
    Array<Real> u_X;
    Array<Real> u_Y;
    Real weight   = std::numeric_limits<Real>::quiet_NaN();
    Real mobility = std::numeric_limits<Real>::quiet_NaN();

    // Gaussian profile
    Real profile_waist          = std::numeric_limits<Real>::quiet_NaN();
    Real profile_duration       = std::numeric_limits<Real>::quiet_NaN();
    Real profile_t_peak         = std::numeric_limits<Real>::quiet_NaN();
    Real profile_focal_distance = std::numeric_limits<Real>::quiet_NaN();

    void ComputeSpacing (Real& Sx, Real& Sy) const;
    void ComputeWeightMobility (Real Sx, Real Sy);
};

#endif
