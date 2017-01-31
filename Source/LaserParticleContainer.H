#ifndef WARPX_LaserParticleContainer_H_
#define WARPX_LaserParticleContainer_H_

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
private:
    laser_t     profile     = laser_t::Null;
    Array<Real> position;
    Array<Real> nvec;
    Array<Real> p_X;
    Array<Real> p_Y;
    Array<Real> u_X;
    Array<Real> u_Y;
    long        pusher_algo = -1;
    Real e_max;
    Real profile_waist;
    Real profile_duration;
    Real profile_t_peak;
    Real profile_focal_distance;
    Real wavelength;
    Real mobility;

    Real ComputeWeight (Real Sx, Real Sy) const;
};

#endif
