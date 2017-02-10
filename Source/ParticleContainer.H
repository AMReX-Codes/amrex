
#ifndef WARPX_ParticleContainer_H_
#define WARPX_ParticleContainer_H_

#include <WarpXParticleContainer.H>
#include <PhysicalParticleContainer.H>
#include <LaserParticleContainer.H>

//
// MultiParticleContainer holds multiple (nspecies or npsecies+1 when
// using laser) WarpXParticleContainer.  WarpXParticlContainer is
// derived from amrex::ParticleContainer, and it is
// polymorphic. PhysicalParticleContainer and LaserParticleContainer are
// concrete class dervied from WarpXParticlContainer.
//

class MultiParticleContainer
{

public:

    MultiParticleContainer (AmrCore* amr_core);

    ~MultiParticleContainer() {}

    WarpXParticleContainer& GetParticleContainer (int ispecies) {
	return *allcontainers[ispecies];
    }

    void AllocData ();

    void InitData ();

    void Evolve (int lev,
		 const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
		 const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
		 MultiFab& jx, MultiFab& jy, MultiFab& jz, Real t, Real dt);

    std::unique_ptr<MultiFab> GetChargeDensity(int lev, bool local = false);

    void Checkpoint (const std::string& dir,
		     const std::string& name,
		     bool is_checkpoint) const;

    void Restart (const std::string& dir, const std::string& name);

    void PostRestart ();

    void ReadHeader (std::istream& is);

    void WriteHeader (std::ostream& os) const;

    void Redistribute (bool where_called=false);

    Array<long> NumberOfParticlesInGrid(int lev) const;

    void Increment (MultiFab& mf, int lev);

    void SetParticleBoxArray (int lev, BoxArray& new_ba);
    void SetParticleDistributionMap (int lev, DistributionMapping& new_dm);

private:

    // physical particles (+ laser)
    Array<std::unique_ptr<WarpXParticleContainer> > allcontainers;

    void ReadParameters ();

    // runtime parameters
    int nspecies = 1;   // physical particles only. If WarpX::use_laser, nspecies+1 == allcontainers.size().
};
#endif /*WARPX_ParticleContainer_H_*/
