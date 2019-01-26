#ifdef NEUTRINO_PARTICLES
#ifndef _NeutrinoParticleContainer_H_
#define _NeutrinoParticleContainer_H_

#include "AMReX_Amr.H"
#include "AMReX_AmrLevel.H"
#include "AMReX_AmrParticles.H"
#include "NyxParticleContainer.H"

class NeutrinoParticleContainer
    : public NyxParticleContainer<2+BL_SPACEDIM>
{

    int         m_relativistic; // if 1 then we weight the mass by gamma in AssignDensity*
    amrex::Real        m_csq;          // the square of the speed of light -- used to compute relativistic effects

public:
    NeutrinoParticleContainer (amrex::Amr* amr)
        : NyxParticleContainer<2+BL_SPACEDIM>(amr)
    { }

    virtual ~NeutrinoParticleContainer () {}

    void SetRelativistic (int relativistic) { m_relativistic = relativistic; }

    void SetCSquared (amrex::Real csq) { m_csq = csq; }

    void AssignDensity (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const;

    void AssignRelativisticDensitySingleLevel (amrex::MultiFab& mf, int level, int ncomp=1, int particle_lvl_offset = 0) const;
    
    void AssignRelativisticDensity (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1, int finest_level = -1) const;

};

#endif /*_NeutrinoParticleContainer_H_*/
#endif /*NEUTRINO_PARTICLES*/
