#ifndef AMREX_FAMRCORE_H_
#define AMREX_FAMRCORE_H_

#include <functional>

#include <AMReX_AmrCore.H>

namespace amrex {

class FAmrCore
    : public AmrCore
{
public:
    FAmrCore ();
    virtual ~FAmrCore ();

    void (*make_new_level_from_scratch) (int lev, Real time, const BoxArray* ba,
					 const DistributionMapping* dm) = nullptr;
    void (*make_new_level_from_coarse) (int lev, Real time, const BoxArray* ba,
					const DistributionMapping* dm) = nullptr;
    void (*remake_level) (int lev, Real time, const BoxArray* ba,
			  const DistributionMapping* dm) = nullptr;
    void (*clear_level) (int lev) = nullptr;
    void (*error_est) (int lev, TagBoxArray* tags, Real time, int ngrow) = nullptr;

    void SetFinestLevel (int new_finest) { finest_level = new_finest; }

protected:

    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
					  const DistributionMapping& dm) override;
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
					 const DistributionMapping& dm) override;
    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba,
			      const DistributionMapping& dm) override;
    virtual void ClearLevel (int lev) override;
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override;
};

}

#endif

