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

    using make_level_funptr_t = void(*)(int lev, Real time, const BoxArray* ba,
				      const DistributionMapping* dm);
    make_level_funptr_t make_new_level_from_scratch = nullptr;
    make_level_funptr_t make_new_level_from_coarse  = nullptr;
    make_level_funptr_t remake_level                = nullptr;

    using clear_level_funptr_t = void(*)(int lev);
    clear_level_funptr_t clear_level = nullptr;

    using error_est_funptr_t = void(*)(int lev, TagBoxArray* tags, Real time, char tagval, char clearval);
    error_est_funptr_t error_est = nullptr;

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

