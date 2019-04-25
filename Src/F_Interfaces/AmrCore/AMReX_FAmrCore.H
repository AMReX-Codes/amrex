#ifndef AMREX_FAMRCORE_H_
#define AMREX_FAMRCORE_H_

#include <functional>

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>

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

    // octree stuff
    Vector<BoxArray> octree_leaf_grids;
    Vector<DistributionMapping> octree_leaf_dmap;
    Vector<std::unique_ptr<MultiFab> > octree_leaf_dummy_mf;
    Vector<Vector<int> > octree_li_full_to_leaf; // li: local index
    Vector<Vector<int> > octree_li_leaf_to_full;

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

