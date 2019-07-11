
#include <AMReX_FAmrCore.H>

using namespace amrex;

extern "C" {

    void amrex_fi_new_amrcore (FAmrCore*& amrcore)
    {
	amrcore = new FAmrCore();
    }

    void amrex_fi_delete_amrcore (FAmrCore* amrcore)
    {
	delete amrcore;
    }

    int amrex_fi_get_max_level (const FAmrCore* amrcore)
    {
	return amrcore->maxLevel();
    }

    void amrex_fi_get_ref_ratio (int* ref_ratio, const FAmrCore* amrcore)
    {
	int n = amrcore->maxLevel();
	for (int i = 0; i < n; ++i) {
	    ref_ratio[i] = amrcore->MaxRefRatio(i);
	}
    }

    int amrex_fi_get_finest_level (const FAmrCore* amrcore)
    {
	return amrcore->finestLevel();
    }

    void amrex_fi_get_boxarray (const BoxArray*& ba, int lev, const FAmrCore* amrcore)
    {
	const BoxArray& ba_ = amrcore->boxArray(lev);
	ba = &ba_;
    }

    void amrex_fi_get_distromap (const DistributionMapping*& dm, int lev, const FAmrCore* amrcore)
    {
	const DistributionMapping& dm_ = amrcore->DistributionMap(lev);
	dm = &dm_;
    }

    void amrex_fi_get_geometry (const Geometry*& geom, int lev, const FAmrCore* amrcore)
    {
	const Geometry& geom_ = amrcore->Geom(lev);
	geom = &geom_;
    }

    void amrex_fi_make_base_grids (BoxArray*& ba, const FAmrCore* amrcore)
    {
	ba = new BoxArray(amrcore->MakeBaseGrids());
    }

    void amrex_fi_set_finest_level (int new_finest_level, FAmrCore* amrcore)
    {
	amrcore->SetFinestLevel(new_finest_level);
    }

    void amrex_fi_set_boxarray (int lev, const BoxArray* ba, FAmrCore* amrcore)
    {
	if (ba)	{
	    amrcore->SetBoxArray(lev, *ba);
	} else {
	    amrcore->ClearBoxArray(lev);
	}
    }

    void amrex_fi_set_distromap (int lev, const DistributionMapping* dm, FAmrCore* amrcore)
    {
	if (dm) {
	    amrcore->SetDistributionMap(lev, *dm);
	} else {
	    amrcore->ClearDistributionMap(lev);
	}
    }

    void amrex_fi_set_geometry (int lev, const Geometry* gm, FAmrCore* amrcore)
    {
        if (gm) {
            amrcore->SetGeometry(lev, *gm);
        }
    }

    void amrex_fi_make_new_grids (int baselev, Real time, int* new_finest, const BoxArray** ba,
				  FAmrCore* amrcore)
    {
	Vector<BoxArray> new_grids(baselev+1);
	new_grids[baselev] = *ba[baselev];
	amrcore->MakeNewGrids(baselev, time, *new_finest, new_grids);
	for (int lev = baselev+1; lev <= *new_finest; ++lev) {
	    delete ba[lev];
	    ba[lev] = new BoxArray(new_grids[lev]);
	}
    }

    void amrex_fi_init_from_scratch (Real t, FAmrCore* amrcore)
    {
	amrcore->InitFromScratch(t);
    }

    void amrex_fi_init_virtual_functions (FAmrCore::make_level_funptr_t mk_lev_scrtch,
                                          FAmrCore::make_level_funptr_t mk_lev_crse,
                                          FAmrCore::make_level_funptr_t mk_lev_re,
					  FAmrCore::clear_level_funptr_t clr_lev,
					  FAmrCore::error_est_funptr_t err_est,
					  FAmrCore* amrcore)
    {
	amrcore->make_new_level_from_scratch = mk_lev_scrtch;
        amrcore->make_new_level_from_coarse  = mk_lev_crse;
        amrcore->remake_level                = mk_lev_re;
	amrcore->clear_level                 = clr_lev;
	amrcore->error_est                   = err_est;
    }

    void amrex_fi_regrid (int baselev, Real t, FAmrCore* amrcore)
    {
	amrcore->regrid(baselev, t);
    }
}
