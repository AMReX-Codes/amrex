
#include <AMReX_FAmrCore.H>

using namespace amrex;

extern "C" {

    void amrex_fi_new_famrcore (FAmrCore*& famrcore)
    {
	famrcore = new FAmrCore();
    }

    void amrex_fi_delete_famrcore (FAmrCore* famrcore)
    {
	delete famrcore;
    }

    int amrex_fi_get_max_level (const FAmrCore* famrcore)
    {
	return famrcore->maxLevel();
    }

    void amrex_fi_get_ref_ratio (int* ref_ratio, const FAmrCore* famrcore)
    {
	int n = famrcore->maxLevel();
	for (int i = 0; i < n; ++i) {
	    ref_ratio[i] = famrcore->MaxRefRatio(i);
	}
    }

    int amrex_fi_get_finest_level (const FAmrCore* famrcore)
    {
	return famrcore->finestLevel();
    }

    void amrex_fi_get_boxarray (const BoxArray*& ba, int lev, const FAmrCore* famrcore)
    {
	const BoxArray& ba_ = famrcore->boxArray(lev);
	ba = &ba_;
    }

    void amrex_fi_get_distromap (const DistributionMapping*& dm, int lev, const FAmrCore* famrcore)
    {
	const DistributionMapping& dm_ = famrcore->DistributionMap(lev);
	dm = &dm_;
    }

    void amrex_fi_get_geometry (const Geometry*& geom, int lev, const FAmrCore* famrcore)
    {
	const Geometry& geom_ = famrcore->Geom(lev);
	geom = &geom_;
    }

    void amrex_fi_make_base_grids (BoxArray*& ba, const FAmrCore* famrcore)
    {
	ba = new BoxArray(famrcore->MakeBaseGrids());
    }

    void amrex_fi_set_finest_level (int new_finest_level, FAmrCore* famrcore)
    {
	famrcore->SetFinestLevel(new_finest_level);
    }

    void amrex_fi_set_boxarray (int lev, const BoxArray* ba, FAmrCore* famrcore)
    {
	if (ba)	{
	    famrcore->SetBoxArray(lev, *ba);
	} else {
	    famrcore->ClearBoxArray(lev);
	}
    }

    void amrex_fi_set_distromap (int lev, const DistributionMapping* dm, FAmrCore* famrcore)
    {
	if (dm) {
	    famrcore->SetDistributionMap(lev, *dm);
	} else {
	    famrcore->ClearDistributionMap(lev);
	}
    }

    void amrex_fi_make_new_grids (int baselev, Real time, int* new_finest, const BoxArray** ba,
				  FAmrCore* famrcore)
    {
	Array<BoxArray> new_grids(baselev+1);
	new_grids[baselev] = *ba[baselev];
	famrcore->MakeNewGrids(baselev, time, *new_finest, new_grids);
	for (int lev = baselev+1; lev <= *new_finest; ++lev) {
	    delete ba[lev];
	    ba[lev] = new BoxArray(new_grids[lev]);
	}
    }

    void amrex_fi_init_from_scratch (Real t, FAmrCore* famrcore)
    {
	famrcore->InitFromScratch(t);
    }
}
