
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
}
