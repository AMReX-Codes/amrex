
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
}
