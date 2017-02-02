
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
}
