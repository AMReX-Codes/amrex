#include <AMReX_FPhysBC.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_new_physbc (FPhysBC*& pbc, FPhysBC::fill_physbc_funptr_t fill)
    {
	pbc = new FPhysBC(fill);
    }

    void amrex_fi_delete_physbc (FPhysBC* pbc)
    {
	delete pbc;
    }
}
