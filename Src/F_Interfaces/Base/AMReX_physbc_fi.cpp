#include <AMReX_FPhysBC.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_new_physbc (FPhysBC*& pbc, FPhysBC::fill_physbc_funptr_t fill, const Geometry* geom)
    {
	pbc = new FPhysBC(fill, geom);
    }

    void amrex_fi_delete_physbc (FPhysBC* pbc)
    {
	delete pbc;
    }
}
