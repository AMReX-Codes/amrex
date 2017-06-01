#include <WarpX.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
WarpX::InitLevelData (int lev, Real time)
{
    for (int i = 0; i < 3; ++i) {
	current_fp[lev][i]->setVal(0.0);
	Efield_fp[lev][i]->setVal(0.0);
	Bfield_fp[lev][i]->setVal(B_init[i]);
    }

    if (lev > 0) {
        for (int i = 0; i < 3; ++i) {
            Efield_aux[lev][i]->setVal(0.0);
            Bfield_aux[lev][i]->setVal(B_init[i]);

            current_cp[lev][i]->setVal(0.0);
            Efield_cp[lev][i]->setVal(0.0);
            Bfield_cp[lev][i]->setVal(B_init[i]);
        }
    }
}
