#include <WarpX.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
WarpX::InitLevelData (int lev, Real time)
{
    Array<Real> B0(3,0.0);
    ParmParse pp("warpx");
    pp.queryarr("B_init", B0, 0, 3);

    for (int i = 0; i < 3; ++i) {
	current_fp[lev][i]->setVal(0.0);
	Efield_fp[lev][i]->setVal(0.0);
	Bfield_fp[lev][i]->setVal(B0[i]);
    }

    if (lev > 0) {
        for (int i = 0; i < 3; ++i) {
            Efield_aux[lev][i]->setVal(0.0);
            Bfield_aux[lev][i]->setVal(B0[i]);

            current_cp[lev][i]->setVal(0.0);
            Efield_cp[lev][i]->setVal(0.0);
            Bfield_cp[lev][i]->setVal(B0[i]);
        }
    }
}
