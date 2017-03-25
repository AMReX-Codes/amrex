#include <WarpX.H>

using namespace amrex;

void
WarpX::InitLevelData (int lev)
{
    for (int i = 0; i < 3; ++i) {
	current[lev][i]->setVal(0.0);
	Efield[lev][i]->setVal(0.0);
	Bfield[lev][i]->setVal(0.0);
    }
}
