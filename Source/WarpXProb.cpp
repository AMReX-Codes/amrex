
//
// Each problem must have its own version of WarpX::InitLevelData(int lev)
// to initialize the field data on this level.
//

#include <WarpX.H>

using namespace amrex;

void
WarpX::InitLevelData (int lev)
{
    static_assert(false,
		  "Each problem must have its own version of WarpX::InitLevelData(int lev)");

#if 0
    for (int i = 0; i < 3; ++i) {
	current[lev][i]->setVal(0.0);
	Efield[lev][i]->setVal(0.0);
	Bfield[lev][i]->setVal(0.0);
    }
#endif
}
