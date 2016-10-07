
#include <WarpX.H>

void
WarpX::InitData ()
{
    BL_PROFILE("WarpX::InitData()");

    MultiFab dummy_mf(ba_arr[0], 1, 0);
    mypc->Init(dummy_mf);

    if (verbose) mypc->WriteAsciiFile("Particles_before");

    for (int i = 0; i < BL_SPACEDIM; ++i) {
	Efield[i]->setVal(0.0);
	Bfield[i]->setVal(0.0);
    }
}
