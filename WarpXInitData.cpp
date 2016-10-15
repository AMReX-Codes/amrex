
#include <WarpX.H>

void
WarpX::InitData ()
{
    BL_PROFILE("WarpX::InitData()");

    MultiFab dummy_mf(ba_arr[0], 1, 0);
    mypc->Init(dummy_mf);

//    if (verbose) mypc->WriteAsciiFile("Particles_before");

    for (int i = 0; i < BL_SPACEDIM; ++i) {
	current[i]->setVal(0.0);
	Efield[i]->setVal(0.0);
	Bfield[i]->setVal(0.0);
    }

    if (plot_int > 0) {
	WritePlotFile(0, 0.0);
    }
}
