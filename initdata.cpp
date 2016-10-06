
#include <WarpX.H>

void
WarpX::InitData ()
{
    BL_PROFILE("WPX::InitData");

    MultiFab dummy_mf(ba_arr[0], 1, 0);
    mypc->Init(dummy_mf);

    if (verbose) mypc->WriteAsciiFile("Particles_before");
}
