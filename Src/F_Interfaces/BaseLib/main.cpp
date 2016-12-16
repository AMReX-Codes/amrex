
#include <AMReX_BoxLib.H>
#include <AMReX_BLProfiler.H>

extern "C" { void fmain(); }

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    fmain();

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
