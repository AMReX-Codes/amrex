
#include <AMReX.H>
#include <AMReX_BLProfiler.H>

extern "C" { 
    void amrex_fmain();
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    amrex_fmain();

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
