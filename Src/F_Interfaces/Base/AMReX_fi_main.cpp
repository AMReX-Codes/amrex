
#include <AMReX.H>
#include <AMReX_BLProfiler.H>

extern "C" { 
    void amrex_fi_init();
    void amrex_fi_finalize();
    void amrex_fmain();
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex_fi_init();

    BL_PROFILE_VAR("main()", pmain);

    amrex_fmain();

    BL_PROFILE_VAR_STOP(pmain);

    amrex_fi_finalize();

    amrex::Finalize();
}
