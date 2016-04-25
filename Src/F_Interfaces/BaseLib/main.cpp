
#include <BoxLib.H>
#include <BLProfiler.H>

extern "C" { void fmain(); }

int main (int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    fmain();

    BL_PROFILE_VAR_STOP(pmain);

    BoxLib::Finalize();
}
