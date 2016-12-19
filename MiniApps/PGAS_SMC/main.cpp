
#include <AMReX_ParallelDescriptor.H>
#include <SMC.H>

using namespace amrex;

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    {
	SMC smc;
	smc.evolve();
    }

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
