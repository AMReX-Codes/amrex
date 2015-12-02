
#include <ParallelDescriptor.H>
#include <SMC.H>

int
main (int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    {
	SMC smc;
	smc.evolve();
    }

    BL_PROFILE_VAR_STOP(pmain);

    BoxLib::Finalize();
}
