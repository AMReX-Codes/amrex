
#include <ParallelDescriptor.H>
#include <SMC.H>

int
main (int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    SMC smc;

    smc.evolve();

    BoxLib::Finalize();
}
