
#include <cstdlib>
#include <string>

#include <VisMF.H>
#include <Utility.H>

int
main (int argc, char** argv)
{
    BoxLib::Initialize(argc, argv);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Successfully initialized BoxLib" << std::endl;

    VisMF vmf(std::string("SD_0_New_MF"));

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Successfully read the VisMF header" << std::endl;

    MultiFab rho(vmf.boxArray(),1,0);

    for (MFIter mfi(rho); mfi.isValid(); ++mfi)
    {
        std::cout << "Attempting to read fab at index: " << mfi.index() << std::endl;

        rho[mfi].copy(vmf.GetFab(mfi.index(),3),0,0,1);

        vmf.clear(mfi.index(),3);
    }

    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Attempting to write Rho ..." << std::endl;

    VisMF::Write(rho, std::string("Rho"));

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Successfully wrote Rho ..." << std::endl;

    BoxLib::Finalize();
}
