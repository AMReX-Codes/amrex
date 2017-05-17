// ----------------------------------
// email push test 0
// ----------------------------------

#include <cstdlib>
#include <string>

#include <AMReX_VisMF.H>
#include <AMReX_Utility.H>

using namespace amrex;

static int nBoxs  = 10;

static char* the_prog_name;


static
void
usage ()
{
    std::cout << "usage: "
              << the_prog_name
              << " [-nboxs N]"
              << std::endl;
    exit(1);
}

static
void
parse_args (char**& argv)
{
    while (*++argv && **argv == '-')
    {
        if (strcmp(*argv, "-nboxs") ==  0)
        {
            if (*++argv)
            {
                nBoxs = atoi(*argv);

                if (nBoxs <= 0)
                {
                    std::cout << "nboxs must be positive" << std::endl;
                    usage();
                }
            }
            else
            {
                std::cout << "No argument to -nboxs supplied.\n";
                usage();
            }
        }
        else
        {
            std::cout << "Exiting, unknown option: " << *argv << std::endl;
            usage();
        }
    }
}

static
void
Write_N_Read (const MultiFab& mf,
              const std::string&  mf_name)
{
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Writing the MultiFab to disk ...\n";
    }

    double start, end;

    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        start = ParallelDescriptor::second();
    }

    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        end = ParallelDescriptor::second();

        std::cout << "\nWallclock time for MF write: " << (end-start) << '\n';

        std::cout << "Reading the MultiFab from disk ...\n";
    }

    VisMF vmf(mf_name);

    BL_ASSERT(vmf.size() == mf.boxArray().size());

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        //const FArrayBox& fab = vmf[mfi.index()];
        const FArrayBox& fab = vmf.GetFab(mfi.index(), 0);

        std::cout << "\tCPU #"
                  << ParallelDescriptor::MyProc()
                  << " read FAB #"
                  << mfi.index()
                  << '\n';
    }

    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Building new MultiFab from disk version ....\n\n";
    }

    MultiFab new_mf;
    
    VisMF::Read(new_mf, mf_name);
}

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    the_prog_name = argv[0];
    parse_args(argv);

    BoxArray ba(nBoxs);

    ba.set(0, Box(IntVect(AMREX_D_DECL(0,0,0)), IntVect(AMREX_D_DECL(2,2,2))));

    for (int i = 1; i < nBoxs; i++)
    {
        ba.set(i,amrex::grow(ba[i-1],2));
    }

    DistributionMapping dm{ba};

    MultiFab mf(ba, dm, 2, 1);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        mf[mfi.index()].setVal(mfi.index()+1);
    }
    //
    // Set cells in ghost region to zero.
    //
    mf.setBndry(0);

    static const std::string mf_name = "Spam-n-Eggs";

    Write_N_Read (mf,
                  mf_name);

    amrex::Finalize();
}
