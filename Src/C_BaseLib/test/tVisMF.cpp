//BL_COPYRIGHT_NOTICE

//
// $Id: tVisMF.cpp,v 1.19 1999-05-10 18:54:24 car Exp $
//

#include <stdlib.h>

#include <VisMF.H>
#include <Utility.H>

static int nProcs = 1;

static int nBoxs  = 10;

static char* the_prog_name;

static aString PerFab("PerFab");

static aString PerCPU("PerCPU");

//
// How defaults to PerCPU.
//
static aString How(PerCPU);

static
void
usage ()
{
    std::cout << "usage: "
              << the_prog_name
              << " [-how PerFab|PerCPU]"
              << " [-nprocs N]"
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
        if (strcmp(*argv, "-nprocs") ==  0)
        {
            if (*++argv)
            {
                nProcs = atoi(*argv);

                if (nProcs <= 0)
                {
                    std::cout << "nprocs must be positive" << std::endl;
                    usage();
                }
            }
            else
            {
                std::cout << "No argument to -nprocs supplied.\n";
                usage();
            }
        }
        else if (strcmp(*argv, "-nboxs") ==  0)
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
                std::cout << "No argument to -nprocs supplied.\n";
                usage();
            }
        }
        else if (strcmp(*argv, "-how") ==  0)
        {
            if (*++argv)
            {
                How = *argv;

                if (!(How == PerCPU || How == PerFab))
                {
                    std::cout << "Invalid value for -how argument\n";
                    usage();
                }
            }
            else
            {
                std::cout << "No argument to -how supplied.\n";
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
              const aString&  mf_name,
              VisMF::How      how)
{
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Writing the MultiFab to disk ...\n";
    }

    double start, end;

    ParallelDescriptor::Synchronize();

    if (ParallelDescriptor::IOProcessor())
    {
        start = Utility::wsecond();
    }

    switch (how)
    {
    case VisMF::OneFilePerCPU:
        VisMF::Write(mf, mf_name, VisMF::OneFilePerCPU); break;
    case VisMF::OneFilePerFab:
        VisMF::Write(mf, mf_name, VisMF::OneFilePerFab); break;
    default:
        BoxLib::Error("Bad case in switch");
    }

    ParallelDescriptor::Synchronize();

    if (ParallelDescriptor::IOProcessor())
    {
        end = Utility::wsecond();

        std::cout << "\nWallclock time for MF write: " << (end-start) << '\n';

        std::cout << "Reading the MultiFab from disk ...\n";
    }

    VisMF vmf(mf_name);

    BL_ASSERT(vmf.length() == mf.boxArray().length());

    for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
    {
        //const FArrayBox& fab = vmf[mfi.index()];
        const FArrayBox& fab = vmf.GetFab(mfi.index(), 0);

        std::cout << "\tCPU #"
                  << ParallelDescriptor::MyProc()
                  << " read FAB #"
                  << mfi.index()
                  << '\n';
    }

    ParallelDescriptor::Synchronize();

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Building new MultiFab from disk version ....\n\n";
    }

    MultiFab new_mf;
    
    VisMF::Read(new_mf, mf_name);
}

int
main (int, char** argv)
{
    the_prog_name = argv[0];

    parse_args(argv);

    ParallelDescriptor::StartParallel(nProcs, &argc, &argv);

    BoxArray ba(nBoxs);

    ba.set(0, Box(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(2,2,2))));

    for (int i = 1; i < nBoxs; i++)
    {
        ba.set(i,grow(ba[i-1],2));
    }

    MultiFab mf(ba, 2, 1);

    for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
    {
        mf[mfi.index()].setVal(mfi.index()+1);
    }
    //
    // Set cells in ghost region to zero.
    //
    mf.setBndry(0);

    static const aString mf_name = "Spam-n-Eggs";

    Write_N_Read (mf,
                  mf_name,
                  (How==PerCPU) ? VisMF::OneFilePerCPU : VisMF::OneFilePerFab);

    ParallelDescriptor::EndParallel();
}
