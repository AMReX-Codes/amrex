
#include <unistd.h>

#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <ArrayView.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <TV_TempWrite.H>



static
void
PrintUsage(int argc, char *argv[])
{
    cout << "Usage: " << endl;
    cout << argv[0] << " infile [options] \n\tOptions:" << endl;
    cout << "\t   -ascii   \t[if set, dump ascii mfab to stdout]" << endl;
    cout << endl;
    exit(0);
}


int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(&argc, &argv);

//
//  Parse the command line
//
    if (argc < 2)
        PrintUsage(argc,argv);

    if (argv[1][0] == '-')
    {
        cerr << "input file must be first argument\n";
        PrintUsage(argc, argv);
    }

    ParmParse pp(argc-2,argv+2);

    if (pp.contains("help"))
        PrintUsage(argc, argv);

    aString iFile = argv[1];

    bool ascii = false;
    if (pp.contains("ascii"))
        ascii = true;
//
//  Read multifab
//
    MultiFab mf;
    readMF(mf,iFile.c_str());

    int ngrow = mf.nGrow();
    MultiFab tmp(mf.boxArray(),mf.nComp(),ngrow,Fab_allocate);

//
// Copy the valid region to the temporary MF
//
    MultiFab::Copy(tmp,mf,0,0,mf.nComp(),ngrow);

    for (MultiFabIterator mfiDst(tmp); mfiDst.isValid(); ++mfiDst)
    {
        Box dstBox = mfiDst.fabbox();

        for (MultiFabIterator mfiSrc(tmp); mfiSrc.isValid(); ++mfiSrc)
        {
            Box srcBox = mfiSrc.validbox();

            if (mfiDst.index() != mfiSrc.index() &&
                dstBox.intersects(srcBox))
            {
                Box intersection = srcBox & dstBox;

                mfiDst().copy(mfiSrc(), intersection, 0,
                                        intersection, 0, mf.nComp());
            }
        }
    }

    for (MultiFabIterator mfi(tmp); mfi.isValid(); ++mfi)
    {
        int indx = mfi.index();

        mfi().minus(mf[indx]);
    }

//
// ASCII output
//
    if (ascii)
    {
        for (MultiFabIterator mfi(tmp); mfi.isValid(); ++mfi)
        {
            cout << "FAB: " << mfi.index() << endl;
            cout << mfi() << endl;
        }
        return true;
    }

    return ArrayViewMultiFab(&tmp);
}
