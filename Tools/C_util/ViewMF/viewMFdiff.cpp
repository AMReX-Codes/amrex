#ifdef BL_ARCH_CRAY
#ifdef BL_USE_DOUBLE
#error "DOUBLE PRECISION NOT ALLOWED ON CRAY"
#endif
#endif

#ifndef        WIN32
#include <unistd.h>
#endif

#include "MultiFab.H"
#include "ArrayView.H"
#include "ParmParse.H"
#include "Utility.H"
#include "ParallelDescriptor.H"
#include "TV_TempWrite.H"
#include "MFNorm.H"

#ifdef BL_USE_NEW_HFILES
#include <new>
using std::setprecision;
#include <iostream>
#ifndef WIN32
using std::set_new_handler;
#endif
#else
#include <new.h>
#include <iostream.h>
#endif


static
void 
PrintUsage(int argc, char *argv[])
{
    cout << "Usage: " << endl;
    cout << argv[0] << endl;
    cout << "     mfab0 = MultiFab 0" << endl;
    cout << "     mfab1 = MultiFab 1" << endl;
    cout << "     comp0 = Starting component for MultiFab 0" << endl
         << "             (default is comp0=0)" << endl;
    cout << "     comp1 = Starting component for MultiFab 1" << endl
         << "             (default is comp1=0)" << endl;
    cout << "     ncomp = Number of components" << endl
         << "             (default is ncomp=mfab0.nComp())" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -validbox : Confine the difference to the valid region" << endl
         << "              of the FABs.  Difference values are set to" << endl
         << "              zero on the ghostcells.  Default behaviour" << endl
         << "              is to compare the entire FAB." << endl;
    exit(0);
}


int main (int   argc,
	  char* argv[])
{
    //
    // Make sure to catch new failures.
    //
#ifndef WIN32
    set_new_handler(Utility::OutOfMemory);
#endif

    ParallelDescriptor::StartParallel(&argc, &argv);
    
//
//  Parse the command line
//
    if (argc == 1)
        PrintUsage(argc, argv);

    ParmParse pp(argc-1,argv+1);

    if (pp.contains("help"))
        PrintUsage(argc, argv);

    aString iFile0, iFile1, diffFile;

    pp.query("mfab0", iFile0);
    if (iFile0.isNull())
        BoxLib::Abort("You must specify `mfab0'");

    pp.query("mfab1", iFile1);
    if (iFile1.isNull())
        BoxLib::Abort("You must specify `mfab1'");

    int comp0 = 0;
    pp.query("comp0", comp0);

    int comp1 = 0;
    pp.query("comp1", comp1);

    int nComp = -1;
    pp.query("ncomp", nComp);

    bool validBoxOnly;
    validBoxOnly = pp.contains("validbox");


//
//  Actually Calculate the Difference
//
    MultiFab mf0, mf1;
    readMF(mf0,iFile0.c_str());
    readMF(mf1,iFile1.c_str());

    if (mf0.boxArray() != mf1.boxArray())
        BoxLib::Abort("BoxArray's incompatible");
    
    if (mf0.nGrow()  != mf1.nGrow())
        BoxLib::Abort("nGrow's incompatible");

    if (nComp == -1) {
        if (mf0.nComp() != mf1.nComp())
            BoxLib::Abort("You must specify `ncomp' if (mfab0.nComp() != mfab1.nComp())");

        nComp = mf0.nComp();
    }

    if (mf0.nComp() < comp0 + nComp  || mf1.nComp() < comp1 + nComp)
    {
	cerr << "nComp's incompatible" << endl;
        cerr << "(need,have): (" << comp0 + nComp << "," << mf0.nComp() 
             << "), (" << comp1 + nComp << "," << mf1.nComp() << ")" << endl;
	return 0;
    }

//
//
//
    BoxArray compBoxes = mf0.boxArray();
    if (!validBoxOnly)
        compBoxes.grow(mf0.nGrow());

    MultiFab diff(mf0.boxArray(), nComp, mf0.nGrow(), Fab_allocate);
    for (MultiFabIterator diffmfi(diff); diffmfi.isValid(); ++diffmfi)
    {
	DependentMultiFabIterator mf0mfi(diffmfi, mf0);
	DependentMultiFabIterator mf1mfi(diffmfi, mf1);

        int mfiIndx = diffmfi.index();
	const Box compBox = compBoxes[mfiIndx];

	diffmfi().setVal(0.0);
	diffmfi().copy(mf0mfi(),compBox,comp0,compBox,0,nComp);
	diffmfi().minus(mf1mfi(),compBox,compBox,comp1,0,nComp);
    }
    Real norm0 = MFNorm(diff, compBoxes, 0, 0, nComp);
    Real norm1 = MFNorm(diff, compBoxes, 1, 0, nComp);
    Real norm2 = MFNorm(diff, compBoxes, 2, 0, nComp);

    if(ParallelDescriptor::IOProcessor())
    {
	cout << "Norms of diff (0,1,2): "
	     << norm0 << ", " << norm1 << ", " << norm2 << endl;
    }
    
    if (!diffFile.isNull())
    {
	writeMF(&diff,diffFile.c_str());
	return 1;

    } else {

	if (norm0 == 0 && norm1 == 0 && norm2 == 0)
	{
	    cout << "MultiFabs equal!" << endl;
	    return 1;

	} else {
	    
	    return ArrayViewMultiFab(&diff);
	}
    }
}

