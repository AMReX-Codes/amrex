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
    cout << argv[0] << " infile0 infile1 [options] \n\tOptions:" << endl;
    cout << "\t   comp0 = Starting component for MultiFab in infile0" << endl
         << "\t           (default is comp0=0)" << endl;
    cout << "\t   comp1 = Starting component for MultiFab in infile1" << endl
         << "\t           (default is comp1=0)" << endl;
    cout << "\t   ncomp = Number of components" << endl
         << "\t           (default is ncomp=mfab0.nComp())" << endl;
    cout << "\t outfile = Name of file to dump result of diff" << endl;
    cout << "\t   ngrow = Number of grow cells to include in result" << endl;
    cout << "\t           (default is ngrow=0)" << endl;
    cout << endl;
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
    if (argc < 3)
        PrintUsage(argc,argv);

    if (argv[1][0] == '-')
    {
        cerr << "input file must be first argument\n";
        PrintUsage(argc, argv);
    }

    ParmParse pp(argc-3,argv+3);
    
    if (pp.contains("help"))
        PrintUsage(argc, argv);
    
    const aString iFile0 = argv[1];
    const aString iFile1 = argv[2];

    int comp0 = 0;
    pp.query("comp0", comp0);

    int comp1 = 0;
    pp.query("comp1", comp1);

    int nComp = -1;
    pp.query("ncomp", nComp);
    
    int ngrow = 0;
    pp.query("ngrow",ngrow);
    assert(ngrow>=0);

    aString outfile;
    pp.query("outfile",outfile);

//
//  Actually Calculate the Difference
//
    MultiFab mf0, mf1;
    readMF(mf0,iFile0.c_str());
    readMF(mf1,iFile1.c_str());

    BoxArray compBoxes = mf0.boxArray();

    if (mf0.boxArray() != mf1.boxArray())
    {
        //
        // For this, assume no grow cells
        //
        BoxList common_bl;
        for (int i=0; i<mf0.boxArray().length(); ++i)
            common_bl.join(BoxList(::intersect(mf1.boxArray(), mf0.boxArray()[i])));
        compBoxes = BoxArray(common_bl);
    }
    
    if (ngrow != Min(ngrow,mf0.nGrow()))
    {
        BoxLib::Warning("Shrinking ngrow to that available in mfab0");
        ngrow = mf0.nGrow();
    }

    if (ngrow != Min(ngrow,mf1.nGrow()))
    {
        BoxLib::Warning("Shrinking ngrow to that available in mfab1");
        ngrow = mf1.nGrow();
    }

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
    // Result may have different processor mapping than src mfabs
    // I think that means you got to things in the following way
    //
    PArray<FArrayBox> fabs(compBoxes.length(),PArrayManage);
    for (int i=0; i<fabs.length(); ++i)
    {
        fabs.set(i,new FArrayBox(compBoxes[i],nComp));
        fabs[i].setVal(0.0);
    
        for (MultiFabIterator mf0_mfi(mf0); mf0_mfi.isValid(); ++mf0_mfi)
        {
            const Box box = ::grow(mf0_mfi.validbox(),ngrow) & fabs[i].box();
            if (box.ok())
                fabs[i].copy(mf0_mfi(),box,comp0,box,0,nComp);
        }
        for (MultiFabIterator mf1_mfi(mf1); mf1_mfi.isValid(); ++mf1_mfi)
        {
            const Box box = ::grow(mf1_mfi.validbox(),ngrow) & fabs[i].box();
            if (box.ok())
                fabs[i].minus(mf1_mfi(),box,box,comp1,0,nComp);
        }
    }

    //
    // Get the result into a viewable MultiFab
    //
    MultiFab diffmfab(compBoxes,nComp,ngrow,Fab_allocate);
    for (MultiFabIterator mfi(diffmfab); mfi.isValid(); ++mfi)
        for (int i=0; i<fabs.length(); ++i)
            mfi().copy(fabs[i]);

    Real norm0 = MFNorm(diffmfab, 0, 0, nComp, ngrow);
    Real norm1 = MFNorm(diffmfab, 1, 0, nComp, ngrow);
    Real norm2 = MFNorm(diffmfab, 2, 0, nComp, ngrow);

    if(ParallelDescriptor::IOProcessor())
    {
	cout << "Norms of diff (0,1,2): "
	     << norm0 << ", " << norm1 << ", " << norm2 << endl;
    }
    
    if (!outfile.isNull())
    {
	writeMF(&diffmfab,outfile.c_str());
	return 1;

    } else {

	if (norm0 == 0 && norm1 == 0 && norm2 == 0)
	{
	    cout << "MultiFabs equal!" << endl;
	    return 1;

	} else {
	    
	    return ArrayViewMultiFab(&diffmfab);
	}
    }
}

