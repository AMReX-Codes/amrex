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

#include <new>
using std::setprecision;
#include <iostream>
#ifndef WIN32
using std::set_new_handler;
#endif

static
void 
PrintUsage(int argc, char *argv[])
{
    cout << "Usage: " << endl;
    cout << argv[0] << " iFile0=mf0 iFile1=mf1 [options] \n\tOptions:" << endl;
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
    BoxLib::Initialize(argc,argv);

    ParmParse pp;
    
    if (pp.contains("help"))
        PrintUsage(argc, argv);
    
    std::string iFile0;
    std::string iFile1;

    pp.query("iFile0", iFile0);
    pp.query("iFile1", iFile1);

    if (iFile0.empty() || iFile1.empty())
        PrintUsage(argc, argv);

    int comp0 = 0;
    pp.query("comp0", comp0);

    int comp1 = 0;
    pp.query("comp1", comp1);

    int nComp = -1;
    pp.query("ncomp", nComp);
    
    int ngrow = 0;
    pp.query("ngrow",ngrow);
    BL_ASSERT(ngrow>=0);

    std::string outfile;
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
        for (int i=0; i<mf0.boxArray().size(); ++i)
            common_bl.join(BoxList(BoxLib::intersect(mf1.boxArray(), mf0.boxArray()[i])));
        compBoxes = BoxArray(common_bl);
    }
    
    if (ngrow != std::min(ngrow,mf0.nGrow()))
    {
        BoxLib::Warning("Shrinking ngrow to that available in mfab0");
        ngrow = mf0.nGrow();
    }

    if (ngrow != std::min(ngrow,mf1.nGrow()))
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
    PArray<FArrayBox> fabs(compBoxes.size(),PArrayManage);
    for (int i=0; i<fabs.size(); ++i)
    {
        fabs.set(i,new FArrayBox(compBoxes[i],nComp));
        fabs[i].setVal(0.0);
    
        for (MFIter mf0_mfi(mf0); mf0_mfi.isValid(); ++mf0_mfi)
        {
            const Box box = BoxLib::grow(mf0_mfi.validbox(),ngrow) & fabs[i].box();
            if (box.ok())
                fabs[i].copy(mf0[mf0_mfi],box,comp0,box,0,nComp);
        }
        for (MFIter mf1_mfi(mf1); mf1_mfi.isValid(); ++mf1_mfi)
        {
            const Box box = BoxLib::grow(mf1_mfi.validbox(),ngrow) & fabs[i].box();
            if (box.ok())
                fabs[i].minus(mf1[mf1_mfi],box,box,comp1,0,nComp);
        }
    }

    //
    // Get the result into a viewable MultiFab
    //
    MultiFab diffmfab(compBoxes,nComp,ngrow,Fab_allocate);
    for (MFIter mfi(diffmfab); mfi.isValid(); ++mfi)
        for (int i=0; i<fabs.size(); ++i)
            diffmfab[mfi].copy(fabs[i]);

    Real norm0 = MFNorm(diffmfab, 0, 0, nComp, ngrow);
    Real norm1 = MFNorm(diffmfab, 1, 0, nComp, ngrow);
    Real norm2 = MFNorm(diffmfab, 2, 0, nComp, ngrow);

    if(ParallelDescriptor::IOProcessor())
    {
	cout << "Norms of diff (0,1,2): "
	     << norm0 << ", " << norm1 << ", " << norm2 << endl;
    }
    
    if (!outfile.empty())
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

    BoxLib::Finalize();
}

