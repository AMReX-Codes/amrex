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
    
    aString outfile;
    pp.query("outfile",outfile);

    int ngrow = 0;

    int lvl = 0;
    pp.query("lvl",lvl);

    int dir = 2;
    pp.query("dir",dir);

    IntVect edge0Small,  edge0Big, edge1Small, edge1Big;
    if (lvl == 0)
    {
       if (dir == 0)
        {
            edge0Small = IntVect(D_DECL( 0,  0,-1));
            edge0Big   = IntVect(D_DECL( 0,119,49));
            edge1Small = IntVect(D_DECL(72,  0,-1));
            edge1Big   = IntVect(D_DECL(72,119,49));
        }
       else if (dir == 2)
        {
            edge0Small = IntVect(D_DECL(-1,  0, 0));
            edge0Big   = IntVect(D_DECL(72,119, 0));
            edge1Small = IntVect(D_DECL(-1,  0,48));
            edge1Big   = IntVect(D_DECL(72,119,48));
        }
    }
    else if (lvl == 1)
    {
       if (dir == 0)
        {
            edge0Small = IntVect(D_DECL(  0,  0, -1));
            edge0Big   = IntVect(D_DECL(  0,479,192));
            edge1Small = IntVect(D_DECL(288,  0, -1));
            edge1Big   = IntVect(D_DECL(288,479,192));
        }
       else if (dir == 2)
        {
            edge0Small = IntVect(D_DECL( 4,232,  0));
            edge0Big   = IntVect(D_DECL(31,243,  0));
            edge1Small = IntVect(D_DECL( 4,232,192));
            edge1Big   = IntVect(D_DECL(31,243,192));
/*
            edge0Small = IntVect(D_DECL( -1,  0,  0));
            edge0Big   = IntVect(D_DECL(288,479,  0));
            edge1Small = IntVect(D_DECL( -1,  0,192));
            edge1Big   = IntVect(D_DECL(288,479,192));
*/
        }
    }

//
//  Actually Calculate the Difference
//
    MultiFab mf0, mf1;
    readMF(mf0,iFile0.c_str());
    readMF(mf1,iFile1.c_str());

    Box pedge0(edge0Small, edge0Big, mf0.box(0).ixType());
    Box pedge1(edge1Small, edge1Big, mf1.box(0).ixType());

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
cout << pedge0 << " " << pedge1 << endl;
    BoxArray pedge0BA(&pedge0, 1);
    BoxArray pedge1BA(&pedge1, 1);

    MultiFab datEdge0(pedge0BA, nComp, ngrow);
    MultiFab datEdge1(pedge1BA, nComp, ngrow);
    datEdge0.setVal(1.0e30);
    datEdge1.setVal(1.0e30);
    
    for (MultiFabIterator mf0_mfi(mf0); mf0_mfi.isValid(); ++mf0_mfi)
    {
        const Box box = ::grow(mf0_mfi.validbox(),mf0.nGrow()) & datEdge0[0].box();
        if (box.ok())
            datEdge0[0].copy(mf0_mfi(),box,comp0,box,0,nComp);
    }

    for (MultiFabIterator mf1_mfi(mf1); mf1_mfi.isValid(); ++mf1_mfi)
    {
        const Box box = ::grow(mf1_mfi.validbox(),mf1.nGrow()) & datEdge1[0].box();
        if (box.ok())
            datEdge1[0].copy(mf1_mfi(),box,comp1,box,0,nComp);
    }


    //
    // Get the result into a viewable MultiFab
    //
    MultiFab diffmfab(datEdge0.boxArray(),nComp,ngrow,Fab_allocate);
    diffmfab.setVal(0.0);

    FArrayBox shftedDatEdge1(pedge0, nComp);
    shftedDatEdge1.copy(datEdge1[0], pedge1, 0, pedge0, 0, nComp);

    diffmfab.copy(datEdge0);
    diffmfab[0].minus(shftedDatEdge1);

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

