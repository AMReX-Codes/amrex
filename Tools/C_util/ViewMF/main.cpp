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

#include <new>
using std::setprecision;
#include <iostream>
#ifndef WIN32
using std::set_new_handler;
#endif

const int NPROCS = 1;

REAL
norm ( const MultiFab& mfab,
       int             exponent );

// Set view to 1 or 2, and define or undefine the BUILD_TESTDATA macro
//  to generate some test data

#define BUILD_TESTDATA
#undef BUILD_TESTDATA

#define VIEW 2

#if (VIEW == 1)

void 
print_usage(int, char *argv[])
{
    cerr << "usage: ";
    cerr << argv[0] << " MultiFabFile" << endl;
    exit(0);
}


#ifdef BUILD_TESTDATA

int
main (int   argc,
      char* argv[])
{
    //
    // Make sure to catch new failures.
    //
#ifndef WIN32
    set_new_handler(Utility::OutOfMemory);
#endif

    // Parse command line
    ParmParse pp(argc-1,argv+1,NULL,NULL); 

    int nprocs = NPROCS; pp.query("nprocs", nprocs);
#ifndef BL_USE_BSP
    if (nprocs > 1)
    {
      cerr << "Error in main:  multiple processors specified with "
           << "code compiled without a parallel library.\n";
      exit(-1);
    }
#endif
    ParallelDescriptor::StartParallel();

    Box box0(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(3,3,3)));
    Box box1(IntVect(D_DECL(4,4,4)),IntVect(D_DECL(7,7,7)));
    BoxList bl;
    bl.append(box0);
    bl.append(box1);
    BoxArray ba(bl);

    int nGrow = 1;
    int sComp = 0;
    int nComp = 5;

    MultiFab tmp0(ba, nComp, nGrow, Fab_allocate);
    tmp0.setVal(5);
    for (MultiFabIterator tmp0mfi(tmp0); tmp0mfi.isValid(); ++tmp0mfi)
    {
	tmp0mfi().setVal(tmp0mfi.index()+1,tmp0mfi.validbox(),sComp,nComp);
    }
    writeMF(&tmp0,"tmp0");
    
    ParallelDescriptor::EndParallel();
}
    


#else  /* BUILD_TESTDATA */

int
main (int   argc,
      char* argv[])
{
    //
    // Make sure to catch new failures.
    //
#ifndef WIN32
    set_new_handler(Utility::OutOfMemory);
#endif

    if (argc != 2) print_usage(argc,argv);

    int nprocs = NPROCS;
#ifndef BL_USE_BSP
    if (nprocs > 1)
    {
      cerr << "Error in main:  multiple processors specified with "
           << "code compiled without a parallel library.\n";
      exit(-1);
    }
#endif
    ParallelDescriptor::StartParallel();

    MultiFab mf;
    readMF(mf,argv[1]);
    
    return ArrayViewMultiFab(&mf);
}

#endif /* BUILD_TESTDATA */

#elif (VIEW == 2)

void 
print_usage(int, char *argv[])
{
    cerr << "usage: ";
    cerr << argv[0] << " MF1 MF2 comp1 comp2 nComp [diffOutfile]" << endl;
    exit(0);
}

#ifdef BUILD_TESTDATA

int
main (int   argc,
      char* argv[])
{
    //
    // Make sure to catch new failures.
    //
#ifndef WIN32
    set_new_handler(Utility::OutOfMemory);
#endif

    // Parse command line
    ParmParse pp(argc-1,argv+1,NULL,NULL); 

    int nprocs = NPROCS; pp.query("nprocs", nprocs);
#ifndef BL_USE_BSP
    if (nprocs > 1)
    {
      cerr << "Error in main:  multiple processors specified with "
           << "code compiled without a parallel library.\n";
      exit(-1);
    }
#endif
    ParallelDescriptor::StartParallel();

    Box box0(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(3,3,3)));
    Box box1(IntVect(D_DECL(4,4,4)),IntVect(D_DECL(7,7,7)));
    BoxList bl;
    bl.append(box0);
    bl.append(box1);
    BoxArray ba(bl);

    int nGrow = 1;
    int sComp = 0;
    int nComp = 5;

    MultiFab tmp0(ba, nComp, nGrow, Fab_allocate);
    tmp0.setVal(5);
    for (MultiFabIterator tmp0mfi(tmp0); tmp0mfi.isValid(); ++tmp0mfi)
    {
	tmp0mfi().setVal(tmp0mfi.index()+1,tmp0mfi.validbox(),sComp,nComp);
    }

    writeMF(&tmp0,"tmp0");

    MultiFab tmp1(ba, nComp, nGrow, Fab_allocate);
    tmp1.setVal(7);
    for (MultiFabIterator tmp1mfi(tmp1); tmp1mfi.isValid(); ++tmp1mfi)
    {
	tmp1mfi().setVal(tmp1mfi.index()+3,tmp1mfi.validbox(),sComp,nComp);
    }

    writeMF(&tmp1,"tmp1");
    
    ParallelDescriptor::EndParallel();
}

#else  /* BUILD_TESTDATA */

int main (int   argc,
	  char* argv[])
{
    //
    // Make sure to catch new failures.
    //
#ifndef WIN32
    set_new_handler(Utility::OutOfMemory);
#endif
    
    if (argc < 6 || argc > 7) print_usage(argc,argv);

    int nprocs = NPROCS;
#ifndef BL_USE_BSP
    if (nprocs > 1)
    {
      cerr << "Error in main:  multiple processors specified with "
           << "code compiled without a parallel library.\n";
      exit(-1);
    }
#endif
    ParallelDescriptor::StartParallel();

    MultiFab mf0,mf1;
    readMF(mf0,argv[1]);
    readMF(mf1,argv[2]);

    if (mf0.boxArray() != mf1.boxArray())
    {
	cerr << "BoxArray's incompatible" << endl;
	return 0;
    }
    
    if (mf0.nGrow()  != mf1.nGrow())
    {
	cerr << "nGrow's incompatible" << endl;
	return 0;
    }

    int comp0 = atoi(argv[3]);
    int comp1 = atoi(argv[4]);
    int nComp = atoi(argv[5]);
    if (mf0.nComp() < comp0 + nComp  || mf1.nComp() < comp1 + nComp)
    {
	cerr << "nComp's incompatible" << endl;
        cerr << "(need,have): (" << comp0 + nComp << "," << mf0.nComp() 
             << "), (" << comp1 + nComp << "," << mf1.nComp() << ")" << endl;
	return 0;
    }

    MultiFab diff(mf0.boxArray(), nComp, mf0.nGrow(), Fab_allocate);
    for (MultiFabIterator diffmfi(diff); diffmfi.isValid(); ++diffmfi)
    {
	DependentMultiFabIterator mf0mfi(diffmfi, mf0);
	DependentMultiFabIterator mf1mfi(diffmfi, mf1);

	const Box& validBox = mf0mfi.validbox();
	diffmfi().setVal(0.0);
	diffmfi().copy(mf0mfi(),validBox,comp0,validBox,0,nComp);
	diffmfi().minus(mf1mfi(),validBox,validBox,comp1,0,nComp);
    }
    
    REAL norm0 = norm(diff,0);
    REAL norm1 = norm(diff,1);
    REAL norm2 = norm(diff,2);

    if(ParallelDescriptor::IOProcessor())
    {
	cout << "Norms of diff (0,1,2): "
	     << norm0 << ", " << norm1 << ", " << norm2 << endl;
    }
    
    if (argc == 7)
    {
	writeMF(&diff,argv[6]);
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

#endif /* BUILD_TESTDATA */

#endif

// What's the slowest way I can think of to compute all the norms??
REAL
norm ( const MultiFab& mfab,
       int             exponent )
{
    int srcComp = 0;
    int numComp = mfab.nComp();
    int numGrow = 0;
    const BoxArray& grids = mfab.boxArray();
    
    // Get a copy of the multifab, zero covered locations
    MultiFab mftmp( grids, numComp, numGrow, Fab_allocate );
    mftmp.copy( mfab );

    Real myNorm = 0;
    if ( exponent == 0 )
    {
        for ( MultiFabIterator mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
        {
            mftmpmfi().abs(mftmpmfi.validbox());
        }
        myNorm = mftmp.max( srcComp );
	ParallelDescriptor::ReduceRealMax(myNorm);

    } else if ( exponent == 1 )
    {
        for ( MultiFabIterator mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
        {
            mftmpmfi().abs();
            myNorm += mftmpmfi().sum( srcComp, numComp );
        }
	ParallelDescriptor::ReduceRealSum(myNorm);

    } else if ( exponent == 2 )
    {
        for ( MultiFabIterator mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
        {
            const Box& grid = mftmpmfi.validbox();
            for ( IntVect idx = grid.smallEnd(); idx <= grid.bigEnd(); grid.next(idx) )
            {
                mftmpmfi()(idx) = mftmpmfi()(idx) * mftmpmfi()(idx);
            }
        }
        int sum_exp = 1;
        for ( MultiFabIterator mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
        {
	    const Box& bx = mftmpmfi.validbox();
            myNorm += mftmpmfi().norm( bx, sum_exp, srcComp, numComp );
        }
	ParallelDescriptor::ReduceRealSum(myNorm);
        myNorm = sqrt( myNorm );

    } else {

        BoxLib::Error("Invalid exponent to norm function");
    }
    
    return myNorm;
}

