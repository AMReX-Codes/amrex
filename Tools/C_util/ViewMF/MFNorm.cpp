#ifdef BL_ARCH_CRAY
#ifdef BL_USE_DOUBLE
#error "DOUBLE PRECISION NOT ALLOWED ON CRAY"
#endif
#endif

#ifndef        WIN32
#include <unistd.h>
#endif

#include "MultiFab.H"
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


//
// What's the slowest way I can think of to compute all the norms??
//
Real
MFNorm (const MultiFab& mfab, 
        const BoxArray& normBoxes, 
        const int       exponent,
        const int       srcComp,
        const int       numComp)
{
    // 
    // Check that the FabBoxes actually contain the normBoxes
    //
    BoxArray fabBoxes = mfab.boxArray();
    fabBoxes.grow(mfab.nGrow());
    assert (fabBoxes.contains(normBoxes));
    
    //
    // Get a copy of the multifab
    //
    MultiFab mftmp(normBoxes, numComp, 0);

    for (MultiFabIterator mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
    {
        DependentMultiFabIterator srcMfi(mftmpmfi, mfab);

        const Box& cpBox = normBoxes[mftmpmfi.index()];
        mftmpmfi().copy(srcMfi(), cpBox, srcComp, cpBox, 0, numComp);
    }


    //
    // Actually Calculate the Norms
    //
    Real myNorm = 0;
    if ( exponent == 0 )
    {
        for ( MultiFabIterator mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
        {
            const Box& cpBox = normBoxes[mftmpmfi.index()];

            mftmpmfi().abs(cpBox, 0, numComp);
 
            myNorm = Max(myNorm, mftmpmfi().norm(0, 0, numComp));
        }
	ParallelDescriptor::ReduceRealMax(myNorm);

    } else if ( exponent == 1 )
    {
        for ( MultiFabIterator mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
        {
            const Box& cpBox = normBoxes[mftmpmfi.index()];

            mftmpmfi().abs(cpBox, 0, numComp);

            myNorm += mftmpmfi().norm(1, 0, numComp);
        }
	ParallelDescriptor::ReduceRealSum(myNorm);

    } else if ( exponent == 2 )
    {
        for ( MultiFabIterator mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
        {
            const Box& cpBox = normBoxes[mftmpmfi.index()];

            mftmpmfi().abs(cpBox, 0, numComp);

            myNorm += pow(mftmpmfi().norm(2, 0, numComp), 2);
        }
	ParallelDescriptor::ReduceRealSum(myNorm);
        myNorm = sqrt( myNorm );

    } else {

        BoxLib::Error("Invalid exponent to norm function");
    }
    
    return myNorm;
}

