
#include <unistd.h>

#include <AMReX_MultiFab.H>
#include <MFNorm.H>

#include <iostream>
//
// What's the slowest way I can think of to compute all the norms??
//
Real
MFNorm (const MultiFab& mfab, 
        const int       exponent,
        const int       srcComp,
        const int       numComp,
        const int       numGrow)
{
    BL_ASSERT (numGrow <= mfab.nGrow());
    BoxArray boxes = mfab.boxArray();
    boxes.grow(numGrow);
    //
    // Get a copy of the multifab
    //
    MultiFab mftmp(mfab.boxArray(), numComp, 0);
    MultiFab::Copy(mftmp,mfab,srcComp,0,numComp,numGrow);
    //
    // Calculate the Norms
    //
    Real myNorm = 0;
    if ( exponent == 0 )
    {
        for ( MFIter mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
        {
            mftmp[mftmpmfi].abs(boxes[mftmpmfi.index()], 0, numComp);
            myNorm = std::max(myNorm, mftmp[mftmpmfi].norm(0, 0, numComp));
        }
	ParallelDescriptor::ReduceRealMax(myNorm);

    } else if ( exponent == 1 )
    {
        for ( MFIter mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
        {
            mftmp[mftmpmfi].abs(boxes[mftmpmfi.index()], 0, numComp);

            myNorm += mftmp[mftmpmfi].norm(1, 0, numComp);
        }
	ParallelDescriptor::ReduceRealSum(myNorm);

    } else if ( exponent == 2 )
    {
        for ( MFIter mftmpmfi(mftmp); mftmpmfi.isValid(); ++mftmpmfi)
        {
            mftmp[mftmpmfi].abs(boxes[mftmpmfi.index()], 0, numComp);

            myNorm += pow(mftmp[mftmpmfi].norm(2, 0, numComp), 2);
        }
	ParallelDescriptor::ReduceRealSum(myNorm);
        myNorm = sqrt( myNorm );

    } else {

        amrex::Error("Invalid exponent to norm function");
    }
    
    return myNorm;
}

