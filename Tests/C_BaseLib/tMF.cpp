
//
// $Id: tMF.cpp,v 1.2 2004-03-05 17:49:38 car Exp $
//
// A test program for MultiFab.
//

#include <winstd.H>

#include <Utility.H>
#include <MultiFab.H>

int
main (int argc, char** argv)
{
    BoxLib::Initialize(argc, argv);

    Box bx1(IntVect::TheZeroVector(),IntVect::TheUnitVector());
    Box bx2 = bx1;
    bx2.shift(0,2);

    BoxArray ba(2);
    ba.set(0, bx1);
    ba.set(1, bx2);

    MultiFab junk(ba,2,1);

    BoxLib::Finalize();

    return 0;
}
