
//
// $Id: tMF.cpp,v 1.1 2004-03-04 18:09:01 car Exp $
//
// A test program for FillBoundary().
//

#include <winstd.H>

#include <Utility.H>
#include <MultiFab.H>

#include <new>
#ifndef WIN32
using std::set_new_handler;
#endif

int
main (int argc, char** argv)
{
#ifndef WIN32
    set_new_handler(BoxLib::OutOfMemory);
#endif
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
