//BL_COPYRIGHT_NOTICE

//
// $Id: tFB.cpp,v 1.3 1999-03-12 03:27:07 lijewski Exp $
//
// A test program for FillBoundary().
//

#include <Utility.H>
#include <MultiFab.H>

#ifdef BL_USE_NEW_HFILES
#include <new>
#ifndef WIN32
using std::set_new_handler;
#endif
#else
#include <new.h>
#endif

static
void
DumpFirstFab (const MultiFab& mf)
{
    assert(mf.length() > 1);

    if (ParallelDescriptor::IOProcessor())
    {
        assert(mf.DistributionMap()[0] == ParallelDescriptor::MyProc());

        cout << mf[0] << endl;
    }
}

static
void
DoIt (MultiFab& mf)
{
    assert(mf.nComp() >= 2);

    mf.setVal(0);
    //
    // Do loop just filling first component.
    //
    for (int i = 0; i < 10; i++)
    {
        for(MultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
        {
            mfi().setVal(mfi.index()+1,0);
        }

        mf.FillBoundary(0,1);

        DumpFirstFab(mf);
    }
    //
    // Do loop just filling second component.
    //
    for (int i = 0; i < 10; i++)
    {
        for(MultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
        {
            mfi().setVal(mfi.index()+10,1);
        }

        mf.FillBoundary(1,1);

        DumpFirstFab(mf);
    }
    //
    // Do loop filling all components.
    //
    for (int i = 0; i < 10; i++)
    {
        for(MultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
        {
            mfi().setVal(mfi.index()+100);
        }

        mf.FillBoundary();

        DumpFirstFab(mf);
    }
}

int
main (int argc, char** argv)
{
#ifndef WIN32
    set_new_handler(Utility::OutOfMemory);
#endif

    ParallelDescriptor::StartParallel(&argc,&argv);

    BoxList bl;

    Box bx1(IntVect::TheZeroVector(),IntVect::TheUnitVector());
    Box bx2 = bx1;
    bx2.shift(0,2);

    bl.append(bx1); bl.append(bx2);

    BoxArray ba(bl);
    
    MultiFab junk(ba,2,1), junky(ba,2,1);

    DoIt(junk); DoIt(junky); DoIt(junk); DoIt(junky);

    ParallelDescriptor::EndParallel();

    return 0;
}
