//BL_COPYRIGHT_NOTICE

//
// $Id: tFB.cpp,v 1.2 1999-03-12 02:55:43 lijewski Exp $
//
// A test program for FillBoundary().
//

#include <Utility.H>
#include <MultiFab.H>

#ifdef BL_USE_NEW_HFILES
#include <new>
using std::setprecision;
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
    // Do loop filling both components.
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

    Box bx(IntVect::TheZeroVector(),IntVect::TheUnitVector());

    bl.append(bx);

    bx.shift(0,2);

    bl.append(bx);

    BoxArray ba(bl);
    
    MultiFab junk(ba,2,1);

    DoIt(junk);

    MultiFab junky(ba,2,1);

    DoIt(junky);
    
    ParallelDescriptor::EndParallel();

    return 0;
}
