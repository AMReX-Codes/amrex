//BL_COPYRIGHT_NOTICE

//
// $Id: tFB.cpp,v 1.1 1998-09-30 17:07:39 lijewski Exp $
//
// A test program for FillBoundary().
//

#include <Utility.H>
#include <MultiFab.H>
#include <new.h>

static
void
dumpMF (const MultiFab& mf,
        const char*     str)
{
    if (ParallelDescriptor::IOProcessor())
    {
	cout << str << endl;
	ConstMultiFabIterator it(mf);
	for (;it.isValid(); ++it)
	{
	    cout << it() << endl;
	}
    }
    if (!ParallelDescriptor::IOProcessor())
    {
	ConstMultiFabIterator it(mf);
	for (;it.isValid(); ++it)
	{
	    cout << it() << endl;
	}
    }
    ParallelDescriptor::Barrier();
}

int
main (int argc, char** argv)
{
    int nprocs = 1;

    ParallelDescriptor::StartParallel(nprocs,&argc,&argv);

    BoxList bl;
    bl.append(Box(IntVect(0,0),IntVect(1,1)));
    bl.append(Box(IntVect(0,2),IntVect(1,3)));
    BoxArray bs(bl);
    
    MultiFab junk(bs, 1, 1, Fab_allocate);

    for(MultiFabIterator mfi(junk); mfi.isValid(); ++mfi)
	mfi().setVal(mfi.index()+1);
    junk.FillBoundary();
    dumpMF(junk,"after 1 time...");
    
    for(MultiFabIterator mfi(junk); mfi.isValid(); ++mfi)
	mfi().setVal(mfi.index()+5);
    junk.FillBoundary();
    dumpMF(junk,"after 2 times...");
    
    ParallelDescriptor::EndParallel();
}
