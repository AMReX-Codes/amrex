//BL_COPYRIGHT_NOTICE

//
// $Id: MultiFab.cpp,v 1.3 1997-09-17 17:46:31 lijewski Exp $
//

#include <iostream.h>
#include <iomanip.h>

#include <Assert.H>
#include <Misc.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>
#include <Utility.H>

#if defined(BL_ARCH_IEEE)
#ifdef BL_USE_DOUBLE
    const Real INFINITY = 1.0e100;
#elif  BL_USE_FLOAT
    const Real INFINITY = 1.0e35;
#endif
#elif defined(BL_ARCH_CRAY)
    const Real INFINITY = 1.0e100;
#endif

//
// This isn't inlined as it's virtual.
//
MultiFab::~MultiFab() {}

ostream &
operator<< (ostream&        os,
            const MultiFab& mf)
{
    cerr << "Error:  MultiFab operator<< not implemented for parallel" << endl;
    ParallelDescriptor::Abort("MultiFab operator<<");

    os << "(MultiFab "
       << mf.length() <<  ' '
       << mf.nGrow()  << '\n';
    for(int i = 0; i < mf.length(); ++i)
        os << mf[i] << '\n';
    os << ")" << flush;

    if (os.fail())
        BoxLib::Error("operator<<(ostream&,MultiFab&) failed");

    return os;
}

ostream &
MultiFab::writeOn (ostream& os) const
{
    assert(boxarray.ready());

    if (ParallelDescriptor::IOProcessor())
    {
        os << n_comp << '\n';
        os << n_grow << '\n';
        boxarray.writeOn(os);
        os << 0 <<  '\n';
    }

    streampos filePosition;
    ParallelDescriptor::ShareVar(&filePosition, sizeof(streampos));
    ParallelDescriptor::Synchronize();

    int myproc = ParallelDescriptor::MyProc();
    int fabProc;
    for (int i = 0; i < length(); ++i)
    {
        fabProc = distributionMap[i];
        if (fabProc == myproc)
        {
            fabparray[i].writeOn(os);
            filePosition = os.tellp();
        }
        ParallelDescriptor::Broadcast(fabProc, &filePosition, &filePosition);
        os.seekp(filePosition);
    }

    ParallelDescriptor::Synchronize();
    ParallelDescriptor::UnshareVar(&filePosition);
    //
    // No need to check os here as itll be done in FArrayBox::writeOn().
    //
    return os;
}

istream &
MultiFab::readFrom (istream& is)
{
    assert(!boxarray.ready());

    is >> n_comp;
    while (is.get() != '\n')
        ;

    is >> n_grow;
    while (is.get() != '\n')
        ;

    boxarray.define(is);
    distributionMap.define(ParallelDescriptor::NProcs(), boxarray);

    int has_ba;
    is >> has_ba;
    if (has_ba)
    {
        int cw;
        is >> cw;
        //
	// Do not do anything with BoxAssoc.
        //
    }
    while (is.get() != '\n')
        ;

    int nbox = boxarray.length();
    fabparray.resize(nbox);

/*
    original code:

    for (int i = 0; i < nbox; i++)
    {
        FArrayBox* tmp = new FArrayBox;
        if (tmp == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);
        tmp->readFrom(is);
        fabparray.set(i,tmp);
    }
*/

    streampos filePosition;
    ParallelDescriptor::ShareVar(&filePosition, sizeof(streampos));
    ParallelDescriptor::Synchronize();

    int myproc = ParallelDescriptor::MyProc();
    int fabProc;
    for (int i = 0; i < nbox; i++)
    {
        fabProc = distributionMap[i];
        if (fabProc == myproc)
        {
            FArrayBox* tmp = new FArrayBox;
            if (tmp == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);

            tmp->readFrom(is);
            fabparray.set(i,tmp);
            filePosition = is.tellg();
        }
        ParallelDescriptor::Broadcast(fabProc, &filePosition, &filePosition);
        is.seekg(filePosition);
    }

    ParallelDescriptor::Synchronize();
    ParallelDescriptor::UnshareVar(&filePosition);

    if (is.fail())
        BoxLib::Error("MultiFab::readFrom(istream&) failed");

    return is;
}

void
MultiFab::probe (ostream& os,
                 IntVect& pt)
{
    Real  dat[20];
    int prec = os.precision(14);
    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        if (mfi.validbox().contains(pt))
        {
            int nv = mfi().nComp();

            assert(nv <= 20);

            mfi().getVal(dat,pt);
            os << "point " << pt << " in box " << mfi.validbox()
               << " data = ";
            for (int i = 0; i < nv; i++)
                os << "  " << setw(20) << dat[i];
            os << '\n';
        }
    }
    os.precision(prec);
    os.flush();

    if (os.fail())
        BoxLib::Error("MultiFab::probe(ostream&,IntVect&) failed");
}

Real
MultiFab::min (int comp,
               int nghost) const
{
    assert(nghost >= 0 && nghost <= n_grow);

    Real mn = INFINITY;

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        mn = Min(mn, mfi().min(b, comp));
    }

    ParallelDescriptor::ReduceRealMin(mn);

    return mn;
}

Real
MultiFab::min (const Box& region,
               int        comp,
               int        nghost) const
{
    assert(nghost >= 0 && nghost <= n_grow);

    Real mn = INFINITY;

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        b &= region;
        if(b.ok())
        {
            mn = Min(mn, mfi().min(b, comp));
        }
    }

    ParallelDescriptor::ReduceRealMin(mn);

    return mn;
}

Real
MultiFab::max (int comp,
               int nghost) const
{
    assert(nghost >= 0 && nghost <= n_grow);

    Real mn = -INFINITY;

    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        mn = Max(mn, mfi().max(b, comp));
    }

    ParallelDescriptor::ReduceRealMax(mn);

    return mn;
}

Real
MultiFab::max (const Box& region,
               int        comp,
               int        nghost) const
{
    assert(nghost >= 0 && nghost <= n_grow);

    Real mn = -INFINITY;

    int first = true;
    for (ConstMultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        b &= region;
        if (b.ok())
        {
            Real mg = mfi().max(b,comp);
            if(first)
            {
                mn = mg;
                first = false;
            }
            mn = Max(mn,mg);
        }
    }

    ParallelDescriptor::ReduceRealMax(mn);

    return mn;
}

void
MultiFab::minus (const MultiFab& mf,
                 int             strt_comp,
                 int             num_comp,
                 int             nghost)
{
    assert(boxarray == mf.boxarray);
    assert(strt_comp >= 0);
#ifndef NDEBUG
    int lst_comp = strt_comp + num_comp - 1;
#endif
    assert(lst_comp < n_comp && lst_comp < mf.n_comp);
    assert(nghost <= n_grow && nghost <= mf.n_grow);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi, mf);
        Box bx(mfi.validbox());
        bx.grow(nghost);
        mfi().minus(dmfi(), bx, strt_comp, strt_comp, num_comp);
    }
}

void
MultiFab::plus (Real val,
                int  comp,
                int  num_comp,
                int  nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        mfi().plus(val, b, comp, num_comp);
    }
}

void
MultiFab::plus (Real       val,
                const Box& region,
                int        comp,
                int        num_comp,
                int        nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        b &= region;
        if (b.ok())
        {
            mfi().plus(val, b, comp, num_comp);
        }
    }
}

void
MultiFab::plus (const MultiFab& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost)
{
    assert(boxarray == mf.boxarray);
    assert(strt_comp >= 0);
#ifndef NDEBUG
    int lst_comp = strt_comp + num_comp - 1;
#endif
    assert(lst_comp < n_comp && lst_comp < mf.n_comp);
    assert(nghost <= n_grow && nghost <= mf.n_grow);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator dmfi(mfi, mf);
        Box bx(mfi.validbox());
        bx.grow(nghost);
        mfi().plus(dmfi(), bx, strt_comp, strt_comp, num_comp);
    }
}

void
MultiFab::mult (Real val,
                int  comp,
                int  num_comp,
                int  nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        mfi().mult(val, b, comp, num_comp);
    }
}

void
MultiFab::mult (Real       val,
                const Box& region,
                int        comp,
                int        num_comp,
                int        nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        b &= region;
        if(b.ok())
        {
            mfi().mult(val, b, comp, num_comp);
        }
    }
}

void
MultiFab::invert (Real numerator,
                  int  comp,
                  int  num_comp,
                  int  nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        mfi().invert(numerator, b, comp, num_comp);
    }
}

void
MultiFab::invert (Real       numerator,
                  const Box& region,
                  int        comp,
                  int        num_comp,
                  int        nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        b &= region;
        if(b.ok())
        {
            mfi().invert(numerator, b, comp, num_comp);
        }
    }
}

void
MultiFab::negate (int comp,
                  int num_comp,
                  int nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        mfi().negate(b, comp, num_comp);
    }
}

void
MultiFab::negate (const Box& region,
                  int        comp,
                  int        num_comp,
                  int        nghost)
{
    assert(nghost >= 0 && nghost <= n_grow);
    assert(comp+num_comp <= n_comp);

    for (MultiFabIterator mfi(*this); mfi.isValid(); ++mfi)
    {
        Box b(grow(mfi.validbox(), nghost));
        b &= region;
        if (b.ok())
        {
            mfi().negate(b, comp, num_comp);
        }
    }
}

void
linInterp (FArrayBox&      dest,
           const Box&      subbox,
           const MultiFab& f1,
           const MultiFab& f2,
           Real            t1,
           Real            t2,
           Real            t,
           bool            extrap)
{
    const Real teps = (t2-t1)/1000.0;

    assert(t>t1-teps && (extrap || t < t2+teps));

    if (t < t1+teps)
        f1.copy(dest,subbox);
    else if (t > t2-teps && t < t2+teps)
        f2.copy(dest,subbox);
    else
    {
        const int       nv  = dest.nComp();
        const BoxArray& boxarray2 = f2.boxArray();

        assert(f1.boxArray() == boxarray2);
        assert(nv == f1.n_comp && nv == f2.n_comp);

        const int dc = 0;
        const int sc = 0;

        cerr << "Error in MultiFab::linInterp 1:  fix for parallel" << endl;
        ParallelDescriptor::Abort("Error: MultiFab::linInterp 1:  fix parallel.");

        for (ConstMultiFabIterator mfi(f2); mfi.isValid(); ++mfi)
        {
	    ConstDependentMultiFabIterator dmfi(mfi, f1);
            if (mfi.validbox().intersects(subbox))
            {
                Box destbox(mfi.validbox());
                destbox &= subbox;
                dest.linInterp(dmfi(), destbox, sc,
                               mfi(),  destbox, sc,
                               t1, t2, t, destbox, dc, nv);
            }
        }
    }
}

void
linInterp (FArrayBox&      dest,
           const Box&      subbox,
           const MultiFab& f1,
           const MultiFab& f2,
           Real            t1,
           Real            t2,
           Real            t,
           int             src_comp,
           int             dest_comp,
           int             num_comp,
           bool            extrap)
{
    const Real teps = (t2-t1)/1000.0;

    assert(t>t1-teps && (extrap || t < t2+teps));

    if (t < t1+teps)
        f1.copy(dest,subbox,src_comp,dest_comp,num_comp);
    else if (t > t2-teps && t < t2+teps)
        f2.copy(dest,subbox,src_comp,dest_comp,num_comp);
    else
    {
        const BoxArray& boxarray2 = f2.boxArray();

        assert(f1.boxArray() == boxarray2);
        assert(f1.n_comp == f2.n_comp);
        assert(src_comp + num_comp <= f1.n_comp);
        assert(dest_comp + num_comp <= dest.nComp());

        const int dc = dest_comp;
        const int sc = src_comp;

        if (ParallelDescriptor::NProcs() > 1)
        {
            ParallelDescriptor::Abort("MultiFab::linInterp 2 not implemented in parallel.");
        }
        else
        {
            cerr << "MultiFab::linInterp 2 not implemented in parallel." << endl;
        }

        for (ConstMultiFabIterator mfi(f2); mfi.isValid(); ++mfi)
        {
            ConstDependentMultiFabIterator dmfi(mfi, f1);
            if (mfi.validbox().intersects(subbox))
            {
                //
                // Restrict copy to domain of validity of source.
                //
                Box destbox(mfi.validbox());
                destbox &= subbox;

                dest.linInterp(dmfi(), destbox, sc,
                               mfi(),  destbox, sc,
                               t1, t2, t, destbox, dc, num_comp);
            }
        }
    }
}

void
linInterpAddBox (MultiFabCopyDescriptor& fabCopyDesc,
                 BoxList&                returnUnfilledBoxes,
                 Array<FillBoxId>&       returnedFillBoxIds,
                 const Box&              subbox,
                 const MultiFabId&       faid1,
                 const MultiFabId&       faid2,
                 Real                    t1,
                 Real                    t2,
                 Real                    t,
                 int                     src_comp,
                 int                     dest_comp,
                 int                     num_comp,
                 bool                    extrap)
{
    const Real teps = (t2-t1)/1000.0;

    assert(t>t1-teps && (extrap || t < t2+teps));

    if (t < t1+teps)
    {
	returnedFillBoxIds.resize(1);
	returnedFillBoxIds[0] = fabCopyDesc.AddBox(faid1, subbox,
                                                   returnUnfilledBoxes,
                                                   src_comp, dest_comp,
                                                   num_comp);
    }
    else if (t > t2-teps && t < t2+teps)
    {
	returnedFillBoxIds.resize(1);
	returnedFillBoxIds[0] = fabCopyDesc.AddBox(faid2, subbox,
                                                   returnUnfilledBoxes,
                                                   src_comp, dest_comp,
                                                   num_comp);
    }
    else
    {
        //assert(f1.boxArray() == boxarray2);
        //assert(f1.n_comp == f2.n_comp);
        //assert(src_comp + num_comp <= f1.n_comp);
        //assert(dest_comp + num_comp <= dest.nComp());

	returnedFillBoxIds.resize(2);
	BoxList tempUnfilledBoxes(subbox.ixType());
	returnedFillBoxIds[0] = fabCopyDesc.AddBox(faid1, subbox,
				       returnUnfilledBoxes,
				       src_comp, dest_comp, num_comp);
	returnedFillBoxIds[1] = fabCopyDesc.AddBox(faid2, subbox,
				       tempUnfilledBoxes,
				       src_comp, dest_comp, num_comp);
	// note:  the boxarrays for faid1 and faid2 should be the
	//        same so only use returnUnfilledBoxes from one AddBox here
    }
}

void
linInterpFillFab (MultiFabCopyDescriptor& fabCopyDesc,
                  const Array<FillBoxId>& fillBoxIds,
                  const MultiFabId&       faid1,
                  const MultiFabId&       faid2,
                  FArrayBox&              dest,
                  Real                    t1,
                  Real                    t2,
                  Real                    t,
                  int                     src_comp,   // these comps need to be removed
                  int                     dest_comp,  // from this routine
                  int                     num_comp,
                  bool                    extrap)
{
    const Real teps = (t2-t1)/1000.0;

    assert(t>t1-teps && (extrap || t < t2+teps));

    if (t < t1+teps)
    {
	fabCopyDesc.FillFab(faid1, fillBoxIds[0], dest);
    }
    else if (t > t2-teps && t < t2+teps)
    {
	fabCopyDesc.FillFab(faid2, fillBoxIds[0], dest);
    }
    else
    {
        assert(dest_comp + num_comp <= dest.nComp());

	FArrayBox dest1(dest.box(), dest.nComp());
	dest1.setVal(1.e30);
	FArrayBox dest2(dest.box(), dest.nComp());
	dest2.setVal(1.e30);
	fabCopyDesc.FillFab(faid1, fillBoxIds[0], dest1);
	fabCopyDesc.FillFab(faid2, fillBoxIds[1], dest2);
        dest.linInterp(dest1, dest1.box(), src_comp,
                       dest2, dest2.box(), src_comp,
                       t1, t2, t, dest.box(), dest_comp, num_comp);
    }
}
