//BL_COPYRIGHT_NOTICE

//
// $Id: StateData.cpp,v 1.14 1998-03-30 20:04:18 lijewski Exp $
//

#include <RunStats.H>
#include <StateData.H>
#include <ParallelDescriptor.H>

#if defined(BL_USE_FLOAT) && !defined(BL_ARCH_CRAY)
const Real INVALID_TIME = -1.0e30;
#else
const Real INVALID_TIME = -1.0e100;
#endif

const int MFNEWDATA = 0;
const int MFOLDDATA = 1;

#ifdef __GNUG__
StateData::StateData(const StateData&)
{
    BoxLib::Error("Don't call this function!");
}
StateData& StateData::operator= (const StateData&)
{
    BoxLib::Error("Don't call this function!");
}
#endif

StateData::StateData () 
{
   desc = 0;
   new_data = old_data = 0;
   new_time.start = INVALID_TIME;
   new_time.stop  = INVALID_TIME;
   old_time.start = INVALID_TIME;
   old_time.stop  = INVALID_TIME;
}

StateData::StateData (const Box&             p_domain,
                      const BoxArray&        grds,
                      const StateDescriptor* d,
                      Real                   time,
                      Real                   dt)
{
    define(p_domain, grds, *d, time, dt);
}

void
StateData::define (const Box&             p_domain,
                   const BoxArray&        grds,
                   const StateDescriptor& d,
                   Real                   time,
                   Real                   dt)
{
    domain = p_domain;
    desc = &d;
    grids.define(grds);
    //
    // Convert to proper type.
    //
    IndexType typ(desc->getType());
    StateDescriptor::TimeCenter t_typ(desc->timeType());
    if (!typ.cellCentered())
    {
        domain.convert(typ);
        grids.convert(typ);
    }
    if (t_typ == StateDescriptor::Point)
    {
        new_time.start = new_time.stop = time;
        old_time.start = old_time.stop = time - dt;
    }
    else
    {
        new_time.start = time;
        new_time.stop  = time+dt;
        old_time.start = time-dt;
        old_time.stop  = time;
    }
    int ncomp = desc->nComp();

    new_data = new MultiFab(grids,ncomp,desc->nExtra(),Fab_allocate);

    old_data = 0;
    buildBC();
}

void
StateData::restart (istream&               is,
                    const StateDescriptor& d,
                    const aString&         chkfile)
{
    desc = &d;

    is >> domain;
    grids.define(is);

    is >> old_time.start;
    is >> old_time.stop;
    is >> new_time.start;
    is >> new_time.stop;

    int nsets;
    is >> nsets;

    old_data = 0;

    new_data = new MultiFab;
    aString mf_name;
    is >> mf_name;
    //
    // Note that mf_name is relative to the Header file.
    // We need to prepend the name of the chkfile directory.
    //
    aString FullPathName = chkfile;
    if (!chkfile.isNull() && chkfile[chkfile.length()-1] != '/')
        FullPathName += '/';
    FullPathName += mf_name;
    VisMF::Read(*new_data, FullPathName);

    if (nsets == 2)
    {
        old_data = new MultiFab;
        is >> mf_name;
        //
        // Note that mf_name is relative to the Header file.
        // We need to prepend the name of the chkfile directory.
        //
        FullPathName = chkfile;
        if (!chkfile.isNull() && chkfile[chkfile.length()-1] != '/')
            FullPathName += '/';
        FullPathName += mf_name;
        VisMF::Read(*old_data, FullPathName);
    }

    buildBC();
}

void
StateData::buildBC ()
{
    int ncomp = desc->nComp();
    bc.resize(ncomp);
    int i;
    for (i = 0; i < ncomp; i++)
    {
        bc[i].resize(grids.length());
        int j;
        for (j = 0; j < grids.length(); j++)
        {
            BCRec bcr;
            setBC(grids[j],domain,desc->getBC(i),bcr);
            bc[i].set(j,bcr);
        }
    }
}

StateData::~StateData()
{
   desc = 0;
   delete new_data;
   delete old_data;
}

void
StateData::setTimeLevel (Real time,
                         Real dt_old,
                         Real dt_new)
{
    if (desc->timeType() == StateDescriptor::Point)
    {
        new_time.start = new_time.stop = time;
        old_time.start = old_time.stop = time - dt_old;
    }
    else
    {
        new_time.start = time;
        new_time.stop  = time+dt_new;
        old_time.start = time-dt_old;
        old_time.stop  = time;
    }
}

void
StateData::swapTimeLevels (Real dt)
{
    old_time = new_time;
    if (desc->timeType() == StateDescriptor::Point)
    {
        new_time.start += dt;
        new_time.stop  += dt;
    }
    else
    {
        new_time.start = new_time.stop;
        new_time.stop += dt;
    }
    Swap(old_data, new_data);
}

void
StateData::allocOldData ()
{
    if (old_data == 0)
    {
          old_data = new MultiFab(grids,desc->nComp(),desc->nExtra(),Fab_allocate);
    }
}

void
StateData::removeOldData()
{
    delete old_data;
    old_data = 0;
}

void
StateData::reset ()
{
    new_time = old_time;
    old_time.start = old_time.stop = INVALID_TIME;
    Swap(old_data, new_data);
}

void
StateData::FillBoundary (const Real*    dx,
                         const RealBox& prob_domain,
                         int            src_comp,
                         int            num_comp,
                         int            do_new)
{
    Real cur_time;
    if (desc->timeType() == StateDescriptor::Point)
    {
        cur_time = new_time.start;
        if (!do_new)
            cur_time = old_time.start;
    }
    else
    {
        cur_time = 0.5*(new_time.start + new_time.stop);
        if (!do_new)
            cur_time = 0.5*(old_time.start + old_time.stop);
    }

    MultiFabIterator new_datamfi(*new_data);
   
    for ( ; new_datamfi.isValid(); ++new_datamfi)
    {
        DependentMultiFabIterator old_datamfi(new_datamfi, *old_data);
        FArrayBox* dest = &(new_datamfi());
        if (!do_new)
            dest = &(old_datamfi());
        const Box& bx = dest->box();
        if (!domain.contains(bx))
        {
            const int* dlo = bx.loVect();
            const int* dhi = bx.hiVect();
            const int* plo = domain.loVect();
            const int* phi = domain.hiVect();
            Real xlo[BL_SPACEDIM];
            BCRec bcr;
            const Real* problo = prob_domain.lo();
            int i;
            for (i = 0; i < BL_SPACEDIM; i++)
            {
                xlo[i] = problo[i] + dx[i]*(dlo[i]-plo[i]);
            }
            for (i = 0; i < num_comp; i++)
            {
                int sc = src_comp+i;
                Real* dat = dest->dataPtr(sc);
                setBC(bx,domain,desc->getBC(sc),bcr);
                desc->bndryFill(sc)(dat,ARLIM(dlo),ARLIM(dhi),
                                    plo,phi,dx,xlo, &cur_time,bcr.vect());
            }                                  
        }
    }
}

void
StateData::FillBoundary (FArrayBox&     dest,
                         Real           time,
                         const Real*    dx,
                         const RealBox& prob_domain,
                         int            dest_comp,
                         int            src_comp,
                         int            num_comp)
{
    Box dbox(dest.box());
    assert( dbox.ixType() == desc->getType() );
   
    if (domain.contains(dbox)) return;

    dbox &= domain;
    const Box& bx = dest.box();
    const int* dlo = dest.loVect();
    const int* dhi = dest.hiVect();
    const int* plo = domain.loVect();
    const int* phi = domain.hiVect();
    Real xlo[BL_SPACEDIM];
    BCRec bcr;
    const Real* problo = prob_domain.lo();
    int i;
    for (i = 0; i < BL_SPACEDIM; i++)
    {
        xlo[i] = problo[i] + dx[i]*(dlo[i]-plo[i]);
    }
    for (i = 0; i < num_comp; i++)
    {
        int dc = dest_comp+i;
        int sc = src_comp+i;
        Real* dat = dest.dataPtr(dc);
        setBC(bx,domain,desc->getBC(sc),bcr);
        desc->bndryFill(sc)(dat,ARLIM(dlo),ARLIM(dhi),
                            plo,phi,dx,xlo,&time,bcr.vect());
    }                                  
}

void
StateData::RegisterData (MultiFabCopyDescriptor& multiFabCopyDesc,
                         Array<MultiFabId>&      mfid)
{
    mfid.resize(2);
    mfid[MFNEWDATA] = multiFabCopyDesc.RegisterFabArray(new_data);
    mfid[MFOLDDATA] = multiFabCopyDesc.RegisterFabArray(old_data);
}

void
StateData::linInterpAddBox (MultiFabCopyDescriptor& multiFabCopyDesc,
                            Array<MultiFabId>&      mfid,
                            BoxList&                unfillableBoxes,
                            Array<FillBoxId>&       returnedFillBoxIds,
                            const Box&              subbox,
                            Real                    time,
                            int                     src_comp,
                            int                     dest_comp,
                            int                     num_comp,
                            bool                    extrap)
{
    Real teps = (new_time.start - old_time.start)/1000.0;
    if (desc->timeType() == StateDescriptor::Point)
    {
        if (old_data == 0)
        {
            teps = new_time.start/10000.0;
            Real dt = time - new_time.start;
            if (dt < 0.0) dt = -dt;
            returnedFillBoxIds.resize(1);
            returnedFillBoxIds[0] = multiFabCopyDesc.AddBox(mfid[MFNEWDATA],
                                                            subbox,
                                                            unfillableBoxes,
                                                            src_comp,
                                                            dest_comp,
                                                            num_comp);
            return;
        }

      ::linInterpAddBox(multiFabCopyDesc, unfillableBoxes, returnedFillBoxIds,
                  subbox, mfid[MFOLDDATA], mfid[MFNEWDATA],
                  old_time.start, new_time.start, time,
                  src_comp, dest_comp, num_comp,
                  extrap);
   }
    else
    {
      if (time > new_time.start-teps && time < new_time.stop+teps)
      {
         returnedFillBoxIds.resize(1);
         returnedFillBoxIds[0] = multiFabCopyDesc.AddBox(mfid[MFNEWDATA],
                                                         subbox,
                                                         unfillableBoxes,
                                                         src_comp,
                                                         dest_comp,
                                                         num_comp);
      } else if (old_data != 0 && time > old_time.start-teps &&
                 time < old_time.stop+teps)
      {
         returnedFillBoxIds.resize(1);
         returnedFillBoxIds[0] = multiFabCopyDesc.AddBox(mfid[MFOLDDATA],
                                                         subbox,
                                                         unfillableBoxes,
                                                         src_comp,
                                                         dest_comp,
                                                         num_comp);
      }
      else
      {
         BoxLib::Error("StateData::linInterp(): cannot interp");
      }
   }
}

void
StateData::linInterpFillFab (MultiFabCopyDescriptor&  multiFabCopyDesc,
                             const Array<MultiFabId>& mfid,
                             const Array<FillBoxId>&  fillBoxIds,
                             FArrayBox&               dest,
                             Real                     time,
                             int                      src_comp,
                             int                      dest_comp,
                             int                      num_comp,
                             bool                     extrap)
{
    Real teps = (new_time.start - old_time.start)/1000.0;
    if (desc->timeType() == StateDescriptor::Point)
    {
        if (old_data == 0)
        {
            teps = new_time.start/10000.0;
            Real dt = time - new_time.start;
            if (dt < 0.0) dt = -dt;
            multiFabCopyDesc.FillFab(mfid[MFNEWDATA], fillBoxIds[0], dest);
            return;
        }

      ::linInterpFillFab(multiFabCopyDesc, fillBoxIds,
                         mfid[MFOLDDATA], mfid[MFNEWDATA],
                         dest,
                         old_time.start, new_time.start, time,
                         src_comp, dest_comp, num_comp,
                         extrap);
   }
    else
    {
      if (time > new_time.start-teps && time < new_time.stop+teps)
      {
         multiFabCopyDesc.FillFab(mfid[MFNEWDATA], fillBoxIds[0], dest);
      }
      else if (old_data != 0 && time > old_time.start-teps &&
                 time < old_time.stop+teps) {
         multiFabCopyDesc.FillFab(mfid[MFOLDDATA], fillBoxIds[0], dest);
      }
      else
      {
         BoxLib::Error("StateData::linInterp(): cannot interp");
      }
   }
}

void
StateData::checkPoint (const aString& name,
                       const aString& fullpathname,
                       ostream&       os,
                       VisMF::How     how,
                       bool           dump_old)
{
    static const aString NewSuffix("_New_MF");
    static const aString OldSuffix("_Old_MF");

    if (dump_old == true && old_data == 0)
    {
        dump_old = false;
    }

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // The relative name gets written to the Header file.
        //
        aString mf_name_old = name; mf_name_old += OldSuffix;
        aString mf_name_new = name; mf_name_new += NewSuffix;

        os << domain << '\n';

        grids.writeOn(os);

        os << old_time.start << '\n'
           << old_time.stop  << '\n'
           << new_time.start << '\n'
           << new_time.stop  << '\n';

        if (dump_old)
        {
            os << 2 << '\n' << mf_name_new << '\n' << mf_name_old << '\n';
        }
        else
        {
            os << 1 << '\n' << mf_name_new << '\n';
        }
    }

    assert(new_data);
    aString mf_fullpath_new = fullpathname; mf_fullpath_new += NewSuffix;
    RunStats::addBytes(VisMF::Write(*new_data, mf_fullpath_new, how));

    if (dump_old)
    {
        assert(old_data);
        aString mf_fullpath_old = fullpathname; mf_fullpath_old += OldSuffix;
        RunStats::addBytes(VisMF::Write(*old_data, mf_fullpath_old, how));
    }
}

void
StateData::printTimeInterval (ostream &os)
{
    os << '['
       << old_time.start
       << ' '
       << old_time.stop
       << "] ["
       << new_time.start
       << ' ' << new_time.stop
       << ']'
       << '\n';
}
