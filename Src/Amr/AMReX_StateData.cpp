
#include <iostream>
#include <algorithm>

#include <unistd.h>

#include <AMReX_RealBox.H>
#include <AMReX_StateData.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

const Real INVALID_TIME = -1.0e200;

const int MFNEWDATA = 0;
const int MFOLDDATA = 1;

Vector<std::string> StateData::fabArrayHeaderNames;
std::map<std::string, Vector<char> > *StateData::faHeaderMap;


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
		      const DistributionMapping& dm,
                      const StateDescriptor* d,
                      Real                   cur_time,
                      Real                   dt,
                      const FabFactory<FArrayBox>& factory)
{
    define(p_domain, grds, dm, *d, cur_time, dt, factory);
}

void
StateData::Initialize (StateData& dest, const StateData& src)
{

    // Define the object with the same properties as src.

    dest.define(src.getDomain(), src.boxArray(), src.DistributionMap(),
		*(src.descriptor()), src.curTime(), src.curTime() - src.prevTime(),
                src.Factory());

    // Now, for both the new data and the old data if it's there,
    // generate a new MultiFab for the dest and remove the previous data.

    if (src.hasOldData()) {
      dest.allocOldData();
      dest.copyOld(src);
    }

    if (src.hasNewData())
      dest.copyNew(src);

}

void
StateData::define (const Box&             p_domain,
                   const BoxArray&        grds,
		   const DistributionMapping& dm,
                   const StateDescriptor& d,
                   Real                   time,
                   Real                   dt,
                   const FabFactory<FArrayBox>& factory)
{
    BL_PROFILE("StateData::define()");
    domain = p_domain;
    desc = &d;
    grids = grds;
    dmap = dm;
    m_factory.reset(factory.clone());
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

    MFInfo info;

#ifdef AMREX_USE_DEVICE
    info.SetDevice(d.deviceCopy());
#endif

    new_data = new MultiFab(grids,dmap,ncomp,desc->nExtra(), info, *m_factory);

    old_data = 0;
}

void
StateData::copyOld (const StateData& state)
{

  BL_ASSERT(state.hasOldData());
  BL_ASSERT(old_data != 0);

  const MultiFab& MF = state.oldData();

  int nc = MF.nComp();
  int ng = MF.nGrow();

  BL_ASSERT(nc == (*old_data).nComp());
  BL_ASSERT(ng == (*old_data).nGrow());

  MultiFab::Copy(*old_data, state.oldData(), 0, 0, nc, ng);

  StateDescriptor::TimeCenter t_typ(desc->timeType());

  if (t_typ == StateDescriptor::Point)
    {
      old_time.start = old_time.stop = state.prevTime();
    }
  else
    {
      Real dt = state.curTime() - state.prevTime();

      old_time.start = state.prevTime() - dt/2.0;
      old_time.stop  = state.prevTime() + dt/2.0;
    }

}

void
StateData::copyNew (const StateData& state)
{

  BL_ASSERT(state.hasNewData());
  BL_ASSERT(new_data != 0);

  const MultiFab& MF = state.newData();

  int nc = MF.nComp();
  int ng = MF.nGrow();

  BL_ASSERT(nc == (*new_data).nComp());
  BL_ASSERT(ng == (*new_data).nGrow());

  MultiFab::Copy(*new_data, state.newData(), 0, 0, nc, ng);

  StateDescriptor::TimeCenter t_typ(desc->timeType());

  if (t_typ == StateDescriptor::Point)
    {
      new_time.start = new_time.stop = state.curTime();
    }
  else
    {
      Real dt = state.curTime() - state.prevTime();

      new_time.start = state.curTime() - dt/2.0;
      new_time.stop  = state.curTime() + dt/2.0;
    }

}

void
StateData::reset ()
{
    new_time = old_time;
    old_time.start = old_time.stop = INVALID_TIME;
    std::swap(old_data, new_data);
}

void
StateData::restart (std::istream&          is,
		    const Box&             p_domain,
		    const BoxArray&        grds,
		    const DistributionMapping& dm,
                    const FabFactory<FArrayBox>& factory,
                    const StateDescriptor& d,
                    const std::string&     chkfile)
{
    desc = &d;
    domain = p_domain;
    grids = grds;
    dmap = dm;
    m_factory.reset(factory.clone());

    // Convert to proper type.
    IndexType typ(desc->getType());
    if (!typ.cellCentered()) {
        domain.convert(typ);
        grids.convert(typ);
    }

    {
	Box domain_in;
	BoxArray grids_in;
	is >> domain_in;
	grids_in.readFrom(is);
	BL_ASSERT(domain_in == domain);
	BL_ASSERT(amrex::match(grids_in,grids));
    }

    restartDoit(is, chkfile);
}

void 
StateData::restartDoit (std::istream& is, const std::string& chkfile)
{
    BL_PROFILE("StateData::restartDoit()");

    is >> old_time.start;
    is >> old_time.stop;
    is >> new_time.start;
    is >> new_time.stop;

    int nsets;
    is >> nsets;

    MFInfo info;
#ifdef AMREX_USE_DEVICE
    info.SetDevice(desc->deviceCopy());
#endif

    old_data = (nsets == 2) ? new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(),
                                           info, *m_factory)
                              : nullptr;
    new_data =                new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(),
                                           info, *m_factory);

    //
    // If no data is written then we just allocate the MF instead of reading it in. 
    // This assumes that the application will do something with it.
    // We set it to zero in case a compiler complains about uninitialized data.
    //
    if (nsets == 0) {
       new_data->setVal(0.0);
    }

    std::string mf_name;
    std::string FullPathName;

    for(int ns(1); ns <= nsets; ++ns) {
      MultiFab *whichMF = nullptr;
      if(ns == 1) {
	whichMF = new_data;
      } else if(ns == 2) {
	whichMF = old_data;
      } else {
        amrex::Abort("**** Error in StateData::restart:  invalid nsets.");
      }

      is >> mf_name;
      //
      // Note that mf_name is relative to the Header file.
      // We need to prepend the name of the chkfile directory.
      //
      FullPathName = chkfile;
      if ( ! chkfile.empty() && chkfile[chkfile.length()-1] != '/') {
          FullPathName += '/';
      }
      FullPathName += mf_name;

      // ---- check for preread header
      std::string FullHeaderPathName(FullPathName + "_H");
      const char *faHeader = 0;
      if(faHeaderMap != 0) {
        std::map<std::string, Vector<char> >::iterator fahmIter;
	fahmIter = faHeaderMap->find(FullHeaderPathName);
	if(fahmIter != faHeaderMap->end()) {
	  faHeader = fahmIter->second.dataPtr();
	}
      }

      VisMF::Read(*whichMF, FullPathName, faHeader);
    }
}

void 
StateData::restart (const StateDescriptor& d,
		    const StateData& rhs)
{
    desc = &d;
    domain = rhs.domain;
    grids = rhs.grids;
    old_time.start = rhs.old_time.start;
    old_time.stop  = rhs.old_time.stop;
    new_time.start = rhs.new_time.start;
    new_time.stop  = rhs.new_time.stop;
    old_data = 0;

    MFInfo info;
#ifdef AMREX_USE_DEVICE
    info.SetDevice(d.deviceCopy());
#endif
    new_data = new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(), info, *m_factory);
    new_data->setVal(0.);
}

StateData::~StateData()
{
   desc = 0;
   delete new_data;
   delete old_data;
}

void
StateData::allocOldData ()
{
    if (old_data == 0)
    {
	MFInfo info;
#ifdef AMREX_USE_DEVICE
	info.SetDevice(desc->deviceCopy());
#endif
        old_data = new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(), info, *m_factory);
    }
}

BCRec
StateData::getBC (int comp, int i) const
{
    BCRec bcr;
    amrex::setBC(grids[i],domain,desc->getBC(comp),bcr);
    return bcr;
}

void
StateData::setOldTimeLevel (Real time)
{
    if (desc->timeType() == StateDescriptor::Point)
    {
        old_time.start = old_time.stop = time;
    }
    else
    {
        amrex::Error("StateData::setOldTimeLevel called with Interval");
    }
}

void
StateData::setNewTimeLevel (Real time)
{
    if (desc->timeType() == StateDescriptor::Point)
    {
        new_time.start = new_time.stop = time;
    }
    else
    {
        amrex::Error("StateData::setNewTimeLevel called with Interval");
    }
}

void
StateData::syncNewTimeLevel (Real time)
{
    Real teps = (new_time.stop - old_time.stop)*1.e-3;
    if (time > new_time.stop-teps && time < new_time.stop+teps)
    {
	if (desc->timeType() == StateDescriptor::Point)
	{
	    new_time.start = new_time.stop = time;
	}
	else
	{
	    new_time.stop = time;
	}
    }
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
    std::swap(old_data, new_data);
}

void
StateData::replaceOldData (MultiFab* mf)
{
    std::swap(old_data, mf);
    delete mf;
}

void
StateData::replaceNewData (MultiFab* mf)
{
    std::swap(new_data, mf);
    delete mf;
}

void
StateData::PrepareForFillBoundary(FArrayBox& dest,
                                  const Real* dx,
                                  const RealBox& prob_domain,
                                  Real* xlo,
                                  int* bc,
                                  int dest_comp,
                                  int src_comp,
                                  int num_comp)
{
    const int* dlo = dest.loVect();
    const int* dhi = dest.hiVect();
    const int* plo = domain.loVect();

    const Real* problo = prob_domain.lo();

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        xlo[i] = problo[i] + dx[i]*(dlo[i]-plo[i]);
    }

    const Box& bx = dest.box();
    BCRec bcr;

    for (int i = 0; i < num_comp; ++i) {
        const int sc = src_comp+i;
        amrex::setBC(bx,domain,desc->getBC(sc),bcr);
        int bcidx = 2 * AMREX_SPACEDIM * i;
        for (int j = 0; j < 2 * AMREX_SPACEDIM; ++j) {
            bc[bcidx + j] = bcr.vect()[j];
        }
    }

}

void
StateData::FillBoundary (FArrayBox&     dest,
                         const Real*    time,
                         const Real*    dx,
                         const Real*    xlo,
                         const int*     bcrs,
                         const RealBox& prob_domain,
                         int            dest_comp,
                         int            src_comp,
                         int            num_comp)
{
    BL_PROFILE("StateData::FillBoundary()");
    BL_ASSERT(dest.box().ixType() == desc->getType());
   
    if (domain.contains(dest.box())) return;

    const Box& bx  = dest.box();
    const int* dlo = dest.loVect();
    const int* dhi = dest.hiVect();
    const int* plo = domain.loVect();
    const int* phi = domain.hiVect();

    for (int i = 0; i < num_comp; )
    {
        const int dc  = dest_comp+i;
        const int sc  = src_comp+i;
        Real*     dat = dest.dataPtr(dc);

        int bcidx = 2 * AMREX_SPACEDIM * i;

        if (desc->master(sc))
        {
            const int groupsize = desc->groupsize(sc);

            BL_ASSERT(groupsize != 0);

#ifdef AMREX_USE_DEVICE
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
            // Ensure that our threadblock size is such that it is
            // evenly divisible by the number of zones in the box,
            // and is at least one larger than the number of ghost zones.
            // This ensures that the corners plus one interior zone
            // are all on the same threadblock.

            const IntVect left = domain.smallEnd() - bx.smallEnd();
            const IntVect rght = bx.bigEnd() - domain.bigEnd();

            int ng[3] = {0, 0, 0};

            for (int n = 0; n < BL_SPACEDIM; ++n)
                ng[n] = std::max(0, std::max(left[n], rght[n]));

            const IntVect size = bx.size();
            IntVect numThreadsMin(ng[0] + 1, ng[1] + 1, ng[2] + 1);

            for (int n = 0; n < BL_SPACEDIM; ++n) {
                while (size[n] % numThreadsMin[n] != 0) {
                    ++numThreadsMin[n];
                }
            }

            if (std::min({numThreadsMin[0], numThreadsMin[1], numThreadsMin[2]}) < 1)
                amrex::Error("Minimum number of CUDA threads must be positive.");

            Device::setNumThreadsMin(numThreadsMin[0], numThreadsMin[1], numThreadsMin[2]);
#endif
            Device::prepare_for_launch(dest.loVect(), dest.hiVect());
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
            Device::setNumThreadsMin(1, 1, 1);
#endif
#endif

            if (groupsize+i <= num_comp)
            {
                //
                // Can do the whole group at once.
                // Use the "group" boundary fill routine.
                //
		desc->bndryFill(sc)(dat,dlo,dhi,plo,phi,dx,xlo,time,&bcrs[bcidx],groupsize);

                i += groupsize;
            }
            else
            {
                desc->bndryFill(sc)(dat,dlo,dhi,plo,phi,dx,xlo,time,&bcrs[bcidx]);

                i++;
            }
        }
        else
        {
            desc->bndryFill(sc)(dat,dlo,dhi,plo,phi,dx,xlo,time,&bcrs[bcidx]);

            i++;
        }

    }
}

void
StateData::RegisterData (MultiFabCopyDescriptor& multiFabCopyDesc,
                         Vector<MultiFabId>&      mfid)
{
    mfid.resize(2);
    mfid[MFNEWDATA] = multiFabCopyDesc.RegisterFabArray(new_data);
    mfid[MFOLDDATA] = multiFabCopyDesc.RegisterFabArray(old_data);
}

void
StateData::InterpAddBox (MultiFabCopyDescriptor& multiFabCopyDesc,
			 Vector<MultiFabId>&      mfid,
			 BoxList*                unfillableBoxes,
			 Vector<FillBoxId>&       returnedFillBoxIds,
			 const Box&              subbox,
			 Real                    time,
			 int                     src_comp,
			 int                     dest_comp,
			 int                     num_comp,
			 bool                    extrap)
{
    if (desc->timeType() == StateDescriptor::Point)
    {
        if (old_data == 0)
        {
            returnedFillBoxIds.resize(1);
            returnedFillBoxIds[0] = multiFabCopyDesc.AddBox(mfid[MFNEWDATA],
                                                            subbox,
                                                            unfillableBoxes,
                                                            src_comp,
                                                            dest_comp,
                                                            num_comp);
        }
        else
        {
            amrex::InterpAddBox(multiFabCopyDesc,
				 unfillableBoxes,
				 returnedFillBoxIds,
				 subbox,
				 mfid[MFOLDDATA],
				 mfid[MFNEWDATA],
				 old_time.start,
				 new_time.start,
				 time,
				 src_comp,
				 dest_comp,
				 num_comp,
				 extrap);
        }
    }
    else
    {
        const Real teps = (new_time.start - old_time.start)*1.e-3;

        if (time > new_time.start-teps && time < new_time.stop+teps)
        {
            returnedFillBoxIds.resize(1);
            returnedFillBoxIds[0] = multiFabCopyDesc.AddBox(mfid[MFNEWDATA],
                                                            subbox,
                                                            unfillableBoxes,
                                                            src_comp,
                                                            dest_comp,
                                                            num_comp);
        }
        else if (old_data != 0              &&
                 time > old_time.start-teps &&
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
            amrex::Error("StateData::Interp(): cannot interp");
        }
   }
}

void
StateData::InterpFillFab (MultiFabCopyDescriptor&  multiFabCopyDesc,
			  const Vector<MultiFabId>& mfid,
			  const Vector<FillBoxId>&  fillBoxIds,
			  FArrayBox&               dest,
			  Real                     time,
			  int                      src_comp,
			  int                      dest_comp,
			  int                      num_comp,
			  bool                     extrap)
{
    BL_PROFILE("StateData::InterpFillFab()");
    if (desc->timeType() == StateDescriptor::Point)
    {
        if (old_data == 0)
        {
            multiFabCopyDesc.FillFab(mfid[MFNEWDATA], fillBoxIds[0], dest);
        }
        else
        {
            amrex::InterpFillFab(multiFabCopyDesc,
				  fillBoxIds,
				  mfid[MFOLDDATA],
				  mfid[MFNEWDATA],
				  dest,
				  old_time.start,
				  new_time.start,
				  time,
				  src_comp,
				  dest_comp,
				  num_comp,
				  extrap);
        }
    }
    else
    {
        const Real teps = (new_time.start - old_time.start)*1.e-3;

        if (time > new_time.start-teps && time < new_time.stop+teps)
        {
            multiFabCopyDesc.FillFab(mfid[MFNEWDATA], fillBoxIds[0], dest);
        }
        else if (old_data != 0              &&
                 time > old_time.start-teps &&
                 time < old_time.stop+teps)
        {
            multiFabCopyDesc.FillFab(mfid[MFOLDDATA], fillBoxIds[0], dest);
        }
        else
        {
            amrex::Error("StateData::Interp(): cannot interp");
        }
    }
}

void
StateData::getData (Vector<MultiFab*>& data,
		    Vector<Real>& datatime,
		    Real time) const
{
    data.clear();
    datatime.clear();

    if (desc->timeType() == StateDescriptor::Point)
    {
	BL_ASSERT(new_data != 0);
        if (old_data == 0)
        {
	    data.push_back(new_data);
	    datatime.push_back(new_time.start);
        }
        else
        {
	    const Real teps = (new_time.start - old_time.start)*1.e-3;
	    if (time > new_time.start-teps && time < new_time.start+teps) {
		data.push_back(new_data);
		datatime.push_back(new_time.start);
	    } else if (time > old_time.start-teps && time < old_time.start+teps) {
	    	    data.push_back(old_data);
		    datatime.push_back(old_time.start);
	    } else {
		data.push_back(old_data);
		data.push_back(new_data);
		datatime.push_back(old_time.start);
		datatime.push_back(new_time.start);
	    }
        }
    }
    else
    {
        const Real teps = (new_time.start - old_time.start)*1.e-3;

        if (time > new_time.start-teps && time < new_time.stop+teps)
        {
	    data.push_back(new_data);
	    datatime.push_back(time);
        }
        else if (old_data != 0              &&
                 time > old_time.start-teps &&
                 time < old_time.stop+teps)
        {
	    data.push_back(old_data);
	    datatime.push_back(time);
        }
        else
        {
            amrex::Error("StateData::getData(): how did we get here?");
        }
    }
}

void
StateData::checkPoint (const std::string& name,
                       const std::string& fullpathname,
                       std::ostream&  os,
                       VisMF::How     how,
                       bool           dump_old)
{
    BL_PROFILE("StateData::checkPoint()");
    static const std::string NewSuffix("_New_MF");
    static const std::string OldSuffix("_Old_MF");

    if (dump_old == true && old_data == 0)
    {
        dump_old = false;
    }

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // The relative name gets written to the Header file.
        //
        std::string mf_name_old(name + OldSuffix);
        std::string mf_name_new(name + NewSuffix);

        os << domain << '\n';

        grids.writeOn(os);

        os << old_time.start << '\n'
           << old_time.stop  << '\n'
           << new_time.start << '\n'
           << new_time.stop  << '\n';

        if (desc->store_in_checkpoint()) 
        {
           if (dump_old)
           {
               os << 2 << '\n' << mf_name_new << '\n' << mf_name_old << '\n';
	       fabArrayHeaderNames.push_back(mf_name_new);
	       fabArrayHeaderNames.push_back(mf_name_old);
           }
           else
           {
               os << 1 << '\n' << mf_name_new << '\n';
	       fabArrayHeaderNames.push_back(mf_name_new);
           }
        }
        else
        {
               os << 0 << '\n';
        }
    }

    if (desc->store_in_checkpoint())
    {
       BL_ASSERT(new_data);
       std::string mf_fullpath_new(fullpathname + NewSuffix);
       VisMF::Write(*new_data,mf_fullpath_new,how);

       if (dump_old)
       {
           BL_ASSERT(old_data);
           std::string mf_fullpath_old(fullpathname + OldSuffix);
           VisMF::Write(*old_data,mf_fullpath_old,how);
       }
    }
}

void
StateData::printTimeInterval (std::ostream &os) const
{
    os << '['
       << old_time.start
       << ' '
       << old_time.stop
       << "] ["
       << new_time.start
       << ' '
       << new_time.stop
       << ']'
       << '\n';
}

StateDataPhysBCFunct::StateDataPhysBCFunct (StateData&sd, int sc, const Geometry& geom_)
    : statedata(&sd),
      src_comp(sc),
      geom(geom_)
{ }

void
StateDataPhysBCFunct::FillBoundary (MultiFab& mf, int dest_comp, int num_comp, Real time)
{
    BL_PROFILE("StateDataPhysBCFunct::FillBoundary");

    const Box&     domain      = statedata->getDomain();
    const int*     domainlo    = domain.loVect();
    const int*     domainhi    = domain.hiVect();
    const Real*    dx          = geom.CellSize();
    const RealBox& prob_domain = geom.ProbDomain();

#ifdef CRSEGRNDOMP
#ifdef _OPENMP
#pragma omp parallel
#endif
#endif
    {
	FArrayBox tmp;

 	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
 	{
	    FArrayBox& dest = mf[mfi];
	    const Box& bx = dest.box();
	    
	    bool has_phys_bc = false;
	    bool is_periodic = false;
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		bool touch = bx.smallEnd(i) < domainlo[i] || bx.bigEnd(i) > domainhi[i];
		if (geom.isPeriodic(i)) {
		    is_periodic = is_periodic || touch;
		} else {
		    has_phys_bc = has_phys_bc || touch;
		}
	    }
	    
	    if (has_phys_bc)
	    {

                Real xlo[BL_SPACEDIM];
                int bcrs[2 * AMREX_SPACEDIM * num_comp];

                statedata->PrepareForFillBoundary(dest, dx, prob_domain, xlo,
                                                  bcrs, dest_comp, src_comp, num_comp);

#ifdef AMREX_USE_CUDA
                Real* time_f = mfi.get_fortran_pointer(&time);
                Real* xlo_f  = mfi.get_fortran_pointer(xlo, 3, AMREX_SPACEDIM);
                int* bcrs_f  = mfi.get_fortran_pointer(bcrs, 2 * AMREX_SPACEDIM * num_comp);
#else
                Real* time_f = &time;
                Real* xlo_f  = xlo;
                int* bcrs_f  = bcrs;
#endif

		statedata->FillBoundary(dest, time_f, dx, xlo_f, bcrs_f, prob_domain, dest_comp, src_comp, num_comp);

		if (is_periodic) // fix up corner
		{
		    Box GrownDomain = domain;
		    
		    for (int dir = 0; dir < BL_SPACEDIM; dir++)
		    {
			if (!geom.isPeriodic(dir))
			{
			    const int lo = domainlo[dir] - bx.smallEnd(dir);
			    const int hi = bx.bigEnd(dir) - domainhi[dir];
			    if (lo > 0) GrownDomain.growLo(dir,lo);
			    if (hi > 0) GrownDomain.growHi(dir,hi);
			}
		    }
		    
		    for (int dir = 0; dir < BL_SPACEDIM; dir++)
		    {
			if (!geom.isPeriodic(dir)) continue;
			
			Box lo_slab = bx;
			Box hi_slab = bx;
			lo_slab.shift(dir, domain.length(dir));
			hi_slab.shift(dir,-domain.length(dir));
			lo_slab &= GrownDomain;
			hi_slab &= GrownDomain;
			
			if (lo_slab.ok())
			{
			    lo_slab.shift(dir,-domain.length(dir));
			    
			    tmp.resize(lo_slab,num_comp);
			    tmp.copy(dest,dest_comp,0,num_comp);
			    tmp.shift(dir,domain.length(dir));
			    
			    statedata->FillBoundary(tmp, time_f, dx, xlo_f, bcrs_f, prob_domain, dest_comp, src_comp, num_comp);
			    
			    tmp.shift(dir,-domain.length(dir));
			    dest.copy(tmp,0,dest_comp,num_comp);
			}
			
			if (hi_slab.ok())
			{
			    hi_slab.shift(dir,domain.length(dir));
			    
			    tmp.resize(hi_slab,num_comp);
			    tmp.copy(dest,dest_comp,0,num_comp);
			    tmp.shift(dir,-domain.length(dir));
			    
			    statedata->FillBoundary(tmp, time_f, dx, xlo_f, bcrs_f, prob_domain, dest_comp, src_comp, num_comp);
			    
			    tmp.shift(dir,domain.length(dir));
			    dest.copy(tmp,0,dest_comp,num_comp);
			}
		    }
		}
 	    }
 	}
    }
}


void StateData::AddProcsToComp(const StateDescriptor &sdPtr,
                               int ioProcNumSCS, int ioProcNumAll,
                               int scsMyId, MPI_Comm scsComm)
{
#if BL_USE_MPI
      // ---- StateDescriptor
      desc = &sdPtr;

      // ---- TimeIntervals
      ParallelDescriptor::Bcast(&new_time.start, 1, ioProcNumSCS, scsComm);
      ParallelDescriptor::Bcast(&new_time.stop,  1, ioProcNumSCS, scsComm);
      ParallelDescriptor::Bcast(&old_time.start, 1, ioProcNumSCS, scsComm);
      ParallelDescriptor::Bcast(&old_time.stop,  1, ioProcNumSCS, scsComm);

      // ---- Boxes
      amrex::BroadcastBox(domain, scsMyId, ioProcNumSCS, scsComm);

      // ---- BoxArrays
      amrex::BroadcastBoxArray(grids, scsMyId, ioProcNumSCS, scsComm);

      // ---- MultiFabs
      int makeNewDataId(-7), makeOldDataId(-7);
      if(new_data != 0) {
	makeNewDataId = new_data->AllocatedFAPtrID();
      }
      ParallelDescriptor::Bcast(&makeNewDataId, 1, ioProcNumSCS, scsComm);
      if(scsMyId != ioProcNumSCS) {
        if(makeNewDataId >= 0) {
          new_data = new MultiFab;
        } else {
          new_data = 0;
	}
      }
      if(new_data != 0) {
        new_data->AddProcsToComp(ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
      }

      if(old_data != 0) {
	makeOldDataId = old_data->AllocatedFAPtrID();
      }
      ParallelDescriptor::Bcast(&makeOldDataId, 1, ioProcNumSCS, scsComm);
      if(scsMyId != ioProcNumSCS) {
        if(makeOldDataId >= 0) {
          old_data = new MultiFab;
        } else {
          old_data = 0;
	}
      }
      if(old_data != 0) {
        old_data->AddProcsToComp(ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
      }

      ParallelDescriptor::Barrier(scsComm);
#endif
}


void StateData::Check() const
{
      if(new_data != 0) {
        new_data->DistributionMap().Check();
      }
      if(old_data != 0) {
        old_data->DistributionMap().Check();
      }
}

}
