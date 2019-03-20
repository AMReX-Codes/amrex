
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

static constexpr Real INVALID_TIME = -1.0e200;

static constexpr int MFNEWDATA = 0;
static constexpr int MFOLDDATA = 1;

Vector<std::string> StateData::fabArrayHeaderNames;
std::map<std::string, Vector<char> > *StateData::faHeaderMap;


StateData::StateData () 
    : desc(nullptr),
      new_time{INVALID_TIME,INVALID_TIME},
      old_time{INVALID_TIME,INVALID_TIME}
{
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

StateData::StateData (StateData&& rhs) noexcept
    : m_factory(std::move(rhs.m_factory)),
      desc(rhs.desc),
      domain(rhs.domain),
      grids(std::move(rhs.grids)),
      dmap(std::move(rhs.dmap)),
      new_time(rhs.new_time),
      old_time(rhs.old_time),
      new_data(std::move(rhs.new_data)),
      old_data(std::move(rhs.old_data))
{   
}

void
StateData::operator= (StateData const& rhs)
{
    m_factory.reset(rhs.m_factory->clone());
    desc = rhs.desc;
    domain = rhs.domain;
    grids = rhs.grids;
    dmap = rhs.dmap;
    new_time = rhs.new_time;
    old_time = rhs.old_time;
    new_data.reset(new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(), MFInfo(), *m_factory));
    MultiFab::Copy(*new_data, *rhs.new_data, 0, 0, desc->nComp(),desc->nExtra());
    if (rhs.old_data) {
        old_data.reset(new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(), MFInfo(), *m_factory));
        MultiFab::Copy(*old_data, *rhs.old_data, 0, 0, desc->nComp(),desc->nExtra());
    } else {
        old_data.reset();
    }
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

    new_data.reset(new MultiFab(grids,dmap,ncomp,desc->nExtra(), MFInfo(), *m_factory));
    old_data.reset();
}

void
StateData::copyOld (const StateData& state)
{
    const MultiFab& MF = state.oldData();
    
    int nc = MF.nComp();
    int ng = MF.nGrow();
    
    BL_ASSERT(nc == (*old_data).nComp());
    BL_ASSERT(ng == (*old_data).nGrow());
    
    MultiFab::Copy(*old_data, state.oldData(), 0, 0, nc, ng);
    
    old_time = state.old_time;
}

void
StateData::copyNew (const StateData& state)
{
    const MultiFab& MF = state.newData();

    int nc = MF.nComp();
    int ng = MF.nGrow();
    
    BL_ASSERT(nc == (*new_data).nComp());
    BL_ASSERT(ng == (*new_data).nGrow());
    
    MultiFab::Copy(*new_data, state.newData(), 0, 0, nc, ng);

    new_time = state.new_time;
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

    new_data.reset(new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(),
                                MFInfo(), *m_factory));
    old_data.reset();
    if (nsets == 2) {
        old_data.reset(new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(),
                                    MFInfo(), *m_factory));
    }
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
	whichMF = new_data.get();
      } else if(ns == 2) {
	whichMF = old_data.get();
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
    old_data.reset();
    new_data.reset(new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(), MFInfo(), *m_factory));
    new_data->setVal(0.);
}

StateData::~StateData()
{
    desc = nullptr;
}

void
StateData::allocOldData ()
{
    if (old_data == nullptr)
    {
        old_data.reset(new MultiFab(grids,dmap,desc->nComp(),desc->nExtra(), MFInfo(), *m_factory));
    }
}

BCRec
StateData::getBC (int comp, int i) const noexcept
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
StateData::replaceOldData (MultiFab&& mf)
{
    old_data.reset(new MultiFab(std::move(mf)));
}

// This version does NOT delete the replaced data.

void
StateData::replaceOldData (StateData& s)
{
    std::swap(old_data, s.old_data);
}

void
StateData::replaceNewData (MultiFab&& mf)
{
    new_data.reset(new MultiFab(std::move(mf)));
}

// This version does NOT delete the replaced data.

void
StateData::replaceNewData (StateData& s)
{
    std::swap(new_data, s.new_data);
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
    BL_PROFILE("StateData::FillBoundary(dx)");
    BL_ASSERT(dest.box().ixType() == desc->getType());
   
    if (domain.contains(dest.box())) return;

    const Box& bx  = dest.box();
    const int* dlo = dest.loVect();
    const int* dhi = dest.hiVect();
    const int* plo = domain.loVect();
    const int* phi = domain.hiVect();

    Vector<int> bcrs;

    Real xlo[AMREX_SPACEDIM];
    BCRec bcr;
    const Real* problo = prob_domain.lo();

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = problo[i] + dx[i]*(dlo[i]-plo[i]);
    }
    for (int i = 0; i < num_comp; )
    {
        const int dc  = dest_comp+i;
        const int sc  = src_comp+i;
        Real*     dat = dest.dataPtr(dc);

        if (desc->master(sc))
        {
            const int groupsize = desc->groupsize(sc);

            BL_ASSERT(groupsize != 0);

            if (groupsize+i <= num_comp)
            {
                //
                // Can do the whole group at once.
                //
                bcrs.resize(2*AMREX_SPACEDIM*groupsize);

                int* bci  = bcrs.dataPtr();

                for (int j = 0; j < groupsize; j++)
                {
                    amrex::setBC(bx,domain,desc->getBC(sc+j),bcr);

                    const int* bc = bcr.vect();

                    for (int k = 0; k < 2*AMREX_SPACEDIM; k++)
                        bci[k] = bc[k];

                    bci += 2*AMREX_SPACEDIM;
                }
                //
                // Use the "group" boundary fill routine.
                //
		desc->bndryFill(sc)(dat,dlo,dhi,plo,phi,dx,xlo,&time,bcrs.dataPtr(),groupsize);
                i += groupsize;
            }
            else
            {
                amrex::setBC(bx,domain,desc->getBC(sc),bcr);
                desc->bndryFill(sc)(dat,dlo,dhi,plo,phi,dx,xlo,&time,bcr.vect());
                i++;
            }
        }
        else
        {
            amrex::setBC(bx,domain,desc->getBC(sc),bcr);
            desc->bndryFill(sc)(dat,dlo,dhi,plo,phi,dx,xlo,&time,bcr.vect());
            i++;
        }
    }

#ifdef AMREX_USE_CUDA
    // Add a synchronize here in case the user code launched kernels
    // to handle the boundary fills.
    AMREX_GPU_SAFE_CALL(cudaDeviceSynchronize());
#endif
}

void
StateData::FillBoundary (Box const&      bx,
                         FArrayBox&      dest,
                         Real            time,
                         const Geometry& geom,
                         int             dest_comp,
                         int             src_comp,
                         int             num_comp)
{
    BL_PROFILE("StateData::FillBoundary(geom)");
    BL_ASSERT(bx.ixType() == desc->getType());
   
    if (domain.contains(bx)) return;

    Vector<BCRec> bcr(num_comp);

    for (int i = 0; i < num_comp; )
    {
        const int dc  = dest_comp+i;
        const int sc  = src_comp+i;

        if (desc->master(sc))
        {
            const int groupsize = desc->groupsize(sc);

            BL_ASSERT(groupsize != 0);

            if (groupsize+i <= num_comp)
            {
                for (int j = 0; j < groupsize; j++)
                {
                    amrex::setBC(bx,domain,desc->getBC(sc+j),bcr[j]);
                }
                //
                // Use the "group" boundary fill routine.
                //
		desc->bndryFill(sc)(bx,dest,dc,groupsize,geom,time,bcr,0,sc);
                i += groupsize;
            }
            else
            {
                amrex::setBC(bx,domain,desc->getBC(sc),bcr[0]);
                desc->bndryFill(sc)(bx,dest,dc,1,geom,time,bcr,0,sc);
                i++;
            }
        }
        else
        {
            amrex::setBC(bx,domain,desc->getBC(sc),bcr[0]);
            desc->bndryFill(sc)(bx,dest,dc,1,geom,time,bcr,0,sc);
            i++;
        }
    }
}

void
StateData::RegisterData (MultiFabCopyDescriptor& multiFabCopyDesc,
                         Vector<MultiFabId>&      mfid)
{
    mfid.resize(2);
    mfid[MFNEWDATA] = multiFabCopyDesc.RegisterFabArray(new_data.get());
    mfid[MFOLDDATA] = multiFabCopyDesc.RegisterFabArray(old_data.get());
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
        if (old_data == nullptr)
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
        else if (old_data != nullptr        &&
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
        if (old_data == nullptr)
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
        else if (old_data != nullptr        &&
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
	BL_ASSERT(new_data != nullptr);
        if (old_data == nullptr)
        {
	    data.push_back(new_data.get());
	    datatime.push_back(new_time.start);
        }
        else
        {
	    const Real teps = (new_time.start - old_time.start)*1.e-3;
	    if (time > new_time.start-teps && time < new_time.start+teps) {
		data.push_back(new_data.get());
		datatime.push_back(new_time.start);
	    } else if (time > old_time.start-teps && time < old_time.start+teps) {
	    	    data.push_back(old_data.get());
		    datatime.push_back(old_time.start);
	    } else {
		data.push_back(old_data.get());
		data.push_back(new_data.get());
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
	    data.push_back(new_data.get());
	    datatime.push_back(time);
        }
        else if (old_data != nullptr        &&
                 time > old_time.start-teps &&
                 time < old_time.stop+teps)
        {
	    data.push_back(old_data.get());
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

    if (dump_old == true && old_data == nullptr)
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
StateDataPhysBCFunct::FillBoundary (MultiFab& mf, int dest_comp, int num_comp, Real time, int /*bccomp*/)
{
    BL_PROFILE("StateDataPhysBCFunct::FillBoundary");

    const Box&     domain      = statedata->getDomain();
    const int*     domainlo    = domain.loVect();
    const int*     domainhi    = domain.hiVect();
    const Real*    dx          = geom.CellSize();
    const RealBox& prob_domain = geom.ProbDomain();

    bool has_bndryfunc_fab = statedata->desc->hasBndryFuncFab();
    bool run_on_gpu = statedata->desc->RunOnGPU() && Gpu::inLaunchRegion();

#if defined(AMREX_CRSEGRNDOMP) || (!defined(AMREX_XSDK) && defined(CRSEGRNDOMP))
#ifdef _OPENMP
#pragma omp parallel if (!run_on_gpu)
#endif
#endif
    {
	FArrayBox tmp;

	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	    FArrayBox* dest = run_on_gpu ? mf.fabPtr(mfi) : &(mf[mfi]);
	    const Box& bx = mfi.fabbox();
    
	    bool has_phys_bc = false;
	    bool is_periodic = false;
	    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
		bool touch = bx.smallEnd(i) < domainlo[i] || bx.bigEnd(i) > domainhi[i];
		if (geom.isPeriodic(i)) {
		    is_periodic = is_periodic || touch;
		} else {
		    has_phys_bc = has_phys_bc || touch;
		}
	    }

	    if (has_phys_bc)
	    {
                if (has_bndryfunc_fab) {
                    statedata->FillBoundary(bx, *dest, time, geom, dest_comp, src_comp, num_comp);
                } else {
                    statedata->FillBoundary(*dest, time, dx, prob_domain, dest_comp, src_comp, num_comp);
                }
		
		if (is_periodic) // fix up corner
		{
		    Box GrownDomain = domain;
		    
		    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
		    {
			if (!geom.isPeriodic(dir))
			{
			    const int lo = domainlo[dir] - bx.smallEnd(dir);
			    const int hi = bx.bigEnd(dir) - domainhi[dir];
			    if (lo > 0) GrownDomain.growLo(dir,lo);
			    if (hi > 0) GrownDomain.growHi(dir,hi);
			}
		    }
		    
		    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
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
                            if (run_on_gpu)
                            {
                                AsyncFab asyncfab(lo_slab,num_comp);
                                FArrayBox* fab = asyncfab.fabPtr();
                                const int ishift = -domain.length(dir);
                                AMREX_LAUNCH_DEVICE_LAMBDA (lo_slab, tbx,
                                {
                                    const Box db = amrex::shift(tbx, dir, ishift);
                                    fab->copy(*dest, db, dest_comp, tbx, 0, num_comp);
                                });
                                if (has_bndryfunc_fab) {
                                    statedata->FillBoundary(lo_slab, *fab, time, geom, 0, src_comp, num_comp);
                                } else {
                                    statedata->FillBoundary(*fab, time, dx, prob_domain, 0, src_comp, num_comp);
                                }
                                AMREX_LAUNCH_DEVICE_LAMBDA (lo_slab, tbx,
                                {
                                    const Box db = amrex::shift(tbx, dir, ishift);
                                    dest->copy(*fab, tbx, 0, db, dest_comp, num_comp);
                                });
                            }
                            else
                            {
                                tmp.resize(lo_slab,num_comp);
                                const Box db = amrex::shift(lo_slab, dir, -domain.length(dir));
                                tmp.copy(*dest, db, dest_comp, lo_slab, 0, num_comp);
                                if (has_bndryfunc_fab) {
                                    statedata->FillBoundary(lo_slab, tmp, time, geom, 0, src_comp, num_comp);
                                } else {
                                    statedata->FillBoundary(tmp, time, dx, prob_domain, 0, src_comp, num_comp);
                                }
                                dest->copy(tmp, lo_slab, 0, db, dest_comp, num_comp);
                            }
			}
			
			if (hi_slab.ok())
			{
                            if (run_on_gpu)
                            {
                                AsyncFab asyncfab(hi_slab,num_comp);
                                FArrayBox* fab = asyncfab.fabPtr();
                                const int ishift = domain.length(dir);
                                AMREX_LAUNCH_DEVICE_LAMBDA (hi_slab, tbx,
                                {
                                    const Box db = amrex::shift(tbx, dir, ishift);
                                    fab->copy(*dest, db, dest_comp, tbx, 0, num_comp);
                                });
                                if (has_bndryfunc_fab) {
                                    statedata->FillBoundary(hi_slab, *fab, time, geom, 0, src_comp, num_comp);
                                } else {
                                    statedata->FillBoundary(*fab, time, dx, prob_domain, 0, src_comp, num_comp);
                                }
                                AMREX_LAUNCH_DEVICE_LAMBDA (hi_slab, tbx,
                                {
                                    const Box db = amrex::shift(tbx, dir, ishift);
                                    dest->copy(*fab, tbx, 0, db, dest_comp, num_comp);
                                });
                            }
                            else
                            {
                                tmp.resize(hi_slab,num_comp);
                                const Box db = amrex::shift(hi_slab, dir, domain.length(dir));
                                tmp.copy(*dest, db, dest_comp, hi_slab, 0, num_comp);
                                if (has_bndryfunc_fab) {
                                    statedata->FillBoundary(hi_slab, tmp, time, geom, 0, src_comp, num_comp);
                                } else {
                                    statedata->FillBoundary(tmp, time, dx, prob_domain, 0, src_comp, num_comp);
                                }
                                dest->copy(tmp, hi_slab, 0, db, dest_comp, num_comp);
                            }
			}
		    }
		}
	    }
	}
    }
}

}
