
#include <sstream>

#include <unistd.h>
#include <memory>
#include <limits>

#include <AMReX_AmrLevel.H>
#include <AMReX_Derive.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#endif

#ifdef USE_PERILLA
#include <WorkerThread.H>
#endif

namespace amrex {
#ifdef USE_PERILLA
using namespace perilla;
#endif

#ifdef AMREX_USE_EB
int AmrLevel::m_eb_basic_grow_cells = 5;
int AmrLevel::m_eb_volume_grow_cells = 4;
int AmrLevel::m_eb_full_grow_cells = 2;
EBSupport AmrLevel::m_eb_support_level = EBSupport::volume;
#endif

DescriptorList AmrLevel::desc_lst;
DeriveList     AmrLevel::derive_lst;

void
AmrLevel::postCoarseTimeStep (Real time)
{
    BL_ASSERT(level == 0);
    // sync up statedata time
    for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
	AmrLevel& amrlevel = parent->getLevel(lev);
	for (int i = 0; i < amrlevel.state.size(); ++i) {
	    amrlevel.state[i].syncNewTimeLevel(time);
	}
    }
}

void
AmrLevel::set_preferred_boundary_values (MultiFab& S,
                                         int       state_index,
                                         int       scomp,
                                         int       dcomp,
                                         int       ncomp,
                                         Real      time) const
{}

DeriveList&
AmrLevel::get_derive_lst () noexcept
{
    return derive_lst;
}

void
AmrLevel::manual_tags_placement (TagBoxArray&    tags,
                                 const Vector<IntVect>& bf_lev)
{}

AmrLevel::AmrLevel () noexcept
{

   BL_PROFILE("AmrLevel::AmrLevel()");
   parent = 0;
   level = -1;
   levelDirectoryCreated = false;
}

AmrLevel::AmrLevel (Amr&            papa,
                    int             lev,
                    const Geometry& level_geom,
                    const BoxArray& ba,
		    const DistributionMapping& dm,
                    Real            time)
    :
    geom(level_geom),
    grids(ba),
    dmap(dm)
{
    BL_PROFILE("AmrLevel::AmrLevel(dm)");
    level  = lev;
    parent = &papa;
    levelDirectoryCreated = false;

    fine_ratio = IntVect::TheUnitVector(); fine_ratio.scale(-1);
    crse_ratio = IntVect::TheUnitVector(); crse_ratio.scale(-1);

    if (level > 0)
    {
        crse_ratio = parent->refRatio(level-1);
    }
    if (level < parent->maxLevel())
    {
        fine_ratio = parent->refRatio(level);
    }

    state.resize(desc_lst.size());

#ifdef AMREX_USE_EB
    m_factory = makeEBFabFactory(geom, ba, dm,
                                 {m_eb_basic_grow_cells,
                                  m_eb_volume_grow_cells,
                                  m_eb_full_grow_cells},
                                 m_eb_support_level);
#else
    m_factory.reset(new FArrayBoxFactory());
#endif

    // Note that this creates a distribution map associated with grids.
    for (int i = 0; i < state.size(); i++)
    {
        state[i].define(geom.Domain(),
                        grids,
			dm,
                        desc_lst[i],
                        time,
                        parent->dtLevel(lev),
                        *m_factory);
    }

    if (parent->useFixedCoarseGrids()) constructAreaNotToTag();

    post_step_regrid = 0;

    finishConstructor();
}

void
AmrLevel::writePlotFile (const std::string& dir,
                         std::ostream&      os,
                         VisMF::How         how)
{
    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
    {
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
	{
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
	    {
                plot_var_map.push_back(std::pair<int,int>(typ,comp));
	    }
	}
    }

    std::vector<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();
    for (std::list<DeriveRec>::const_iterator it = dlist.begin();
	 it != dlist.end();
	 ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
        {
            derive_names.push_back(it->name());
	}
    }

    int n_data_items = plot_var_map.size() + derive_names.size();

    // get the time from the first State_Type
    // if the State_Type is ::Interval, this will get t^{n+1/2} instead of t^n
    Real cur_time = state[0].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            amrex::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables
	//
	for (i =0; i < static_cast<int>(plot_var_map.size()); i++)
        {
	    int typ = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

        // derived
        for (auto const& dname : derive_names) {
            os << derive_lst.get(dname)->variableName(0) << '\n';
        }

        os << AMREX_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < AMREX_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < AMREX_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < AMREX_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string sLevel = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if ( ! FullPath.empty() && FullPath[FullPath.size()-1] != '/')
    {
        FullPath += '/';
    }
    FullPath += sLevel;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if ( ! levelDirectoryCreated) {
      if (ParallelDescriptor::IOProcessor()) {
        if ( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
            amrex::CreateDirectoryFailed(FullPath);
	}
      }
      // Force other processors to wait until directory is built.
      ParallelDescriptor::Barrier();
    }

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (n = 0; n < AMREX_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = sLevel;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    int       cnt   = 0;
    const int nGrow = 0;
    MultiFab  plotMF(grids,dmap,n_data_items,nGrow,MFInfo(),Factory());
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < static_cast<int>(plot_var_map.size()); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
	cnt++;
    }

    // derived
    if (derive_names.size() > 0)
    {
	for (auto const& dname : derive_names)
	{
            derive(dname, cur_time, plotMF, cnt);
	    cnt++;
	}
    }

    amrex::prefetchToHost(plotMF);

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);

    amrex::prefetchToDevice(plotMF);

    levelDirectoryCreated = false;  // ---- now that the plotfile is finished
}


void
AmrLevel::writePlotFilePre (const std::string& dir,
                            std::ostream&      os)
{
}


void
AmrLevel::writePlotFilePost (const std::string& dir,
                             std::ostream&      os)
{
}


void
AmrLevel::restart (Amr&          papa,
                   std::istream& is,
		   bool          bReadSpecial)
{
    BL_PROFILE("AmrLevel::restart()");
    parent = &papa;

    is >> level;
    is >> geom;

    fine_ratio = IntVect::TheUnitVector(); fine_ratio.scale(-1);
    crse_ratio = IntVect::TheUnitVector(); crse_ratio.scale(-1);

    if (level > 0)
    {
        crse_ratio = parent->refRatio(level-1);
    }
    if (level < parent->maxLevel())
    {
        fine_ratio = parent->refRatio(level);
    }

    if (bReadSpecial)
    {
        amrex::readBoxArray(grids, is, bReadSpecial);
    }
    else
    {
        grids.readFrom(is);
    }

    int nstate;
    is >> nstate;
    int ndesc = desc_lst.size();

    Vector<int> state_in_checkpoint(ndesc, 1);
    if (ndesc > nstate) {
	set_state_in_checkpoint(state_in_checkpoint);
    } else {
	BL_ASSERT(nstate == ndesc);
    }

    dmap.define(grids);

    parent->SetBoxArray(level, grids);
    parent->SetDistributionMap(level, dmap);

#ifdef AMREX_USE_EB
    m_factory = makeEBFabFactory(geom, grids, dmap,
                                 {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                                 m_eb_support_level);
#else
    m_factory.reset(new FArrayBoxFactory());
#endif

    state.resize(ndesc);
    for (int i = 0; i < ndesc; ++i)
    {
	if (state_in_checkpoint[i]) {
	    state[i].restart(is, geom.Domain(), grids, dmap, *m_factory,
			     desc_lst[i], papa.theRestartFile());
	}
    }
 
    if (parent->useFixedCoarseGrids()) constructAreaNotToTag();

    post_step_regrid = 0;

    finishConstructor();
}

void
AmrLevel::set_state_in_checkpoint (Vector<int>& state_in_checkpoint)
{
    amrex::Error("Class derived AmrLevel has to handle this!");
}

void
AmrLevel::finishConstructor () {}

void
AmrLevel::setTimeLevel (Real time,
                        Real dt_old,
                        Real dt_new)
{
    for (int k = 0; k < desc_lst.size(); k++)
    {
        state[k].setTimeLevel(time,dt_old,dt_new);
    }
}

bool
AmrLevel::isStateVariable (const std::string& name, int& typ, int& n)
{
    for (typ = 0; typ < desc_lst.size(); typ++)
    {
        const StateDescriptor& desc = desc_lst[typ];

        for (n = 0; n < desc.nComp(); n++)
        {
            if (desc.name(n) == name)
                return true;
        }
    }
    return false;
}

long
AmrLevel::countCells () const noexcept
{
    return grids.numPts();
}

void
AmrLevel::checkPoint (const std::string& dir,
                      std::ostream&      os,
                      VisMF::How         how,
                      bool               dump_old)
{
    BL_PROFILE("AmrLevel::checkPoint()");
    int ndesc = desc_lst.size(), i;
    //
    // Build directory to hold the MultiFabs in the StateData at this level.
    // The directory is relative the the directory containing the Header file.
    //
    std::string LevelDir, FullPath;
    LevelDirectoryNames(dir, LevelDir, FullPath);
    if( ! levelDirectoryCreated) {
      CreateLevelDirectory(dir);
      // ---- Force other processors to wait until directory is built.
      ParallelDescriptor::Barrier("AmrLevel::checkPoint::dir");
    }

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << '\n' << geom  << '\n';
        grids.writeOn(os);
        os << ndesc << '\n';
    }
    //
    // Output state data.
    //

    for (i = 0; i < ndesc; i++)
    {
        //
        // Now build the full relative pathname of the StateData.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        std::string PathNameInHdr = amrex::Concatenate(LevelDir + "/SD_", i, 1);
        std::string FullPathName  = amrex::Concatenate(FullPath + "/SD_", i, 1);

        state[i].checkPoint(PathNameInHdr, FullPathName, os, how, dump_old);
    }

    levelDirectoryCreated = false;  // ---- now that the checkpoint is finished
}


void
AmrLevel::checkPointPre (const std::string& dir,
                         std::ostream&      os)
{
    BL_PROFILE("AmrLevel::checkPointPre()");
}


void
AmrLevel::checkPointPost (const std::string& dir,
                          std::ostream&      os)
{
    BL_PROFILE("AmrLevel::checkPointPost()");
}


AmrLevel::~AmrLevel ()
{
    parent = 0;
}

void
AmrLevel::allocOldData ()
{
    for (int i = 0; i < desc_lst.size(); i++)
    {
        state[i].allocOldData();
    }
}

void
AmrLevel::removeOldData ()
{
    for (int i = 0; i < desc_lst.size(); i++)
    {
        state[i].removeOldData();
    }
}

void
AmrLevel::reset ()
{
    for (int i = 0; i < desc_lst.size(); i++)
    {
        state[i].reset();
    }
}

MultiFab&
AmrLevel::get_data (int  state_indx, Real time) noexcept
{
    const Real old_time = state[state_indx].prevTime();
    const Real new_time = state[state_indx].curTime();
    const Real eps = 0.001*(new_time - old_time);

    if (time > old_time-eps && time < old_time+eps)
    {
        return get_old_data(state_indx);
    }
    else if (time > new_time-eps && time < new_time+eps)
    {
        return get_new_data(state_indx);
    }

    amrex::Error("get_data: invalid time");
    static MultiFab bogus;
    return bogus;
}

const BoxArray&
AmrLevel::getEdgeBoxArray (int dir) const noexcept
{
    BL_ASSERT(dir >=0 && dir < AMREX_SPACEDIM);
    if (edge_grids[dir].empty()) {
	edge_grids[dir] = grids;
	edge_grids[dir].surroundingNodes(dir);
    }
    return edge_grids[dir];
}

const BoxArray&
AmrLevel::getNodalBoxArray () const noexcept
{
    if (nodal_grids.empty()) {
	nodal_grids = grids;
	nodal_grids.surroundingNodes();
    }
    return nodal_grids;
}

void
AmrLevel::setPhysBoundaryValues (FArrayBox& dest,
                                 int        state_indx,
                                 Real       time,
                                 int        dest_comp,
                                 int        src_comp,
                                 int        num_comp)
{
    state[state_indx].FillBoundary(dest,time,geom.CellSize(),
                                   geom.ProbDomain(),dest_comp,src_comp,num_comp);
}

FillPatchIteratorHelper::FillPatchIteratorHelper (AmrLevel& amrlevel,
                                                  MultiFab& leveldata)
    :
    m_amrlevel(amrlevel),
    m_leveldata(leveldata),
    m_mfid(m_amrlevel.level+1)
{}

FillPatchIterator::FillPatchIterator (AmrLevel& amrlevel,
                                      MultiFab& leveldata)
    :
    MFIter(leveldata),
    m_amrlevel(amrlevel),
    m_leveldata(leveldata),
    m_ncomp(0)
{}

FillPatchIteratorHelper::FillPatchIteratorHelper (AmrLevel&     amrlevel,
                                                  MultiFab&     leveldata,
                                                  int           boxGrow,
                                                  Real          time,
                                                  int           index,
                                                  int           scomp,
                                                  int           ncomp,
                                                  Interpolater* mapper)
    :
    m_amrlevel(amrlevel),
    m_leveldata(leveldata),
    m_mfid(m_amrlevel.level+1),
    m_time(time),
    m_growsize(boxGrow),
    m_index(index),
    m_scomp(scomp),
    m_ncomp(ncomp)
{
    Initialize(boxGrow,time,index,scomp,ncomp,mapper);
}

FillPatchIterator::FillPatchIterator (AmrLevel& amrlevel,
                                      MultiFab& leveldata,
                                      int       boxGrow,
                                      Real      time,
                                      int       idx,
                                      int       scomp,
                                      int       ncomp)
    :
    MFIter(leveldata),
    m_amrlevel(amrlevel),
    m_leveldata(leveldata),
    m_ncomp(ncomp)
{
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(AmrLevel::desc_lst[idx].inRange(scomp,ncomp));
    BL_ASSERT(0 <= idx && idx < AmrLevel::desc_lst.size());

    Initialize(boxGrow,time,idx,scomp,ncomp);

#ifdef BL_USE_TEAM
    ParallelDescriptor::MyTeam().MemoryBarrier();
#endif
}

static
bool
NeedToTouchUpPhysCorners (const Geometry& geom)
{
    return geom.isAnyPeriodic() && !geom.isAllPeriodic();
}

void
FillPatchIteratorHelper::Initialize (int           boxGrow,
                                     Real          time,
                                     int           idx,
                                     int           scomp,
                                     int           ncomp,
                                     Interpolater* mapper)
{
    BL_PROFILE("FillPatchIteratorHelper::Initialize()");

    BL_ASSERT(mapper);
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(AmrLevel::desc_lst[idx].inRange(scomp,ncomp));
    BL_ASSERT(0 <= idx && idx < AmrLevel::desc_lst.size());

    m_map          = mapper;
    m_time         = time;
    m_growsize     = boxGrow;
    m_index        = idx;
    m_scomp        = scomp;
    m_ncomp        = ncomp;
    m_FixUpCorners = NeedToTouchUpPhysCorners(m_amrlevel.geom);

    const int         MyProc     = ParallelDescriptor::MyProc();
    auto&             amrLevels  = m_amrlevel.parent->getAmrLevels();
    const AmrLevel&   topLevel   = *amrLevels[m_amrlevel.level];
    const Box&        topPDomain = topLevel.state[m_index].getDomain();
    const IndexType&  boxType    = m_leveldata.boxArray().ixType();
    const bool        extrap     = AmrLevel::desc_lst[m_index].extrap();
    //
    // Check that the interpolaters are identical.
    //
    BL_ASSERT(AmrLevel::desc_lst[m_index].identicalInterps(scomp,ncomp));

    for (int l = 0; l <= m_amrlevel.level; ++l)
    {
        amrLevels[l]->state[m_index].RegisterData(m_mfcd, m_mfid[l]);
    }
    for (int i = 0, N = m_leveldata.boxArray().size(); i < N; ++i)
    {
        //
        // A couple typedefs we'll use in the next code segment.
        //
        typedef std::map<int,Vector<Vector<Box> > >::value_type IntAABoxMapValType;

        typedef std::map<int,Vector<Vector<Vector<FillBoxId> > > >::value_type IntAAAFBIDMapValType;

        if (m_leveldata.DistributionMap()[i] != MyProc) continue;
        //
        // Insert with a hint since the indices are ordered lowest to highest.
        //
        IntAAAFBIDMapValType v1(i,Vector<Vector<Vector<FillBoxId> > >());

        m_fbid.insert(m_fbid.end(),v1)->second.resize(m_amrlevel.level+1);

        IntAABoxMapValType v2(i,Vector<Vector<Box> >());

        m_fbox.insert(m_fbox.end(),v2)->second.resize(m_amrlevel.level+1);
        m_cbox.insert(m_cbox.end(),v2)->second.resize(m_amrlevel.level+1);

        m_ba.insert(m_ba.end(),std::map<int,Box>::value_type(i,amrex::grow(m_leveldata.boxArray()[i],m_growsize)));
    }

    BoxList        tempUnfillable(boxType);
    BoxList        unfillableThisLevel(boxType);
    Vector<Box>     unfilledThisLevel;
    Vector<Box>     crse_boxes;
    Vector<IntVect> pshifts(27);

    for (std::map<int,Box>::const_iterator it = m_ba.begin(), End = m_ba.end();
         it != End;
         ++it)
    {
        const int  bxidx = it->first;
        const Box& box   = it->second;

        unfilledThisLevel.clear();
        unfilledThisLevel.push_back(box);

        if (!topPDomain.contains(box))
        {
            unfilledThisLevel.back() &= topPDomain;

            if (topLevel.geom.isAnyPeriodic())
            {
                //
                // May need to add additional unique pieces of valid region
                // in order to do periodic copies into ghost cells.
                //
                topLevel.geom.periodicShift(topPDomain,box,pshifts);

                for (const auto& iv : pshifts)
                {
                    Box shbox = box + iv;
                    shbox    &= topPDomain;

                    if (boxType.nodeCentered())
                    {
                        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
                        {
                            if (iv[dir] > 0)
                            {
                                shbox.growHi(dir,-1);
                            }
                            else if (iv[dir] < 0)
                            {
                                shbox.growLo(dir,-1);
                            }
                        }
                    }

                    if (shbox.ok())
                    {
                        BoxList bl = amrex::boxDiff(shbox,box);

                        unfilledThisLevel.insert(unfilledThisLevel.end(), bl.begin(), bl.end());
                    }
                }
            }
        }

	// cells outside physical boundaries are not included in unfilledThisLevel

        bool Done = false;

        Vector< Vector<Box> >&                TheCrseBoxes = m_cbox[bxidx];
        Vector< Vector<Box> >&                TheFineBoxes = m_fbox[bxidx];
        Vector< Vector< Vector<FillBoxId> > >& TheFBIDs     = m_fbid[bxidx];

        for (int l = m_amrlevel.level; l >= 0 && !Done; --l)
        {
            unfillableThisLevel.clear();

            AmrLevel&       theAmrLevel = *amrLevels[l];
            StateData&      theState    = theAmrLevel.state[m_index];
            const Box&      thePDomain  = theState.getDomain();
            const Geometry& theGeom     = theAmrLevel.geom;
            const bool      is_periodic = theGeom.isAnyPeriodic();
            const IntVect&  fine_ratio  = theAmrLevel.fine_ratio;
            Vector<Box>&     FineBoxes   = TheFineBoxes[l];
            //
            // These are the boxes on this level contained in thePDomain
            // that need to be filled in order to directly fill at the
            // highest level or to interpolate up to the next higher level.
            //
            FineBoxes = unfilledThisLevel;
            //
            // Now build coarse boxes needed to interpolate to fine.
            //
            // If we're periodic and we're not at the finest level, we may
            // need to get some additional data at this level in order to
            // properly fill the CoarseBox()d versions of the fineboxes.
            //
            crse_boxes.clear();

            for (const auto& fbx : FineBoxes)
            {
                crse_boxes.push_back(fbx);

                if (l != m_amrlevel.level)
                {
                    const Box& cbox = m_map->CoarseBox(fbx,fine_ratio);

		    crse_boxes.back() = cbox;

                    if (is_periodic && !thePDomain.contains(cbox))
                    {
                        theGeom.periodicShift(thePDomain,cbox,pshifts);

                        for (const auto& iv : pshifts)
                        {
                            Box shbox = cbox + iv;
                            shbox    &= thePDomain;

                            if (boxType.nodeCentered())
                            {
                                for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
                                {
                                    if (iv[dir] > 0)
                                    {
                                        shbox.growHi(dir,-1);
                                    }
                                    else if (iv[dir] < 0)
                                    {
                                        shbox.growLo(dir,-1);
                                    }
                                }
                            }

                            if (shbox.ok())
                            {
                                crse_boxes.push_back(shbox);
                            }
                        }
                    }
                }
            }

            Vector< Vector<FillBoxId> >& FBIDs     = TheFBIDs[l];
            Vector<Box>&                CrseBoxes = TheCrseBoxes[l];

            FBIDs.resize(crse_boxes.size());
            CrseBoxes.resize(crse_boxes.size());
            //
            // Now attempt to get as much coarse data as possible.
            //
            for (int i = 0, M = CrseBoxes.size(); i < M; i++)
            {
                BL_ASSERT(tempUnfillable.isEmpty());

                CrseBoxes[i] = crse_boxes[i];

                BL_ASSERT(CrseBoxes[i].intersects(thePDomain));

                theState.InterpAddBox(m_mfcd,
				      m_mfid[l],
				      &tempUnfillable,
				      FBIDs[i],
				      CrseBoxes[i],
				      m_time,
				      m_scomp,
				      0,
				      m_ncomp,
				      extrap);

                unfillableThisLevel.catenate(tempUnfillable);
            }

            unfillableThisLevel.intersect(thePDomain);

            if (unfillableThisLevel.isEmpty())
            {
                Done = true;
            }
            else
            {
                unfilledThisLevel.clear();

                unfilledThisLevel.insert(unfilledThisLevel.end(),
                                         unfillableThisLevel.begin(),
                                         unfillableThisLevel.end());
            }
        }
    }

    m_mfcd.CollectData();
}

void
FillPatchIterator::Initialize (int  boxGrow,
                               Real time,
                               int  idx,
                               int  scomp,
                               int  ncomp)
{
    BL_PROFILE("FillPatchIterator::Initialize");

    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(0 <= idx && idx < AmrLevel::desc_lst.size());

    const StateDescriptor& desc = AmrLevel::desc_lst[idx];

    m_ncomp = ncomp;
    m_range = desc.sameInterps(scomp,ncomp);

    m_fabs.define(m_leveldata.boxArray(),m_leveldata.DistributionMap(),
		  m_ncomp,boxGrow,MFInfo(),m_leveldata.Factory());

    const Geometry& geom = m_amrlevel.Geom();

    m_fabs.setDomainBndry(std::numeric_limits<Real>::quiet_NaN(), geom);

    const IndexType& boxType = m_leveldata.boxArray().ixType();
    const int level = m_amrlevel.level;

    for (int i = 0, DComp = 0; i < static_cast<int>(m_range.size()); i++)
    {
        const int SComp = m_range[i].first;
        const int NComp = m_range[i].second;

	if (level == 0)
	{
	    FillFromLevel0(time, idx, SComp, DComp, NComp);
	}
	else
	{
	    if (level == 1 || 
		amrex::ProperlyNested(m_amrlevel.crse_ratio,
				       m_amrlevel.parent->blockingFactor(m_amrlevel.level),
				       boxGrow, boxType, desc.interp(SComp)))
	    {
		FillFromTwoLevels(time, idx, SComp, DComp, NComp);
	    } else {

#ifdef AMREX_USE_EB
                amrex::Abort("Grids must be properly nested for EB");
#endif

		static bool first = true;
		if (first) {
		    first = false;
		    if (ParallelDescriptor::IOProcessor() && amrex::Verbose()) {
			IntVect new_blocking_factor = m_amrlevel.parent->blockingFactor(m_amrlevel.level);
                        new_blocking_factor *= 2;
			for (int j = 0; j < 10; ++j) {
			    if (amrex::ProperlyNested(m_amrlevel.crse_ratio,
						       new_blocking_factor,
						       boxGrow, boxType, desc.interp(SComp))) {
				break;
			    } else {
				new_blocking_factor *= 2;
			    }
			}
                        amrex::Print() << "WARNING: Grids are not properly nested.  We might have to use\n"
                                       << "         two coarse levels to do fillpatch.  Consider using\n";
                        if (new_blocking_factor < IntVect{AMREX_D_DECL(128,128,128)}) {
                            amrex::Print() << "         amr.blocking_factor=" << new_blocking_factor << "\n";
                        } else {
                            amrex::Print() << "         larger amr.blocking_factor.\n";
                        }
		    }
		}

		FillPatchIteratorHelper* fph = 0;
		fph = new FillPatchIteratorHelper(m_amrlevel,
						  m_leveldata,
						  boxGrow,
						  time,
						  idx,
						  SComp,
						  NComp,
						  desc.interp(SComp));
	
#if defined(AMREX_CRSEGRNDOMP) || (!defined(AMREX_XSDK) && defined(CRSEGRNDOMP))
#ifdef _OPENMP
#pragma omp parallel
#endif
#endif
		for (MFIter mfi(m_fabs); mfi.isValid(); ++mfi)
		{
		    fph->fill(*m_fabs.fabPtr(mfi),DComp,mfi.index());
		}
		
		delete fph;
	    }
	}

        DComp += NComp;
    }
    //
    // Call hack to touch up fillPatched data.
    //
    m_amrlevel.set_preferred_boundary_values(m_fabs,
                                             idx,
                                             scomp,
                                             0,
                                             ncomp,
                                             time);
}

void
FillPatchIterator::FillFromLevel0 (Real time, int idx, int scomp, int dcomp, int ncomp)
{
    BL_ASSERT(m_amrlevel.level == 0);

    StateData& statedata = m_amrlevel.state[idx];

    Vector<MultiFab*> smf;
    Vector<Real> stime;
    statedata.getData(smf,stime,time);

    const Geometry& geom = m_amrlevel.geom;

    StateDataPhysBCFunct physbcf(statedata,scomp,geom);

    amrex::FillPatchSingleLevel (m_fabs, time, smf, stime, scomp, dcomp, ncomp, geom, physbcf, scomp);
}

void
FillPatchIterator::FillFromTwoLevels (Real time, int idx, int scomp, int dcomp, int ncomp)
{
    int ilev_fine = m_amrlevel.level;
    int ilev_crse = ilev_fine-1;

    BL_ASSERT(ilev_crse >= 0);

    AmrLevel& fine_level = m_amrlevel;
    AmrLevel& crse_level = m_amrlevel.parent->getLevel(ilev_crse);

    const Geometry& geom_fine = fine_level.geom;
    const Geometry& geom_crse = crse_level.geom;
    
    Vector<MultiFab*> smf_crse;
    Vector<Real> stime_crse;
    StateData& statedata_crse = crse_level.state[idx];
    statedata_crse.getData(smf_crse,stime_crse,time);
    StateDataPhysBCFunct physbcf_crse(statedata_crse,scomp,geom_crse);

    Vector<MultiFab*> smf_fine;
    Vector<Real> stime_fine;
    StateData& statedata_fine = fine_level.state[idx];
    statedata_fine.getData(smf_fine,stime_fine,time);
    StateDataPhysBCFunct physbcf_fine(statedata_fine,scomp,geom_fine);

    const StateDescriptor& desc = AmrLevel::desc_lst[idx];

    amrex::FillPatchTwoLevels(m_fabs, time, 
                              smf_crse, stime_crse, 
                              smf_fine, stime_fine,
                              scomp, dcomp, ncomp, 
                              geom_crse, geom_fine,
                              physbcf_crse, scomp,
                              physbcf_fine, scomp,
                              crse_level.fineRatio(), 
                              desc.interp(scomp),
                              desc.getBCs(),scomp);
}

static
bool
HasPhysBndry (const Box&      b,
              const Box&      dmn,
              const Geometry& geom)
{
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (!geom.isPeriodic(i))
        {
            if (b.smallEnd(i) < dmn.smallEnd(i) || b.bigEnd(i) > dmn.bigEnd(i))
            {
                return true;
            }
        }
    }

    return false;
}

static
void
FixUpPhysCorners (FArrayBox&      fab,
                  AmrLevel&       TheLevel,
                  int             state_indx,
                  Real            time,
                  int             scomp,
                  int             dcomp,
                  int             ncomp)
{
    StateData&      TheState   = TheLevel.get_state_data(state_indx);
    const Geometry& TheGeom    = TheLevel.Geom();
    const Box&      ProbDomain = TheState.getDomain();

    if (!HasPhysBndry(fab.box(),ProbDomain,TheGeom)) return;

    FArrayBox tmp;

    Box GrownDomain = ProbDomain;

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        if (!TheGeom.isPeriodic(dir))
        {
            const int lo = ProbDomain.smallEnd(dir) - fab.box().smallEnd(dir);
            const int hi = fab.box().bigEnd(dir)    - ProbDomain.bigEnd(dir);
            if (lo > 0) GrownDomain.growLo(dir,lo);
            if (hi > 0) GrownDomain.growHi(dir,hi);
        }
    }

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        if (!TheGeom.isPeriodic(dir)) continue;

        Box lo_slab = fab.box();
        Box hi_slab = fab.box();
        lo_slab.shift(dir, ProbDomain.length(dir));
        hi_slab.shift(dir,-ProbDomain.length(dir));
        lo_slab &= GrownDomain;
        hi_slab &= GrownDomain;

        if (lo_slab.ok())
        {
            lo_slab.shift(dir,-ProbDomain.length(dir));

            BL_ASSERT(fab.box().contains(lo_slab));
            BL_ASSERT(HasPhysBndry(lo_slab,ProbDomain,TheGeom));

            tmp.resize(lo_slab,ncomp);
            tmp.copy(fab,dcomp,0,ncomp);
            tmp.shift(dir,ProbDomain.length(dir));
            TheLevel.setPhysBoundaryValues(tmp,
                                           state_indx,
                                           time,
                                           0,
                                           scomp,
                                           ncomp);
            tmp.shift(dir,-ProbDomain.length(dir));
            fab.copy(tmp,0,dcomp,ncomp);
        }

        if (hi_slab.ok())
        {
            hi_slab.shift(dir,ProbDomain.length(dir));

            BL_ASSERT(fab.box().contains(hi_slab));
            BL_ASSERT(HasPhysBndry(hi_slab,ProbDomain,TheGeom));

            tmp.resize(hi_slab,ncomp);
            tmp.copy(fab,dcomp,0,ncomp);
            tmp.shift(dir,-ProbDomain.length(dir));
            TheLevel.setPhysBoundaryValues(tmp,
                                           state_indx,
                                           time,
                                           0,
                                           scomp,
                                           ncomp);
            tmp.shift(dir,ProbDomain.length(dir));
            fab.copy(tmp,0,dcomp,ncomp);
        }
    }
}

void
FillPatchIteratorHelper::fill (FArrayBox& fab,
                               int        dcomp,
                               int        idx)
{
    BL_PROFILE("FillPatchIteratorHelper::fill()");

    BL_ASSERT(fab.box() == m_ba[idx]);
    BL_ASSERT(fab.nComp() >= dcomp + m_ncomp);

    Vector< Vector<std::unique_ptr<FArrayBox> > > cfab(m_amrlevel.level+1);
    Vector< Vector<Box> >&                TheCrseBoxes = m_cbox[idx];
    Vector< Vector<Box> >&                TheFineBoxes = m_fbox[idx];
    Vector< Vector< Vector<FillBoxId> > >& TheFBIDs     = m_fbid[idx];
    const bool                          extrap       = AmrLevel::desc_lst[m_index].extrap();
    auto&                               amrLevels    = m_amrlevel.parent->getAmrLevels();
    //
    // Build all coarse fabs from which we'll interpolate and
    // fill them with coarse data as best we can.
    //
    for (int l = 0; l <= m_amrlevel.level; l++)
    {
        StateData&                       TheState  = amrLevels[l]->state[m_index];
        const Vector<Box>&                CrseBoxes = TheCrseBoxes[l];
        auto&                            CrseFabs  = cfab[l];
        const Vector< Vector<FillBoxId> >& FBIDs     = TheFBIDs[l];
        const int                        NC        = CrseBoxes.size();

        CrseFabs.resize(NC);

        Box domain_box = amrex::convert(amrLevels[l]->Geom().Domain(), fab.box().ixType());
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (Geometry::isPeriodic(idim)) {
                int n = domain_box.length(idim);
                domain_box.grow(idim, n);
            }
        }

        for (int i = 0; i < NC; i++)
        {
            BL_ASSERT(CrseBoxes[i].ok());
            CrseFabs[i].reset(new FArrayBox(CrseBoxes[i],m_ncomp));
            CrseFabs[i]->setComplement(std::numeric_limits<Real>::quiet_NaN(), domain_box, 0, m_ncomp);
	}

        for (int i = 0; i < NC; i++)
        {
            //
            // Set to special value we'll later check
            // to ensure we've filled the FABs at the coarse level.
            //
            TheState.InterpFillFab(m_mfcd,
				   m_mfid[l],
				   FBIDs[i],
				   *CrseFabs[i],
				   m_time,
				   0,
				   0,
				   m_ncomp,
				   extrap);
        }
    }
    //
    // Now work from the bottom up interpolating to next higher level.
    //
    for (int l = 0; l < m_amrlevel.level; l++)
    {
        auto&              CrseFabs   = cfab[l];
        AmrLevel&          TheLevel   = *amrLevels[l];
        StateData&         TheState   = TheLevel.state[m_index];
        const Box&         ThePDomain = TheState.getDomain();
        const int          NC         = CrseFabs.size();

        if (TheLevel.geom.isAnyPeriodic())
        {
            //
            // Fill CrseFabs with periodic data in preparation for interp().
            //
            for (int i = 0; i < NC; i++)
            {
                FArrayBox& dstfab = *CrseFabs[i];

                if (ThePDomain.contains(dstfab.box())) continue;

                Vector<IntVect> pshifts(27);

                TheLevel.geom.periodicShift(ThePDomain,dstfab.box(),pshifts);

                for (const auto& iv : pshifts)
                {
                    Box fullsrcbox = dstfab.box() + iv;
                    fullsrcbox    &= ThePDomain;

                    for (int j = 0; j < NC; j++)
                    {
                        const FArrayBox& srcfab = *CrseFabs[j];
                        const Box&       srcbox = fullsrcbox & srcfab.box();

                        if (srcbox.ok())
                        {
                            const Box& dstbox = srcbox - iv;

                            dstfab.copy(srcfab,srcbox,0,dstbox,0,m_ncomp);
                        }
                    }
                }
            }
        }
        //
        // Set non-periodic BCs in coarse data -- what we interpolate with.
        // This MUST come after the periodic fill mumbo-jumbo.
        for (int i = 0; i < NC; ++i)
        {
            if ( ! ThePDomain.contains(CrseFabs[i]->box()))
            {
                TheLevel.setPhysBoundaryValues(*CrseFabs[i],
                                               m_index,
                                               m_time,
                                               0,
                                               m_scomp,
                                               m_ncomp);
            }
        }

        if (m_FixUpCorners)
        {
            for (int i = 0; i < NC; ++i)
            {
                FixUpPhysCorners(*CrseFabs[i],TheLevel,m_index,m_time,m_scomp,0,m_ncomp);
            }
        }
        //
        // Interpolate up to next level.
        //
        AmrLevel&           crseAmrLevel  = *amrLevels[l];
        AmrLevel&           fineAmrLevel  = *amrLevels[l+1];
        const IntVect&      fine_ratio    = crseAmrLevel.fine_ratio;
        const Vector<Box>&   FineBoxes     = TheFineBoxes[l];
        StateData&          fState        = fineAmrLevel.state[m_index];
        const Box&          fDomain       = fState.getDomain();
        auto&               FinerCrseFabs = cfab[l+1];
        const Vector<BCRec>& theBCs        = AmrLevel::desc_lst[m_index].getBCs();
        const int           NF            = FineBoxes.size();

        for (int ifine = 0; ifine < NF; ++ifine)
        {
            Vector<BCRec> bcr(m_ncomp);
            FArrayBox    finefab(FineBoxes[ifine],m_ncomp);
            FArrayBox    crsefab(m_map->CoarseBox(finefab.box(),fine_ratio),m_ncomp);
            //
            // Fill crsefab from m_cbox via copy on intersect.
            //
            for (int j = 0; j < NC; j++) {
                crsefab.copy(*CrseFabs[j]);
	    }
            //
            // Get boundary conditions for the fine patch.
            //
            amrex::setBC(finefab.box(),
                          fDomain,
                          m_scomp,
                          0,
                          m_ncomp,
                          theBCs,
                          bcr);
            //
            // Interpolate up to fine patch.
            //
            m_map->interp(crsefab,
                          0,
                          finefab,
                          0,
                          m_ncomp,
                          finefab.box(),
                          fine_ratio,
                          crseAmrLevel.geom,
                          fineAmrLevel.geom,
                          bcr,
                          m_scomp,
                          m_index);
            //
            // Copy intersect finefab into next level m_cboxes.
            //
	    for (int j = 0, K = FinerCrseFabs.size(); j < K; ++j) {
		FinerCrseFabs[j]->copy(finefab);
	    }
        }

        CrseFabs.clear();
    }
    //
    // Now for the finest level stuff.
    //
    StateData&         FineState      = m_amrlevel.state[m_index];
    const Box&         FineDomain     = FineState.getDomain();
    const Geometry&    FineGeom       = m_amrlevel.geom;
    auto&              FinestCrseFabs = cfab[m_amrlevel.level];
    //
    // Copy intersect coarse into destination fab.
    //
    for (int i = 0, N = FinestCrseFabs.size(); i < N; ++i) {
        fab.copy(*FinestCrseFabs[i],0,dcomp,m_ncomp);
    }

    if (FineGeom.isAnyPeriodic() && !FineDomain.contains(fab.box()))
    {
        Vector<IntVect> pshifts(27);

        FineGeom.periodicShift(FineDomain,fab.box(),pshifts);

        for (int i = 0, N = FinestCrseFabs.size(); i < N; i++)
        {
            for (const auto& iv : pshifts)
            {
                fab.shift(iv);

                Box src_dst = FinestCrseFabs[i]->box() & fab.box();
                src_dst    &= FineDomain;

                if (src_dst.ok())
                    fab.copy(*FinestCrseFabs[i],src_dst,0,src_dst,dcomp,m_ncomp);

                fab.shift(-iv);
            }
        }
    }
    //
    // No longer need coarse data at finest level.
    //
    FinestCrseFabs.clear();
    //
    // Final set of non-periodic BCs.
    //
    if (! FineState.getDomain().contains(fab.box()))
    {
        m_amrlevel.setPhysBoundaryValues(fab,
                                         m_index,
                                         m_time,
                                         dcomp,
                                         m_scomp,
                                         m_ncomp);
    }

    if (m_FixUpCorners)
    {
        FixUpPhysCorners(fab,m_amrlevel,m_index,m_time,m_scomp,dcomp,m_ncomp);
    }
}

FillPatchIteratorHelper::~FillPatchIteratorHelper () {}

FillPatchIterator::~FillPatchIterator () {
#ifdef USE_PERILLA
        while(regionList.size()){
          RegionGraph* tmp= regionList.front();
          delete tmp;
          regionList.pop_front();
        }

        while(mfList.size()){
          MultiFab *tmp= mfList.front();
          delete tmp;
          mfList.pop_front();
        }

        while(stateDataList.size()){
          StateDataPhysBCFunct *tmp= stateDataList.front();
          delete tmp;
          stateDataList.pop_front();
        }
#endif
    }

void
AmrLevel::FillCoarsePatch (MultiFab& mf,
                           int       dcomp,
                           Real      time,
                           int       idx,
                           int       scomp,
                           int       ncomp,
			   int       nghost)
{
    BL_PROFILE("AmrLevel::FillCoarsePatch()");

    //
    // Must fill this region on crse level and interpolate.
    //
    BL_ASSERT(level != 0);
    BL_ASSERT(ncomp <= (mf.nComp()-dcomp));
    BL_ASSERT(nghost <= mf.nGrow());
    BL_ASSERT(0 <= idx && idx < desc_lst.size());

    int                     DComp   = dcomp;
    const StateDescriptor&  desc    = desc_lst[idx];
    const Box&              pdomain = state[idx].getDomain();
    const BoxArray&         mf_BA   = mf.boxArray();
    const DistributionMapping& mf_DM = mf.DistributionMap();
    AmrLevel&               clev    = parent->getLevel(level-1);
    const Geometry&         cgeom   = clev.geom;

    Box domain_g = pdomain;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
	if (geom.isPeriodic(i)) {
	    domain_g.grow(i,nghost);
	}
    }

    std::vector< std::pair<int,int> > ranges  = desc.sameInterps(scomp,ncomp);

    BL_ASSERT(desc.inRange(scomp, ncomp));

    for (int i = 0; i < static_cast<int>(ranges.size()); i++)
    {
        const int     SComp  = ranges[i].first;
        const int     NComp  = ranges[i].second;
        Interpolater* mapper = desc.interp(SComp);

        BoxArray crseBA(mf_BA.size());
        
        for (int j = 0, N = crseBA.size(); j < N; ++j)
        {
            BL_ASSERT(mf_BA[j].ixType() == desc.getType());
	    const Box& bx = amrex::grow(mf_BA[j],nghost) & domain_g;
            crseBA.set(j,mapper->CoarseBox(bx, crse_ratio));
        }

#ifdef AMREX_USE_EB
        auto cfactory = makeEBFabFactory(cgeom, crseBA, mf_DM, {0,0,0}, EBSupport::basic);
        MultiFab crseMF(crseBA,mf_DM,NComp,0,MFInfo(),*cfactory);
#else
	MultiFab crseMF(crseBA,mf_DM,NComp,0);
#endif

	if ( level == 1 
	     || amrex::ProperlyNested(crse_ratio, parent->blockingFactor(level),
				       nghost, mf_BA.ixType(), mapper) )
	{
	    StateData& statedata = clev.state[idx];
	    
	    Vector<MultiFab*> smf;
	    Vector<Real> stime;
	    statedata.getData(smf,stime,time);

	    StateDataPhysBCFunct physbcf(statedata,SComp,cgeom);

            crseMF.setDomainBndry(std::numeric_limits<Real>::quiet_NaN(), cgeom);
	    amrex::FillPatchSingleLevel(crseMF,time,smf,stime,SComp,0,NComp,cgeom,physbcf,SComp);
	}
	else
	{
	    FillPatch(clev,crseMF,0,time,idx,SComp,NComp,0);
	}

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
            const Box& dbx = amrex::grow(mfi.validbox(),nghost) & domain_g;
	    
            Vector<BCRec> bcr(ncomp);
	    
            amrex::setBC(dbx,pdomain,SComp,0,NComp,desc.getBCs(),bcr);

            FArrayBox const* crsefab = crseMF.fabPtr(mfi);
            FArrayBox* finefab = mf.fabPtr(mfi);
	    mapper->interp(*crsefab,
			   0,
			   *finefab,
			   DComp,
			   NComp,
			   dbx,
			   crse_ratio,
			   cgeom,
			   geom,
			   bcr,
			   SComp,
			   idx);
	}

	StateDataPhysBCFunct physbcf(state[idx],SComp,geom);
	physbcf.FillBoundary(mf, DComp, NComp, time, SComp);

        DComp += NComp;
    }
}

std::unique_ptr<MultiFab>
AmrLevel::derive (const std::string& name, Real time, int ngrow)
{
    BL_ASSERT(ngrow >= 0);

    std::unique_ptr<MultiFab> mf;

    int index, scomp, ncomp;

    if (isStateVariable(name, index, scomp))
    {
        mf.reset(new MultiFab(state[index].boxArray(), dmap, 1, ngrow, MFInfo(), *m_factory));
        FillPatch(*this,*mf,ngrow,time,index,scomp,1,0);
    }
    else if (const DeriveRec* rec = derive_lst.get(name))
    {
        rec->getRange(0, index, scomp, ncomp);

        const BoxArray& srcBA = state[index].boxArray();

        BoxArray dstBA(srcBA);
        dstBA.convert(rec->deriveType());

	int ngrow_src = ngrow;
	{
	    Box bx0 = srcBA[0];
	    Box bx1 = rec->boxMap()(bx0);
	    int g = bx0.smallEnd(0) - bx1.smallEnd(0);
	    ngrow_src += g;
	}

        MultiFab srcMF(srcBA, dmap, rec->numState(), ngrow_src, MFInfo(), *m_factory);

        for (int k = 0, dc = 0; k < rec->numRange(); k++, dc += ncomp)
        {
            rec->getRange(k, index, scomp, ncomp);
            FillPatch(*this,srcMF,ngrow_src,time,index,scomp,ncomp,dc);
        }

        const int ncomp = rec->numDerive();
        mf.reset(new MultiFab(dstBA, dmap, ncomp, ngrow, MFInfo(), *m_factory));

        if (rec->derFuncFab() != nullptr)
        {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(ngrow);
                FArrayBox* derfab = mf->fabPtr(mfi);
                FArrayBox const* datafab = srcMF.fabPtr(mfi);
                rec->derFuncFab()(bx, *derfab, 0, ncomp, *datafab, geom, time, rec->getBC(), level);
            }
        }
        else
        {
#if defined(AMREX_CRSEGRNDOMP) || (!defined(AMREX_XSDK) && defined(CRSEGRNDOMP))
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
        {
            int         grid_no = mfi.index();
            Real*       ddat    = (*mf)[mfi].dataPtr();
            const int*  dlo     = (*mf)[mfi].loVect();
            const int*  dhi     = (*mf)[mfi].hiVect();
	    const Box&  gtbx    = mfi.growntilebox();
	    const int*  lo      = gtbx.loVect();
	    const int*  hi      = gtbx.hiVect();
            int         n_der   = rec->numDerive();
            Real*       cdat    = srcMF[mfi].dataPtr();
            const int*  clo     = srcMF[mfi].loVect();
            const int*  chi     = srcMF[mfi].hiVect();
            int         n_state = rec->numState();
            const int*  dom_lo  = state[index].getDomain().loVect();
            const int*  dom_hi  = state[index].getDomain().hiVect();
            const Real* dx      = geom.CellSize();
            const int*  bcr     = rec->getBC();
            const RealBox temp    (gtbx,geom.CellSize(),geom.ProbLo());
            const Real* xlo     = temp.lo();
            Real        dt      = parent->dtLevel(level);

	    if (rec->derFunc() != static_cast<DeriveFunc>(0)){
		rec->derFunc()(ddat,AMREX_ARLIM(dlo),AMREX_ARLIM(dhi),&n_der,
			       cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),&n_state,
			       lo,hi,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
			       &level,&grid_no);
	    } else if (rec->derFunc3D() != static_cast<DeriveFunc3D>(0)){
		rec->derFunc3D()(ddat,AMREX_ARLIM_3D(dlo),AMREX_ARLIM_3D(dhi),&n_der,
				 cdat,AMREX_ARLIM_3D(clo),AMREX_ARLIM_3D(chi),&n_state,
				 AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),
				 AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
				 AMREX_ZFILL(dx),AMREX_ZFILL(xlo),
				 &time,&dt,
				 AMREX_BCREC_3D(bcr),
				 &level,&grid_no);
	    } else {
		amrex::Error("AmrLevel::derive: no function available");
	    }
        }
#else
        for (MFIter mfi(srcMF); mfi.isValid(); ++mfi)
        {
            int         grid_no = mfi.index();
            const RealBox gridloc(grids[grid_no],geom.CellSize(),geom.ProbLo());
            Real*       ddat    = (*mf)[mfi].dataPtr();
            const int*  dlo     = (*mf)[mfi].loVect();
            const int*  dhi     = (*mf)[mfi].hiVect();
            int         n_der   = rec->numDerive();
            Real*       cdat    = srcMF[mfi].dataPtr();
            const int*  clo     = srcMF[mfi].loVect();
            const int*  chi     = srcMF[mfi].hiVect();
            int         n_state = rec->numState();
            const int*  dom_lo  = state[index].getDomain().loVect();
            const int*  dom_hi  = state[index].getDomain().hiVect();
            const Real* dx      = geom.CellSize();
            const int*  bcr     = rec->getBC();
            const Real* xlo     = gridloc.lo();
            Real        dt      = parent->dtLevel(level);

	    if (rec->derFunc() != static_cast<DeriveFunc>(0)){
		rec->derFunc()(ddat,AMREX_ARLIM(dlo),AMREX_ARLIM(dhi),&n_der,
			       cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),&n_state,
			       dlo,dhi,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
			       &level,&grid_no);
	    } else if (rec->derFunc3D() != static_cast<DeriveFunc3D>(0)){
		rec->derFunc3D()(ddat,AMREX_ARLIM_3D(dlo),AMREX_ARLIM_3D(dhi),&n_der,
				 cdat,AMREX_ARLIM_3D(clo),AMREX_ARLIM_3D(chi),&n_state,
				 AMREX_ARLIM_3D(dlo),AMREX_ARLIM_3D(dhi),
				 AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
				 AMREX_ZFILL(dx),AMREX_ZFILL(xlo),
				 &time,&dt,
				 AMREX_BCREC_3D(bcr),
				 &level,&grid_no);
	    } else {
		amrex::Error("AmrLevel::derive: no function available");
	    }
        }
#endif
        }
    }
    else
    {
        //
        // If we got here, cannot derive given name.
        //
        std::string msg("AmrLevel::derive(MultiFab*): unknown variable: ");
        msg += name;
        amrex::Error(msg.c_str());
    }

    return mf;
}

void
AmrLevel::derive (const std::string& name, Real time, MultiFab& mf, int dcomp)
{
    BL_ASSERT(dcomp < mf.nComp());

    const int ngrow = mf.nGrow();

    int index, scomp, ncomp;

    if (isStateVariable(name,index,scomp))
    {
        FillPatch(*this,mf,ngrow,time,index,scomp,1,dcomp);
    }
    else if (const DeriveRec* rec = derive_lst.get(name))
    {
        rec->getRange(0,index,scomp,ncomp);

        const BoxArray& srcBA = state[index].boxArray();

	int ngrow_src = ngrow;
	{
	    Box bx0 = srcBA[0];
	    Box bx1 = rec->boxMap()(bx0);
	    int g = bx0.smallEnd(0) - bx1.smallEnd(0);
	    ngrow_src += g;
	}

        MultiFab srcMF(srcBA,dmap,rec->numState(),ngrow_src, MFInfo(), *m_factory);

        for (int k = 0, dc = 0; k < rec->numRange(); k++, dc += ncomp)
        {
            rec->getRange(k,index,scomp,ncomp);

            FillPatch(*this,srcMF,ngrow_src,time,index,scomp,ncomp,dc);
        }

        if (rec->derFuncFab() != nullptr)
        {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox();
                FArrayBox* derfab = mf.fabPtr(mfi);
                FArrayBox const* datafab = srcMF.fabPtr(mfi);
                rec->derFuncFab()(bx, *derfab, dcomp, ncomp, *datafab, geom, time, rec->getBC(), level);
            }
        }
        else
        {
#if defined(AMREX_CRSEGRNDOMP) || (!defined(AMREX_XSDK) && defined(CRSEGRNDOMP))
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
        {
            int         idx     = mfi.index();
            Real*       ddat    = mf[mfi].dataPtr(dcomp);
            const int*  dlo     = mf[mfi].loVect();
            const int*  dhi     = mf[mfi].hiVect();
	    const Box&  gtbx    = mfi.growntilebox();
	    const int*  lo      = gtbx.loVect();
	    const int*  hi      = gtbx.hiVect();
            int         n_der   = rec->numDerive();
            Real*       cdat    = srcMF[mfi].dataPtr();
            const int*  clo     = srcMF[mfi].loVect();
            const int*  chi     = srcMF[mfi].hiVect();
            int         n_state = rec->numState();
            const int*  dom_lo  = state[index].getDomain().loVect();
            const int*  dom_hi  = state[index].getDomain().hiVect();
            const Real* dx      = geom.CellSize();
            const int*  bcr     = rec->getBC();
            const RealBox& temp = RealBox(gtbx,geom.CellSize(),geom.ProbLo());
            const Real* xlo     = temp.lo();
            Real        dt      = parent->dtLevel(level);

	    if (rec->derFunc() != static_cast<DeriveFunc>(0)){
		rec->derFunc()(ddat,AMREX_ARLIM(dlo),AMREX_ARLIM(dhi),&n_der,
			       cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),&n_state,
			       lo,hi,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
			       &level,&idx);
	    } else if (rec->derFunc3D() != static_cast<DeriveFunc3D>(0)){
		rec->derFunc3D()(ddat,AMREX_ARLIM_3D(dlo),AMREX_ARLIM_3D(dhi),&n_der,
				 cdat,AMREX_ARLIM_3D(clo),AMREX_ARLIM_3D(chi),&n_state,
				 AMREX_ARLIM_3D(lo),AMREX_ARLIM_3D(hi),
				 AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
				 AMREX_ZFILL(dx),AMREX_ZFILL(xlo),
				 &time,&dt,
				 AMREX_BCREC_3D(bcr),
				 &level,&idx);
	    } else {
		amrex::Error("AmrLevel::derive: no function available");
	    }
        }
#else
        for (MFIter mfi(srcMF); mfi.isValid(); ++mfi)
        {
            int         idx     = mfi.index();
            Real*       ddat    = mf[mfi].dataPtr(dcomp);
            const int*  dlo     = mf[mfi].loVect();
            const int*  dhi     = mf[mfi].hiVect();
            int         n_der   = rec->numDerive();
            Real*       cdat    = srcMF[mfi].dataPtr();
            const int*  clo     = srcMF[mfi].loVect();
            const int*  chi     = srcMF[mfi].hiVect();
            int         n_state = rec->numState();
            const int*  dom_lo  = state[index].getDomain().loVect();
            const int*  dom_hi  = state[index].getDomain().hiVect();
            const Real* dx      = geom.CellSize();
            const int*  bcr     = rec->getBC();
            const RealBox& temp = RealBox(mf[mfi].box(),geom.CellSize(),geom.ProbLo());
            const Real* xlo     = temp.lo();
            Real        dt      = parent->dtLevel(level);

	    if (rec->derFunc() != static_cast<DeriveFunc>(0)){
		rec->derFunc()(ddat,AMREX_ARLIM(dlo),AMREX_ARLIM(dhi),&n_der,
			       cdat,AMREX_ARLIM(clo),AMREX_ARLIM(chi),&n_state,
			       dlo,dhi,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
			       &level,&idx);
	    } else if (rec->derFunc3D() != static_cast<DeriveFunc3D>(0)){
		rec->derFunc3D()(ddat,AMREX_ARLIM_3D(dlo),AMREX_ARLIM_3D(dhi),&n_der,
				 cdat,AMREX_ARLIM_3D(clo),AMREX_ARLIM_3D(chi),&n_state,
				 AMREX_ARLIM_3D(dlo),AMREX_ARLIM_3D(dhi),
				 AMREX_ARLIM_3D(dom_lo),AMREX_ARLIM_3D(dom_hi),
				 AMREX_ZFILL(dx),AMREX_ZFILL(xlo),
				 &time,&dt,
				 AMREX_BCREC_3D(bcr),
				 &level,&idx);
	    } else {
		amrex::Error("AmrLevel::derive: no function available");
	    }
        }
#endif
        }
    }
    else
    {
        //
        // If we got here, cannot derive given name.
        //
        std::string msg("AmrLevel::derive(MultiFab*): unknown variable: ");
        msg += name;
        amrex::Error(msg.c_str());
    }
}

//! Update the distribution maps in StateData based on the size of the map
void
AmrLevel::UpdateDistributionMaps ( DistributionMapping& update_dmap )
{
    long mapsize = update_dmap.size();

    if (dmap.size() == mapsize)
    { dmap = update_dmap; }

    for (int i = 0; i < state.size(); ++i)
    {
       if (state[i].DistributionMap().size() == mapsize)
          { state[i].setDistributionMap(update_dmap); }
    }
}



Vector<int>
AmrLevel::getBCArray (int State_Type,
                      int gridno,
                      int strt_comp,
                      int ncomp)
{
    Vector<int> bc(2*AMREX_SPACEDIM*ncomp);

    BCRec bcr;

    for (int n = 0; n < ncomp; n++)
    {
        bcr = state[State_Type].getBC(strt_comp+n,gridno);
        const int* b_rec = bcr.vect();
        for (int m = 0; m < 2*AMREX_SPACEDIM; m++)
            bc[2*AMREX_SPACEDIM*n + m] = b_rec[m];
    }

    return bc;
}

int
AmrLevel::okToRegrid ()
{
    return true;
}

void
AmrLevel::setPlotVariables ()
{
    ParmParse pp("amr");

    if (pp.contains("plot_vars"))
    {
        std::string nm;
      
        int nPltVars = pp.countval("plot_vars");
      
        for (int i = 0; i < nPltVars; i++)
        {
            pp.get("plot_vars", nm, i);

            if (nm == "ALL") 
                parent->fillStatePlotVarList();
            else if (nm == "NONE")
                parent->clearStatePlotVarList();
            else
                parent->addStatePlotVar(nm);
        }
    }
    else 
    {
        //
        // The default is to add them all.
        //
        parent->fillStatePlotVarList();
    }
  
    if (pp.contains("derive_plot_vars"))
    {
        std::string nm;
      
        int nDrvPltVars = pp.countval("derive_plot_vars");
      
        for (int i = 0; i < nDrvPltVars; i++)
        {
            pp.get("derive_plot_vars", nm, i);

            if (nm == "ALL") 
                parent->fillDerivePlotVarList();
            else if (nm == "NONE")
                parent->clearDerivePlotVarList();
            else
                parent->addDerivePlotVar(nm);
        }
    }
    else 
    {
        //
        // The default is to add none of them.
        //
        parent->clearDerivePlotVarList();
    }
}

void
AmrLevel::setSmallPlotVariables ()
{
    ParmParse pp("amr");

    if (pp.contains("small_plot_vars"))
    {
        std::string nm;
      
        int nPltVars = pp.countval("small_plot_vars");
      
        for (int i = 0; i < nPltVars; i++)
        {
            pp.get("small_plot_vars", nm, i);

            if (nm == "ALL")
                parent->fillStateSmallPlotVarList();
            else if (nm == "NONE")
                parent->clearStateSmallPlotVarList();
            else
                parent->addStateSmallPlotVar(nm);
        }
    }
    else 
    {
        //
        // The default is to use none.
        //
        parent->clearStateSmallPlotVarList();
    }

    if (pp.contains("derive_small_plot_vars"))
    {
        std::string nm;
      
        int nDrvPltVars = pp.countval("derive_small_plot_vars");
      
        for (int i = 0; i < nDrvPltVars; i++)
        {
            pp.get("derive_small_plot_vars", nm, i);

            if (nm == "ALL") 
                parent->fillDeriveSmallPlotVarList();
            else if (nm == "NONE")
                parent->clearDeriveSmallPlotVarList();
            else
                parent->addDeriveSmallPlotVar(nm);
        }
    }
    else 
    {
        //
        // The default is to add none of them.
        //
        parent->clearDeriveSmallPlotVarList();
    }
  
}

AmrLevel::TimeLevel
AmrLevel::which_time (int  indx, Real time) const noexcept
{
    const Real oldtime = state[indx].prevTime();
    const Real newtime = state[indx].curTime();
    const Real haftime = .5 * (oldtime + newtime);
    const Real qtime = oldtime + 0.25*(newtime-oldtime);
    const Real tqtime = oldtime + 0.75*(newtime-oldtime);
    const Real epsilon = 0.001 * (newtime - oldtime);

    BL_ASSERT(time >= oldtime-epsilon && time <= newtime+epsilon);
    
    if (time >= oldtime-epsilon && time <= oldtime+epsilon)
    {
        return AmrOldTime;
    }
    else if (time >= newtime-epsilon && time <= newtime+epsilon)
    {
        return AmrNewTime;
    }
    else if (time >= haftime-epsilon && time <= haftime+epsilon)
    {
        return AmrHalfTime;
    }
    else if (time >= qtime-epsilon && time <= qtime+epsilon)
    {
        return Amr1QtrTime;
    }
    else if (time >= tqtime-epsilon && time <= tqtime+epsilon)
    {
        return Amr3QtrTime;
    }
    return AmrOtherTime;
}

Real
AmrLevel::estimateWork ()
{
    return 1.0*countCells();
}

bool
AmrLevel::writePlotNow ()
{
    return false;
}

bool
AmrLevel::writeSmallPlotNow ()
{
    return false;
}

const BoxArray& AmrLevel::getAreaNotToTag () noexcept
{
    return m_AreaNotToTag;
}

const Box& AmrLevel::getAreaToTag () noexcept
{
    return m_AreaToTag;
}

void AmrLevel::setAreaNotToTag (BoxArray& ba) noexcept
{
    m_AreaNotToTag = ba;
}

void AmrLevel::constructAreaNotToTag ()
{
    if (level == 0 || !parent->useFixedCoarseGrids() || parent->useFixedUpToLevel()>level)
        return;

    // We are restricting the tagging on the finest fixed level
    if (parent->useFixedUpToLevel()==level)
    {
        // We use the next coarser level shrunk by one blockingfactor
        //    as the region in which we allow tagging. 
        // Why level-1? Because we always use the full domain at level 0 
        //    and therefore level 0 in initialba is level 1 in the AMR hierarchy, etc.
        const Vector<BoxArray>& initialba = parent->getInitialBA();
        Box tagarea(initialba[level-1].minimalBox());
        tagarea.grow(-parent->blockingFactor(level));
        m_AreaToTag = tagarea;

        // We disallow tagging in the remaining part of the domain.
        BoxArray tagba = amrex::boxComplement(parent->Geom(level).Domain(),m_AreaToTag);
        m_AreaNotToTag = tagba;

        BoxArray bxa(parent->Geom(level).Domain());
        BL_ASSERT(bxa.contains(m_AreaNotToTag));
    }

    if (parent->useFixedUpToLevel()<level)
    {
        Box tagarea = parent->getLevel(level-1).getAreaToTag();
        tagarea.refine(parent->refRatio(level-1));
        tagarea.grow(-parent->blockingFactor(level));
        m_AreaToTag = tagarea;
        BoxArray tagba = amrex::boxComplement(parent->Geom(level).Domain(),m_AreaToTag);
        m_AreaNotToTag = tagba;
    }
}

void
AmrLevel::FillPatch (AmrLevel& amrlevel,
		     MultiFab& leveldata,
                     int       boxGrow,
                     Real      time,
                     int       index,
                     int       scomp,
                     int       ncomp,
                     int       dcomp)
{
    BL_ASSERT(dcomp+ncomp-1 <= leveldata.nComp());
    BL_ASSERT(boxGrow <= leveldata.nGrow());
    FillPatchIterator fpi(amrlevel, leveldata, boxGrow, time, index, scomp, ncomp);
    const MultiFab& mf_fillpatched = fpi.get_mf();
    MultiFab::Copy(leveldata, mf_fillpatched, 0, dcomp, ncomp, boxGrow);
}

void
AmrLevel::FillPatchAdd (AmrLevel& amrlevel,
                        MultiFab& leveldata,
                        int       boxGrow,
                        Real      time,
                        int       index,
                        int       scomp,
                        int       ncomp,
                        int       dcomp)
{
    BL_ASSERT(dcomp+ncomp-1 <= leveldata.nComp());
    BL_ASSERT(boxGrow <= leveldata.nGrow());
    FillPatchIterator fpi(amrlevel, leveldata, boxGrow, time, index, scomp, ncomp);
    const MultiFab& mf_fillpatched = fpi.get_mf();
    MultiFab::Add(leveldata, mf_fillpatched, 0, dcomp, ncomp, boxGrow);
}

void
AmrLevel::LevelDirectoryNames (const std::string &dir,
                               std::string &LevelDir,
                               std::string &FullPath)
{
    LevelDir = amrex::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    FullPath = dir;
    if( ! FullPath.empty() && FullPath.back() != '/') {
        FullPath += '/';
    }
    FullPath += LevelDir;
}

void
AmrLevel::CreateLevelDirectory (const std::string &dir)
{
    // Build directory to hold the MultiFabs in the StateData at this level.
    // The directory is relative the the directory containing the Header file.

    std::string LevelDir, FullPath;
    LevelDirectoryNames(dir, LevelDir, FullPath);

    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
        amrex::CreateDirectoryFailed(FullPath);
      }
    }
    levelDirectoryCreated = true;
}


#ifdef USE_PERILLA
    void FillPatchIterator::FillPatchSingleLevelPush (Amr& amr, MultiFab& mf, Real time,
                                   Vector<MultiFab*>& smf, Vector<Real>& stime,
                                   RegionGraph* destGraph, RegionGraph* srcGraph, int f,
                                   MultiFab *dmf,
                                   int scomp, int dcomp, int ncomp,
                                   const Geometry& geom, StateDataPhysBCFunct& physbcf, bool singleT)
      {

        BL_PROFILE("FillPatchSingleLevel");

        BL_ASSERT(scomp+ncomp <= smf[0]->nComp());
        BL_ASSERT(dcomp+ncomp <= mf.nComp());
        BL_ASSERT(smf.size() == stime.size());
        BL_ASSERT(smf.size() != 0);

        int tg = WorkerThread::perilla_wid();
        int nt = WorkerThread::perilla_wtid();

        if (smf.size() == 1)
          {
          //mf.copy(smf[0], scomp, dcomp, ncomp, 0, mf.nGrow(), geom.periodicity());        
            Perilla::multifabCopyPushAsync( destGraph, srcGraph, &mf, smf[0],  f, dcomp, scomp, ncomp, mf.nGrow(), 0, singleT);
          }
        else if (smf.size() == 2)
          {
            BL_ASSERT(smf[0]->boxArray() == smf[1]->boxArray());
            //PArray<MultiFab> raii(PArrayManage);
            //MultiFab * dmf;
            int destcomp;
            bool sameba;
            if (mf.boxArray() == smf[0]->boxArray())
              {
                //dmf = &mf;
                destcomp = dcomp;
                sameba = true;

                int fis = smf[0]->IndexArray()[f];
                int fid = mf.IndexArray()[f];
                const Box& bx = mf[fid].box();
                mf[fid].linInterp((*smf[0])[fis],
                                scomp,
                                (*smf[1])[fis],
                                scomp,
                                stime[0],
                                stime[1],
                                time,
                                bx,
                                destcomp,
                                ncomp);
                Perilla::fillBoundaryPush(destGraph, &mf, f);
              }
            else
              {

                //dmf = raii.push_back(new MultiFab(smf[0].boxArray(), ncomp, 0));
                //MultiFab dmf(smf[0].boxArray(), ncomp, 0);
                destcomp = 0;
                sameba = false;

                assert(smf[0]);
                assert(smf[0]->IndexArray().size()>f);
                assert(dmf);
                assert(dmf->IndexArray().size()>f);
                int fis = smf[0]->IndexArray()[f];
                int fid = dmf->IndexArray()[f];

                for(int t=0; t<srcGraph->fabTiles[f]->numTiles; t++)
                  if( singleT || t % (perilla::NUM_THREADS_PER_TEAM-1) == nt)
                    {
                      const Box& bx = *(srcGraph->fabTiles[f]->tileBx[t]);

                      //const Box& bx = (*dmf)[fid].box();
                      (*dmf)[fid].linInterp((*smf[0])[fis],
                                            scomp,
                                            (*smf[1])[fis],
                                            scomp,
                                            stime[0],
                                            stime[1],
                                            time,
                                            bx,
                                            destcomp,
                                            ncomp);
                    }
                if(!singleT)
                  srcGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1); // Barrier to synchronize team threads
                int src_ngrow = 0;
                int dst_ngrow = mf.nGrow();
                Perilla::multifabCopyPushAsync( destGraph, srcGraph, &mf, dmf, f, dcomp, 0, ncomp, mf.nGrow(), 0, singleT);
              }
          }
        else
          {
            amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
          }
      }


    void FillPatchIterator::FillPatchTwoLevelsPush (Amr& amr, MultiFab& mf, Real time,
                                 Vector<MultiFab*>& cmf, Vector<Real>& ct,
                                 Vector<MultiFab*>& fmf, Vector<Real>& ft,
                                 RegionGraph* destGraph, RegionGraph* csrcGraph, RegionGraph* fsrcGraph, int f,
                                 FillPatchIterator* fpIter,
                                 MultiFab *dmf,
                                 MultiFab *dmff,
                                 int scomp, int dcomp, int ncomp,
                                 const Geometry& cgeom, const Geometry& fgeom,
                                 StateDataPhysBCFunct& cbc, StateDataPhysBCFunct& fbc,
                                 const IntVect& ratio,
                                 Interpolater* mapper, const Vector<BCRec>& bcs, unsigned char pushLevel,  bool singleT)
    {
        BL_PROFILE("FillPatchTwoLevels");

        int ngrow = mf.nGrow();

        if(f>=0){//fill only this fab
            if(pushLevel & 0x01 )
            {
                if (ngrow > 0 || mf.getBDKey() != fmf[0]->getBDKey())
                {

                    if (!fpIter->m_fpc->ba_crse_patch.empty())
                    {
                        FillPatchSingleLevelPush(amr, *(fpIter->m_mf_crse_patch), time, cmf, ct, fpIter->m_rg_crse_patch, csrcGraph, f, dmf, scomp, 0, ncomp, cgeom, cbc, singleT);
                    }
                }
            }
            if((pushLevel & 0x02) && (pushLevel != 0x03))
            {
                 FillPatchSingleLevelPush(amr, mf, time, fmf, ft, destGraph, fsrcGraph, f, dmff, scomp, dcomp, ncomp, fgeom, fbc, singleT);
            }
        }else{ //fill the whole multifab
            if(pushLevel & 0x01 && pushLevel & 0x02)
            {
                int tg = perilla::wid();
                for(int fi=0; fi < fmf[0]->IndexArray().size(); fi++)
                {
                    if(WorkerThread::isMyRegion(tg,fi))
                    {
                      FillPatchSingleLevelPush(amr, mf, time, fmf, ft, destGraph, fsrcGraph, fi, dmff, scomp, dcomp, ncomp, fgeom, fbc, singleT);
                    } 
                }
            }       
            if(pushLevel & 0x04)
            {
                int tg = perilla::wid();
                for(int fi=0; fi < fmf[0]->IndexArray().size(); fi++)
                {
                    if(WorkerThread::isMyRegion(tg,fi))
                    {
                        FillPatchSingleLevelPush(amr, mf, time, fmf, ft, destGraph, fsrcGraph, fi, dmff, scomp, dcomp, ncomp, fgeom, fbc, singleT);
                    }   
                }
            }       
        }       
    }       

    void FillPatchIterator::FillPatchSingleLevelPull (MultiFab& mf, Real time,
                                   Vector<MultiFab*>& smf, Vector<Real>& stime,
                                   RegionGraph* destGraph, RegionGraph* srcGraph, int f,
                                   int scomp, int dcomp, int ncomp,
                                   const Geometry& geom, StateDataPhysBCFunct& physbcf, bool singleT)
      {

        BL_PROFILE("FillPatchSingleLevel");

        BL_ASSERT(scomp+ncomp <= smf[0]->nComp());
        BL_ASSERT(dcomp+ncomp <= mf.nComp());
        BL_ASSERT(smf.size() == stime.size());
        BL_ASSERT(smf.size() != 0);

        if (smf.size() == 1)
        {
          //mf.copy(smf[0], scomp, dcomp, ncomp, 0, mf.nGrow(), geom.periodicity());      
          Perilla::multifabCopyPull( destGraph, srcGraph, &mf, smf[0], f, dcomp, scomp, ncomp, mf.nGrow(), 0, singleT);
        }
        else if (smf.size() == 2)
        {
            BL_ASSERT(smf[0]->boxArray() == smf[1]->boxArray());
            //Vector<MultiFab> raii(PArrayManage);
            MultiFab * dmf;
            int destcomp;
            bool sameba;
            if (mf.boxArray() == smf[0]->boxArray()) {
              dmf = &mf;
              destcomp = dcomp;
              sameba = true;
            } else {
              //dmf = srcGraph->assocMF;              
              destcomp = 0;
              sameba = false;
            }
            if (sameba)
            {
                // Note that when sameba is true mf's BoxArray is nonoverlapping.
                // So FillBoundary is safe.
                //mf.FillBoundary(dcomp,ncomp,geom.periodicity());
              Perilla::fillBoundaryPull(destGraph, dmf, f, singleT);
            }
            else
            {
                int src_ngrow = 0;
                int dst_ngrow = mf.nGrow();
                MultiFab* dummyMF;
                //mf.copy(*dmf, 0, dcomp, ncomp, src_ngrow, dst_ngrow, geom.periodicity());
                Perilla::multifabCopyPull( destGraph, srcGraph, &mf, dummyMF, f, dcomp, 0, ncomp, mf.nGrow(), 0, singleT);
            }
        }
        else {
            amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
        }
#if 0
        physbcf.doit_fab(mf, f, dcomp, ncomp, time);
#endif
    }

    void FillPatchIterator::FillPatchTwoLevelsPull (MultiFab& mf, Real time,
                                 Vector<MultiFab*>& cmf, Vector<Real>& ct,
                                 Vector<MultiFab*>& fmf, Vector<Real>& ft,
                                 RegionGraph* destGraph, RegionGraph* csrcGraph, RegionGraph* fsrcGraph, int f,
                                 FillPatchIterator* fpIter,
                                 int scomp, int dcomp, int ncomp,
                                 const Geometry& cgeom, const Geometry& fgeom,
                                 StateDataPhysBCFunct& cbc, StateDataPhysBCFunct& fbc,
                                 const IntVect& ratio,
                                 Interpolater* mapper, const Vector<BCRec>& bcs, bool singleT)
    {
        BL_PROFILE("FillPatchTwoLevels");

        int ngrow = mf.nGrow();

        int tg = WorkerThread::perilla_wid();
        int nt = WorkerThread::perilla_wtid();

        if (ngrow > 0 || mf.getBDKey() != fmf[0]->getBDKey())
        {

            if ( ! fpIter->m_fpc->ba_crse_patch.empty())
            {

                int idummy1=0, idummy2=0;
                bool cc = fpIter->m_fpc->ba_crse_patch.ixType().cellCentered();
                  {
                    int gi = mf.IndexArray()[f];
                    for(int i=0; i<destGraph->task[f]->depTaskIDs.size();i++)
                      {
                        int li = destGraph->task[f]->depTaskIDs[i];
                        int mfi = fpIter->m_mf_crse_patch[0].IndexArray()[li];
                        FillPatchSingleLevelPull(*(fpIter->m_mf_crse_patch), time, cmf, ct, fpIter->m_rg_crse_patch, csrcGraph, li, scomp, 0, ncomp, cgeom, cbc, singleT);
                      }
                    if(!singleT)
                      destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);

                    int nt = WorkerThread::perilla_wtid();
                    Box fdomain = fgeom.Domain();
                    for(int i=0; i<destGraph->task[f]->depTaskIDs.size();i++)
                      {
                        int li = destGraph->task[f]->depTaskIDs[i];
                        int mfi = fpIter->m_mf_crse_patch[0].IndexArray()[li];
                        if(singleT)
                          {
                            const Box& dbx = fpIter->m_fpc->dst_boxes[li];
                            //Array<BCRec> bcr(ncomp);
		            Vector<BCRec> bcr(ncomp);
                            amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);

                            mapper->interp(fpIter->m_mf_crse_patch[0][mfi],
                                           0,
                                           mf[gi],
                                           dcomp,
                                           ncomp,
                                           dbx,
                                           ratio,
                                           cgeom,
                                           fgeom,
                                           bcr,
                                           idummy1, idummy2);
                          }
                        else
                          {
                          if(!cc)
                            {
                              if(WorkerThread::perilla_isMasterWorkerThread())
                                {
                                  const Box& dbx = fpIter->m_fpc->dst_boxes[li];
                                  //Box fdomain = fgeom.Domain();

                                  Vector<BCRec> bcr(ncomp);
                                  amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);

                                mapper->interp(fpIter->m_mf_crse_patch[0][mfi],
                                               0,
                                               mf[gi],
                                               dcomp,
                                               ncomp,
                                               dbx,
                                               ratio,
                                               cgeom,
                                               fgeom,
                                               bcr,
                                               idummy1, idummy2);

                                }
                            }
                          else
                            {
                              if(i % (perilla::NUM_THREADS_PER_TEAM-1) == nt-1)
                                {

                                  const Box& dbx = fpIter->m_fpc->dst_boxes[li];

                                  Vector<BCRec> bcr(ncomp);
                                  amrex::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);

                                  mapper->interp(fpIter->m_mf_crse_patch[0][mfi],
                                                 0,
                                                 mf[gi],
                                                 dcomp,
                                                 ncomp,
                                                 dbx,
                                                 ratio,
                                                 cgeom,
                                                 fgeom,
                                                 bcr,
                                                 idummy1, idummy2);

                                }
                            }
                          }
                        //destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
                      }
                    if(!singleT)
                      destGraph->worker[tg]->barr->sync(perilla::NUM_THREADS_PER_TEAM-1);
                  }
            }
        }


        FillPatchSingleLevelPull(mf, time, fmf, ft, destGraph, fsrcGraph, f, scomp, dcomp, ncomp, fgeom, fbc, singleT);
    }

void FillPatchIterator::FillFromTwoLevelsPush (Real time,
                                          int index,
                                          int scomp,
                                          int dcomp,
                                          int ncomp,
                                          int f,
                                          unsigned char pushLevel,
                                          bool singleT)
{
    int ilev_fine = m_amrlevel.level;
    int ilev_crse = ilev_fine-1;
    
    BL_ASSERT(ilev_crse >= 0);
    
    AmrLevel& fine_level = m_amrlevel;
    AmrLevel& crse_level = m_amrlevel.parent->getLevel(ilev_crse);
        
    Geometry* tgeom_fine = &fine_level.geom;
    Geometry* tgeom_crse = &crse_level.geom;
    
    Vector<MultiFab*> tsmf_crse;
    Vector<MultiFab*> tsmf_fine;
    Vector<Real> tstime_crse;
    Vector<Real> tstime_fine;
    StateData& statedata_crse = crse_level.state[index];
    statedata_crse.getData(tsmf_crse,tstime_crse,time);
    StateDataPhysBCFunct* tphysbcf_crse = new StateDataPhysBCFunct(statedata_crse,scomp,*geom_crse);
        
    StateData& statedata_fine = fine_level.state[index];
    statedata_fine.getData(tsmf_fine,tstime_fine,time);
    StateDataPhysBCFunct* tphysbcf_fine = new StateDataPhysBCFunct(statedata_fine,scomp,*geom_fine);
        
    const StateDescriptor& desc = AmrLevel::desc_lst[index];
    
    FillPatchTwoLevelsPush(*(m_amrlevel.parent), m_fabs, time,
                                   tsmf_crse, tstime_crse,
                                   tsmf_fine, tstime_fine,
                                   destGraph, csrcGraph, fsrcGraph, f,
                                   this,
                                   dmf,
                                   dmff,
                                   scomp, dcomp, ncomp, 
                                   *tgeom_crse, *tgeom_fine,
                                   *tphysbcf_crse, *tphysbcf_fine,
                                   crse_level.fineRatio(), 
                                   desc.interp(scomp), desc.getBCs(), pushLevel, singleT);
}

void
FillPatchIterator::FillFromTwoLevelsPushOnly (Real time, int index, int scomp, int dcomp, int ncomp, int f, unsigned char pushLevel, bool singleT)
{
    int ilev_fine = m_amrlevel.level;
    int ilev_crse = ilev_fine-1;

    BL_ASSERT(ilev_crse >= 0);

    AmrLevel& fine_level = m_amrlevel;
    AmrLevel& crse_level = m_amrlevel.parent->getLevel(ilev_crse);

    //if(physbcf_fine == NULL && physbcf_crse == NULL)
    //{

    Geometry* tgeom_fine = &fine_level.geom;
    Geometry* tgeom_crse = &crse_level.geom;

    Vector<MultiFab*> tsmf_crse;
    Vector<MultiFab*> tsmf_fine;
    Vector<Real> tstime_crse;
    Vector<Real> tstime_fine;
    StateData& statedata_crse = crse_level.state[index];
    statedata_crse.getData(tsmf_crse,tstime_crse,time);
    StateDataPhysBCFunct* tphysbcf_crse = new StateDataPhysBCFunct(statedata_crse,scomp,*geom_crse);

    StateData& statedata_fine = fine_level.state[index];
    statedata_fine.getData(tsmf_fine,tstime_fine,time);
    StateDataPhysBCFunct* tphysbcf_fine = new StateDataPhysBCFunct(statedata_fine,scomp,*geom_fine);
        //}

    const StateDescriptor& desc = AmrLevel::desc_lst[index];

    FillPatchTwoLevelsPush(*(m_amrlevel.parent), m_fabs, time,
                                   tsmf_crse, tstime_crse,
                                   tsmf_fine, tstime_fine,
                                   destGraph, csrcGraph, fsrcGraph, f,
                                   this,
                                   dmf,
                                   dmff,
                                   scomp, dcomp, ncomp,
                                   *tgeom_crse, *tgeom_fine,
                                   *tphysbcf_crse, *tphysbcf_fine,
                                   crse_level.fineRatio(),
                                   desc.interp(scomp), desc.getBCs(), pushLevel, singleT);
}

void FillPatchIterator::FillFromTwoLevelsPull (Real time, int index, int scomp, int dcomp, int ncomp, int f, bool singleT)
{

    int ilev_fine = m_amrlevel.level;
    int ilev_crse = ilev_fine-1;

    BL_ASSERT(ilev_crse >= 0);

    AmrLevel& fine_level = m_amrlevel;
    AmrLevel& crse_level = m_amrlevel.parent->getLevel(ilev_crse);

    Geometry* tgeom_fine = &fine_level.geom;
    Geometry* tgeom_crse = &crse_level.geom;

    Vector<MultiFab*> tsmf_crse;
    Vector<MultiFab*> tsmf_fine;
    Vector<Real> tstime_crse;
    Vector<Real> tstime_fine;
    StateData& statedata_crse = crse_level.state[index];
    statedata_crse.getData(tsmf_crse,tstime_crse,time);
    StateDataPhysBCFunct* tphysbcf_crse = new StateDataPhysBCFunct(statedata_crse,scomp,*geom_crse);

    StateData& statedata_fine = fine_level.state[index];
    statedata_fine.getData(tsmf_fine,tstime_fine,time);
    StateDataPhysBCFunct* tphysbcf_fine = new StateDataPhysBCFunct(statedata_fine,scomp,*geom_fine);


    const StateDescriptor& desc = AmrLevel::desc_lst[index];

    FillPatchTwoLevelsPull(m_fabs, time,
                                   tsmf_crse, tstime_crse,
                                   tsmf_fine, tstime_fine,
                                   destGraph, csrcGraph, fsrcGraph, f,
                                   this,
                                   scomp, dcomp, ncomp,
                                   *tgeom_crse, *tgeom_fine,
                                   *tphysbcf_crse, *tphysbcf_fine,
                                   crse_level.fineRatio(),
                                   desc.interp(scomp), desc.getBCs(), singleT);
}

void
FillPatchIterator::FillPatchPush (int  boxGrow,
                             Real time,
                             int  index,
                             int  scomp,
                             int  ncomp,
                             int f,
                             unsigned char pushLevel,
                             bool singleT)
{
    BL_PROFILE("FillPatchIterator::InitializePush");

    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

    //const IndexType& boxType = m_leveldata.boxArray().ixType();
    const int level = m_amrlevel.level;

    for (int i = 0, DComp = 0; i < m_range.size(); i++)
    {
      if(i>0)
        amrex::Abort("**** Error in FillPatchIterator::Initialize:  non contigeous components not implemented");

      const int SComp = m_range[i].first;
      const int NComp = m_range[i].second;

      if (level == 0)
        {
          FillPatchSingleLevelPush (*(m_amrlevel.parent), m_fabs, time, smf, stime, destGraph, fsrcGraph, f, dmf, SComp, DComp, NComp, *geom, *physbcf, singleT);
        }
      else
        {
          if (level == 1 || isProperlyNested)
            {
              FillFromTwoLevelsPushOnly(time, index, SComp, DComp, NComp, f, pushLevel, singleT);
            } else {
            amrex::Abort("**** Error in FillPatchIterator::Initialize:  !ProperlyNested not implemented");
          }
        }
      DComp += NComp;
    }
}

void
FillPatchIterator::FillPatchPull (int  boxGrow,
                                   Real time,
                                   int  index,
                                   int  scomp,
                                   int  ncomp,
                                   int f,
                                   bool singleT)
{
    BL_PROFILE("FillPatchIterator::InitializePull");

    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

    //const IndexType& boxType = m_leveldata.boxArray().ixType();
    const int level = m_amrlevel.level;


    for (int i = 0, DComp = 0; i < m_range.size(); i++)
      {
        if(i>0)
          amrex::Abort("**** Error in FillPatchIterator::Initialize:  non contigeous components not implemented");

        const int SComp = m_range[i].first;
        const int NComp = m_range[i].second;

        if (level == 0)
          {
            FillPatchSingleLevelPull (m_fabs, time, smf, stime, destGraph, fsrcGraph, f, SComp, DComp, NComp, *geom, *physbcf, singleT);
          }
        else
          {
            if (level == 1 || isProperlyNested)
              {
                FillFromTwoLevelsPull(time, index, SComp, DComp, NComp, f, singleT);
              } else {
              amrex::Abort("**** Error in FillPatchIterator::Initialize:  !ProperlyNested not implemented");
            }
          }
        //if(WorkerThread::isTeamMasterThread(tid))
          {
            const MultiFab& mf_fillpatched = m_fabs;
            
            if(singleT)
              {
                for(int t=0; t<destGraph->fabTiles_gtbx[f]->numTiles; t++)
                  {
                    const Box& bx = *(destGraph->fabTiles_gtbx[f]->tileBx[t]);
                    MultiFab::Copy(m_leveldata, mf_fillpatched, f, 0, DComp, ncomp, bx);
                  } 
              }   
            else
              {
                perilla::syncAllWorkerThreads();
                int nt = WorkerThread::perilla_wtid();
                for(int t=0; t<destGraph->fabTiles_gtbx[f]->numTiles; t++)
                  if(t % (perilla::NUM_THREADS_PER_TEAM-1) == nt-1)
                    {
                      const Box& bx = *(destGraph->fabTiles_gtbx[f]->tileBx[t]);
                      MultiFab::Copy(m_leveldata, mf_fillpatched, f, 0, DComp, ncomp, bx);
                    } 
                perilla::syncAllWorkerThreads();
              } 
          }   
        DComp += NComp;
    }   
    //
    // Call hack to touch up fillPatched data.
    //
    /*m_amrlevel.set_preferred_boundary_values(m_fabs,
                                             index,
                                             scomp,
                                             0,
                                             ncomp,
                                             time);*/

    }

void
FillPatchIterator::initFillPatch(int boxGrow, int time, int index, int scomp, int ncomp, int iter)
{
#ifdef USE_PERILLA_PTHREADS
//    perilla::syncAllThreads();
#endif
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());


    const StateDescriptor& desc = AmrLevel::desc_lst[index];
#ifdef USE_PERILLA_PTHREADS
//    if(perilla::isMasterThread())
#endif
    {
        m_ncomp = ncomp;
        m_range = desc.sameInterps(scomp,ncomp);

        m_fabs.define(m_leveldata.boxArray(),m_leveldata.DistributionMap(), m_ncomp,boxGrow);

        BL_ASSERT(m_leveldata.DistributionMap() == m_fabs.DistributionMap());
    }
#ifdef USE_PERILLA_PTHREADS
//    perilla::syncAllThreads();
#endif

    const IndexType& boxType = m_leveldata.boxArray().ixType();
    const int level = m_amrlevel.level;

    for (int i = 0, DComp = 0; i < m_range.size(); i++)
    {
        const int SComp = m_range[i].first;
        const int NComp = m_range[i].second;
        int dcomp = DComp;
        if (level == 0)
        {
#if 1
#ifdef USE_PERILLA_PTHREADS
//            perilla::syncAllThreads();
//            if(perilla::isMasterThread())
#endif

            {
                BL_ASSERT(m_amrlevel.level == 0);
                StateData& statedata = m_amrlevel.state[index];
                statedata.getData(smf,stime,time);
                geom = &m_amrlevel.geom;
                physbcf = new StateDataPhysBCFunct(statedata,scomp,*geom);
		stateDataList.push_back(physbcf);
                BL_ASSERT(scomp+ncomp <= smf[0]->nComp());
                BL_ASSERT(dcomp+ncomp <= m_fabs.nComp());
                BL_ASSERT(smf.size() == stime.size());
                BL_ASSERT(smf.size() != 0);
            }
#ifdef USE_PERILLA_PTHREADS
//            perilla::syncAllThreads();
#endif
            if (smf.size() == 1)
            {
#ifdef USE_PERILLA_PTHREADS
//                if(perilla::isMasterThread())
#endif
                {
                    destGraph = new RegionGraph(m_fabs.IndexArray().size());
                    fsrcGraph = new RegionGraph(smf[0]->IndexArray().size());
		    regionList.push_back(destGraph);
		    regionList.push_back(fsrcGraph);
                }
#ifdef USE_PERILLA_PTHREADS
//                perilla::syncAllThreads();
#endif

#if 1
                Perilla::multifabExtractCopyAssoc( destGraph, fsrcGraph, m_fabs, *(smf[0]), (const int) ncomp, m_fabs.nGrow(), 0, geom->periodicity());
#endif
#ifdef USE_PERILLA_PTHREADS
//                perilla::syncAllThreads();
//                if(perilla::isMasterThread())
#endif
                {
                    m_amrlevel.parent->graphArray[level].push_back(destGraph);
                    m_amrlevel.parent->graphArray[level].push_back(fsrcGraph);
                }
#ifdef USE_PERILLA_PTHREADS
//              perilla::syncAllThreads();
#endif
            }
            else if (smf.size() == 2)
            {
                BL_ASSERT(smf[0]->boxArray() == smf[1]->boxArray());

                if (m_fabs.boxArray() == smf[0]->boxArray())
                {
#ifdef USE_PERILLA_PTHREADS
//                    if(perilla::isMasterThread())
#endif
                    {
                        dmf = &m_fabs;
                        destGraph = new RegionGraph(m_fabs.IndexArray().size());
		        regionList.push_back(destGraph);
                    }
#ifdef USE_PERILLA_PTHREADS
 //                   perilla::syncAllThreads();
#endif
                        Perilla::multifabBuildFabCon(destGraph, m_fabs, geom->periodicity());
#ifdef USE_PERILLA_PTHREADS
//                    perilla::syncAllThreads();
//                    if(perilla::isMasterThread())
#endif
                    {
                        m_amrlevel.parent->graphArray[level].push_back(destGraph);
                    }
                }
                else
                {
#ifdef USE_PERILLA_PTHREADS
//                    if(perilla::isMasterThread())
#endif
                    {
                        dmf = new MultiFab(smf[0]->boxArray(), smf[0]->DistributionMap(), ncomp, 0);
                        //dmf->initVal(); // for Perilla NUMA
                        destGraph = new RegionGraph(m_fabs.IndexArray().size());
                        fsrcGraph = new RegionGraph(dmf->IndexArray().size());
                        fsrcGraph->buildTileArray(*dmf);
		        regionList.push_back(destGraph);
		        regionList.push_back(fsrcGraph);
		        mfList.push_back(dmf);
                    }
#ifdef USE_PERILLA_PTHREADS
//                    perilla::syncAllThreads();
#endif

                    Perilla::multifabExtractCopyAssoc(destGraph, fsrcGraph, m_fabs, *dmf, ncomp, m_fabs.nGrow(), 0, geom->periodicity());

#ifdef USE_PERILLA_PTHREADS
//                    perilla::syncAllThreads();
//                    if(perilla::isMasterThread())
#endif
                    {
                        m_amrlevel.parent->graphArray[level].push_back(destGraph);
                        m_amrlevel.parent->graphArray[level].push_back(fsrcGraph);
                    }
#ifdef USE_PERILLA_PTHREADS
//                    perilla::syncAllThreads();
#endif
                }
             }
             else
             {
                //BoxLib::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
             }
#endif
            //-------------------------------------------------- FillFromLevel0 initialization completed
          }
        else
          {
#ifdef USE_PERILLA_PTHREADS
//             perilla::syncAllThreads();
//             if(perilla::isMasterThread())
#endif
             {
                isProperlyNested = amrex::ProperlyNested(m_amrlevel.crse_ratio,
                                                      m_amrlevel.parent->blockingFactor(m_amrlevel.level),
                                                      boxGrow, boxType, desc.interp(SComp));
             }
#ifdef USE_PERILLA_PTHREADS
//           perilla::syncAllThreads();
#endif
             if (level == 1 || isProperlyNested)
                {
                  int ilev_fine = m_amrlevel.level;
                  int ilev_crse = ilev_fine-1;
                  BL_ASSERT(ilev_crse >= 0);
                  AmrLevel& fine_level = m_amrlevel;
                  AmrLevel& crse_level = m_amrlevel.parent->getLevel(ilev_crse);
#ifdef USE_PERILLA_PTHREADS
//                  if(perilla::isMasterThread())
#endif
                  {
                      geom_fine = &fine_level.geom;
                      geom_crse = &crse_level.geom;
                  }
                  StateData& statedata_crse = crse_level.state[index];
                  StateData& statedata_fine = fine_level.state[index];
#ifdef USE_PERILLA_PTHREADS
//                  perilla::syncAllThreads();
//                  if(perilla::isMasterThread())
#endif
                  {
                      statedata_crse.getData(smf_crse,stime_crse,time);
                      statedata_fine.getData(smf_fine,stime_fine,time);
                      physbcf_crse = new StateDataPhysBCFunct(statedata_crse,scomp,*geom_crse);
                      physbcf_fine = new StateDataPhysBCFunct(statedata_fine,scomp,*geom_fine);
		      stateDataList.push_back(physbcf_crse);
		      stateDataList.push_back(physbcf_fine);
                  }
#ifdef USE_PERILLA_PTHREADS
//                perilla::syncAllThreads();
#endif
                  const StateDescriptor& desc = AmrLevel::desc_lst[index];
                  int ngrow = m_fabs.nGrow();
                  if (ngrow > 0 || m_fabs.getBDKey() != smf_fine[0]->getBDKey())
                  {
#ifdef USE_PERILLA_PTHREADS
//                    if(perilla::isMasterThread())
#endif
                      {
                          InterpolaterBoxCoarsener coarsener = desc.interp(scomp)->BoxCoarsener(crse_level.fineRatio());
                          Box fdomain = geom_fine->Domain();
                          fdomain.convert(m_fabs.boxArray().ixType());
                          Box fdomain_g(fdomain);
                          for (int i = 0; i < BL_SPACEDIM; ++i) {
                             if (geom_fine->isPeriodic(i)) {
                                 fdomain_g.grow(i,ngrow);
                             }
                          }
			  Box c_dom= amrex::coarsen(geom_fine->Domain(), m_amrlevel.crse_ratio);
                          m_fpc = &FabArrayBase::TheFPinfo(*(smf_fine[0]), m_fabs, fdomain_g, IntVect(ngrow), coarsener, c_dom);
                      }
#ifdef USE_PERILLA_PTHREADS
//                    perilla::syncAllThreads();
#endif
                      if (!m_fpc->ba_crse_patch.empty())
                      {
#ifdef USE_PERILLA_PTHREADS
                          if(perilla::isMasterThread())
#endif
                          {
                             m_mf_crse_patch = new MultiFab(m_fpc->ba_crse_patch, m_fpc->dm_crse_patch, ncomp, 0);
		             mfList.push_back(m_mf_crse_patch);
                             //m_mf_crse_patch->initVal(); // for Perilla NUMA
                             BL_ASSERT(scomp+ncomp <= smf_crse[0]->nComp());
                             BL_ASSERT(dcomp+ncomp <= m_mf_crse_patch->nComp());
                             BL_ASSERT(smf_crse.size() == stime_crse.size());
                             BL_ASSERT(smf_crse.size() != 0);
                          }
#ifdef USE_PERILLA_PTHREADS
//                        perilla::syncAllThreads();
#endif
                          if (iter == 1)
                          {
#ifdef USE_PERILLA_PTHREADS
//                            if(perilla::isMasterThread())
#endif
                              {
                                  m_rg_crse_patch = new RegionGraph(m_mf_crse_patch->IndexArray().size());
                                  m_rg_crse_patch->isDepGraph = true;
                                  csrcGraph = new RegionGraph(smf_crse[0]->IndexArray().size());
		    		  regionList.push_back(m_rg_crse_patch);
		    		  regionList.push_back(csrcGraph);
                              }
#ifdef USE_PERILLA_PTHREADS
//                            perilla::syncAllThreads();
#endif
#if 1
                              Perilla::multifabExtractCopyAssoc( m_rg_crse_patch, csrcGraph, *m_mf_crse_patch, *(smf_crse[0]), ncomp, m_mf_crse_patch->nGrow(), 0,geom_crse->periodicity());
#endif
#ifdef USE_PERILLA_PTHREADS
//                              perilla::syncAllThreads();
//                            if(perilla::isMasterThread())
#endif
                              {
                                  m_amrlevel.parent->graphArray[level].push_back(m_rg_crse_patch);
                                  m_amrlevel.parent->graphArray[level].push_back(csrcGraph);
                              }
                          }
                          else if (iter > 1)
                          {
#if 1
                              //BL_ASSERT(smf_crse[0].boxArray() == smf_crse[1].boxArray());
                              //PArray<MultiFab> raii(PArrayManage);
                              //MultiFab * dmf;

                              if (m_mf_crse_patch->boxArray() == smf_crse[0]->boxArray())
                                {
                                  //dmf = m_mf_crse_patch;
                                  m_rg_crse_patch = new RegionGraph(m_mf_crse_patch->IndexArray().size());

                                  //std::cout<< " level " << level  << " rg_crs_ptch ID " << m_rg_crse_patch->graphID << std::endl;

                                  Perilla::multifabBuildFabCon(m_rg_crse_patch, *m_mf_crse_patch, geom->periodicity());
                                  m_amrlevel.parent->graphArray[level].push_back(m_rg_crse_patch);
		    		  regionList.push_back(m_rg_crse_patch);
                                }
                              else
                              {
#ifdef USE_PERILLA_PTHREADS
//                                if(perilla::isMasterThread())
#endif
                                  {
                                      //dmf = raii.push_back(new MultiFab(smf_crse[0].boxArray(), ncomp, 0));
                                      dmf = new MultiFab(smf_crse[0]->boxArray(), smf_crse[0]->DistributionMap(), ncomp, 0);
                                      //dmf->initVal(); // for Perilla NUMA
                                      m_rg_crse_patch = new RegionGraph(m_mf_crse_patch->IndexArray().size());
                                      m_rg_crse_patch->isDepGraph = true;
                                      csrcGraph = new RegionGraph(dmf->IndexArray().size());
                                      csrcGraph->buildTileArray(*dmf);
		    		      regionList.push_back(m_rg_crse_patch);
		    		      regionList.push_back(csrcGraph);
		    		      mfList.push_back(dmf);
                                  }
#ifdef USE_PERILLA_PTHREADS
//                                perilla::syncAllThreads();
#endif

#if 1
                                  Perilla::multifabExtractCopyAssoc( m_rg_crse_patch, csrcGraph, *m_mf_crse_patch, *dmf, ncomp, m_mf_crse_patch->nGrow(), 0, geom_crse->periodicity());
#endif
#ifdef USE_PERILLA_PTHREADS
//                                perilla::syncAllThreads();
//                                if(perilla::isMasterThread())
#endif
                                  {
                                      m_amrlevel.parent->graphArray[level].push_back(m_rg_crse_patch);
                                      m_amrlevel.parent->graphArray[level].push_back(csrcGraph);
                                  }
                              }
#endif
                           }
                           else
                           {
                             // BoxLib::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
                           }
                        }
                    }
#if 1
                  BL_ASSERT(scomp+ncomp <= smf_fine[0]->nComp());
                  BL_ASSERT(dcomp+ncomp <= m_fabs.nComp());
                  BL_ASSERT(smf_fine.size() == stime_fine.size());
                  BL_ASSERT(smf_fine.size() != 0);

                  if(true) // it will always be the case because same level comm and time will be available
                  {
#ifdef USE_PERILLA_PTHREADS
//                    if(perilla::isMasterThread())
#endif
                      {
                          dmff = new MultiFab(smf_fine[0]->boxArray(), smf_fine[0]->DistributionMap(), ncomp, 0);
                          destGraph = new RegionGraph(m_fabs.IndexArray().size());
                          //fsrcGraph = new RegionGraph(smf_fine[0]->IndexArray().size());
                          fsrcGraph = new RegionGraph(dmff->IndexArray().size());
			  regionList.push_back(destGraph);
			  regionList.push_back(fsrcGraph);
			  mfList.push_back(dmff);

                          if(m_rg_crse_patch != 0)
                          {
                             destGraph->srcLinkGraph = m_rg_crse_patch;
                             //for(int lfi=0; lfi < destGraph->numTasks; lfi++ )
                             {
                                for (MFIter mfi(*(m_mf_crse_patch),false); mfi.isValid(); ++mfi)
                                {
                                  int li = mfi.LocalIndex();
                                  int gi = m_fpc->dst_idxs[li];
                                  //if(gi == m_mf_crse_patch->IndexArray()[li])
                                    {
                                      int lfi = m_fabs.localindex(gi);
                                      destGraph->task[lfi]->depTasksCompleted = false;
                                      destGraph->task[lfi]->depTaskIDs.push_back(li);
                                    }
                                }
                             }
                          }
                      }
#ifdef USE_PERILLA_PTHREADS
//                    perilla::syncAllThreads();
#endif

                      //if(level == 2)
                      //std::cout<< "Sending In "<<destGraph<<" "<< fsrcGraph << " myP " << ParallelDescriptor::MyProc()<<std::endl;                                      
#if 1
                      Perilla::multifabExtractCopyAssoc( destGraph, fsrcGraph, m_fabs, *(dmff), ncomp, m_fabs.nGrow(), 0, geom_fine->periodicity());
#endif
#ifdef USE_PERILLA_PTHREADS
//                    perilla::syncAllThreads();
//                    if(perilla::isMasterThread())
#endif
                      {
                          m_amrlevel.parent->graphArray[level].push_back(destGraph);
                          m_amrlevel.parent->graphArray[level].push_back(fsrcGraph);
                      }
                  }
                  else if (smf_fine.size() == 2)
                  {
                      if (m_fabs.boxArray() == smf_fine[0]->boxArray())
                      {
                          //dmf = &m_fabs;
                          destGraph = new RegionGraph(m_fabs.IndexArray().size());
                          Perilla::multifabBuildFabCon(destGraph, m_fabs, geom->periodicity());
                          m_amrlevel.parent->graphArray[level].push_back(destGraph);
			  regionList.push_back(destGraph);
                      }
                      else
                      {
#ifdef USE_PERILLA_PTHREADS
//                        if(perilla::isMasterThread())
#endif
                          {
                              //dmf = raii.push_back(new MultiFab(smf_fine[0].boxArray(), ncomp, 0));
                              dmff = new MultiFab(smf_fine[0]->boxArray(), smf_fine[0]->DistributionMap(), ncomp, 0);
                              //dmff->initVal(); // for Perilla NUMA
                              destGraph = new RegionGraph(m_fabs.IndexArray().size());
                              fsrcGraph = new RegionGraph(dmff->IndexArray().size());
                              fsrcGraph->buildTileArray(*dmff);
			      regionList.push_back(destGraph);
			      regionList.push_back(fsrcGraph);
			      mfList.push_back(dmff);
                          }
#ifdef USE_PERILLA_PTHREADS
//                        perilla::syncAllThreads();
#endif
#if 1
                          Perilla::multifabExtractCopyAssoc( destGraph, fsrcGraph, m_fabs, *dmff, ncomp, m_fabs.nGrow(), 0, geom_fine->periodicity());
#endif
#ifdef USE_PERILLA_PTHREADS
//                        perilla::syncAllThreads();
//                        if(perilla::isMasterThread())
#endif
                          {
                              m_amrlevel.parent->graphArray[level].push_back(destGraph);
                              m_amrlevel.parent->graphArray[level].push_back(fsrcGraph);
                          }
                       }
                  }
                  else
                  {
                      amrex::Abort("FillPatchSingleLevel: high-order interpolation in time not implemented yet");
                  }
#endif
                  //-------------------- FillFromTwoLevels initialization completed
                } // if(level==1 OR ProperlyNested)
              else
                {
                  amrex::Abort("initFillPatch: level is not properly nested");
                }
          }

        DComp += NComp;
    }
#if 1
#ifdef USE_PERILLA_PTHREADS
//    perilla::syncAllThreads();
//    if(perilla::isMasterThread())
#endif
    {
        destGraph->buildTileArray(m_fabs);
        destGraph->buildTileArray_gtbx(m_leveldata,boxGrow);
    }
#ifdef USE_PERILLA_PTHREADS
//    perilla::syncAllThreads();
#endif
#endif
}

FillPatchIterator::FillPatchIterator (AmrLevel& amrlevel,
                                      MultiFab& leveldata,
                                      int       boxGrow,
                                      Real      time,
                                      int       index,
                                      int       scomp,
                                      int       ncomp,
                                      int       f)
    :
    MFIter(leveldata),
    m_amrlevel(amrlevel),
    m_leveldata(leveldata),
    m_ncomp(ncomp),
    physbcf(0),
    physbcf_crse(0),
    physbcf_fine(0),
    destGraph(0),
    fsrcGraph(0),
    csrcGraph(0),
    m_rg_crse_patch(NULL),
    //raii(PArrayManage)
    dmf(NULL),
    dmff(NULL)
{
#if 1
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(AmrLevel::desc_lst[index].inRange(scomp,ncomp));
    BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

    //InitializePush(boxGrow,time,index,scomp,ncomp,f,tid);

#ifdef BL_USE_TEAM
    ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

#endif
}

//end USE_PERILLA
#endif

}

