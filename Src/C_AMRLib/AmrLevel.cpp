
#include <winstd.H>
#include <sstream>

#include <unistd.h>

#include <AmrLevel.H>
#include <Derive.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <FillPatchUtil.H>
#include <ParmParse.H>
#include <BLProfiler.H>

DescriptorList AmrLevel::desc_lst;
DeriveList     AmrLevel::derive_lst;

#ifdef USE_SLABSTAT
SlabStatList   AmrLevel::slabstat_lst;
#endif

void
AmrLevel::postCoarseTimeStep (Real time)
{}

#ifdef USE_SLABSTAT
SlabStatList&
AmrLevel::get_slabstat_lst ()
{
    return slabstat_lst;
}
#endif

void
AmrLevel::set_preferred_boundary_values (MultiFab& S,
                                         int       state_index,
                                         int       scomp,
                                         int       dcomp,
                                         int       ncomp,
                                         Real      time) const
{}

DeriveList&
AmrLevel::get_derive_lst ()
{
    return derive_lst;
}

void
AmrLevel::manual_tags_placement (TagBoxArray&    tags,
                                 Array<IntVect>& bf_lev)
{}

AmrLevel::AmrLevel ()
{
   parent = 0;
   level = -1;
}

AmrLevel::AmrLevel (Amr&            papa,
                    int             lev,
                    const Geometry& level_geom,
                    const BoxArray& ba,
                    Real            time)
    :
    geom(level_geom),
    grids(ba)
#ifdef USE_PARTICLES
    ,particle_grids(ba)
#endif
{
    BL_PROFILE("AmrLevel::AmrLevel()");
    level  = lev;
    parent = &papa;

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

    // Note that this creates a distribution map associated with grids.
    for (int i = 0; i < state.size(); i++)
    {
        state[i].define(geom.Domain(),
                        grids,
                        desc_lst[i],
                        time,
                        parent->dtLevel(lev));
    }

    if (Amr::useFixedCoarseGrids) constructAreaNotToTag();

#ifdef USE_PARTICLES
    // Note: it is important to call make_particle_dmap *after* the state
    //       has been defined because it makes use of the state's DistributionMap
    make_particle_dmap();
#endif

    finishConstructor();
}

AmrLevel::AmrLevel (Amr&            papa,
                    int             lev,
                    const Geometry& level_geom,
                    const BoxArray& ba,
		    const DistributionMapping& dm,
                    Real            time)
    :
    geom(level_geom),
    grids(ba)
#ifdef USE_PARTICLES
    ,particle_grids(ba)
#endif
{
    BL_PROFILE("AmrLevel::AmrLevel(dm)");
    level  = lev;
    parent = &papa;

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

    // Note that this creates a distribution map associated with grids.
    for (int i = 0; i < state.size(); i++)
    {
        state[i].define(geom.Domain(),
                        grids,
			dm,
                        desc_lst[i],
                        time,
                        parent->dtLevel(lev));
    }

    if (Amr::useFixedCoarseGrids) constructAreaNotToTag();

#ifdef USE_PARTICLES
    // Note: it is important to call make_particle_dmap *after* the state
    //       has been defined because it makes use of the state's DistributionMap
    make_particle_dmap();
#endif

    finishConstructor();
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
        BoxLib::readBoxArray(grids, is, bReadSpecial);
    }
    else
    {
        grids.readFrom(is);
    }

    int nstate;
    is >> nstate;
    int ndesc = desc_lst.size();

    Array<int> state_in_checkpoint(ndesc, 1);
    if (ndesc > nstate) {
	set_state_in_checkpoint(state_in_checkpoint);
    } else {
	BL_ASSERT(nstate == ndesc);
    }

    state.resize(ndesc);
    for (int i = 0; i < ndesc; i++)
    {
	if (state_in_checkpoint[i]) {
	    state[i].restart(is, desc_lst[i], papa.theRestartFile(), bReadSpecial);
	}
    }
 
    if (Amr::useFixedCoarseGrids) constructAreaNotToTag();

#ifdef USE_PARTICLES
    // Note: it is important to call make_particle_dmap *after* the state
    //       has been defined because it makes use of the state's DistributionMap
    particle_grids = grids;
    make_particle_dmap();
#endif

    finishConstructor();
}

void
AmrLevel::set_state_in_checkpoint (Array<int>& state_in_checkpoint)
{
    BoxLib::Error("Class derived AmrLevel has to handle this!");
}

void
AmrLevel::finishConstructor () {}

#ifdef USE_PARTICLES
void
AmrLevel::make_particle_dmap ()
{
    // Here we create particle_grids and make a distribution map for it
    // Right now particle_grids is identical to grids, but the whole point is that 
    // particle_grids can be optimized for distributing the particle work. 

    // take a shortcut if possible
    if (grids == particle_grids) {
	particle_dmap = get_new_data(0).DistributionMap();
	return;
    }

    Array<int> ParticleProcMap;
    ParticleProcMap.resize(particle_grids.size()+1); // +1 is a historical thing

    for (int i = 0; i <= particle_grids.size(); i++)
        ParticleProcMap[i] = -1;

    // Warning: O(N^2)!
    for (int j = 0; j < grids.size(); j++)
    {
        const int who = get_new_data(0).DistributionMap()[j];
        for (int i = 0; i < particle_grids.size(); i++)
        {
            if (grids[j].contains(particle_grids[i]))
            {
                ParticleProcMap[i] = who;
            }
        }
    }

    // Don't forget the last entry!
    ParticleProcMap[particle_grids.size()] = ParallelDescriptor::MyProc();

    // Sanity check that all grids got assigned to processors
    for (int i = 0; i <= particle_grids.size(); i++)
        if (ParticleProcMap[i] == -1)
            BoxLib::Error("Didn't assign every particle_grids box to a processor!!");

    // Different DistributionMappings must have different numbers of boxes 
    if (grids.size() != particle_grids.size())
    {
        const bool put_in_cache = true;
        particle_dmap = DistributionMapping(ParticleProcMap,put_in_cache);
    }
    else
    {
        if (grids != particle_grids)
        {
           // Oops -- can't handle this 
            BoxLib::Error("grids != particle_grids but they have the same number of boxes");
        }
        else
        {
           // Just copy the grids distribution map to the particle_grids distribution map
           particle_dmap = get_new_data(0).DistributionMap();
        }
    }
}
#endif

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
AmrLevel::isStateVariable (const std::string& name,
                           int&           typ,
                           int&            n)
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
AmrLevel::countCells () const
{
    const int N = grids.size();

    long cnt = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:cnt)
#endif
    for (int i = 0; i < N; i++)
    {
        cnt += grids[i].numPts();
    }

    return cnt;
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
    std::string Level = BoxLib::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
    {
        FullPath += '/';
    }
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier("AmrLevel::checkPoint::dir");

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
        // There is only one MultiFab written out at each level in HyperCLaw.
        //
        std::string PathNameInHdr = BoxLib::Concatenate(Level    + "/SD_", i, 1);
        std::string FullPathName  = BoxLib::Concatenate(FullPath + "/SD_", i, 1);

        state[i].checkPoint(PathNameInHdr, FullPathName, os, how, dump_old);
    }
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
AmrLevel::get_data (int  state_indx,
                    Real time)
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

    BoxLib::Error("get_data: invalid time");
    static MultiFab bogus;
    return bogus;
}

const BoxArray&
AmrLevel::getEdgeBoxArray (int dir) const
{
    BL_ASSERT(dir >=0 && dir < BL_SPACEDIM);
    if (edge_grids[dir].empty()) {
	edge_grids[dir] = grids;
	edge_grids[dir].surroundingNodes(dir);
    }
    return edge_grids[dir];
}

const BoxArray&
AmrLevel::getNodalBoxArray () const
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
                                      int       index,
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
    BL_ASSERT(AmrLevel::desc_lst[index].inRange(scomp,ncomp));
    BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

    Initialize(boxGrow,time,index,scomp,ncomp);

#if BL_USE_TEAM
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
                                     int           index,
                                     int           scomp,
                                     int           ncomp,
                                     Interpolater* mapper)
{
    BL_PROFILE("FillPatchIteratorHelper::Initialize()");

    BL_ASSERT(mapper);
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(AmrLevel::desc_lst[index].inRange(scomp,ncomp));
    BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

    m_map          = mapper;
    m_time         = time;
    m_growsize     = boxGrow;
    m_index        = index;
    m_scomp        = scomp;
    m_ncomp        = ncomp;
    m_FixUpCorners = NeedToTouchUpPhysCorners(m_amrlevel.geom);

    const int         MyProc     = ParallelDescriptor::MyProc();
    PArray<AmrLevel>& amrLevels  = m_amrlevel.parent->getAmrLevels();
    const AmrLevel&   topLevel   = amrLevels[m_amrlevel.level];
    const Box&        topPDomain = topLevel.state[m_index].getDomain();
    const IndexType&  boxType    = m_leveldata.boxArray().ixType();
    const bool        extrap     = AmrLevel::desc_lst[m_index].extrap();
    //
    // Check that the interpolaters are identical.
    //
    BL_ASSERT(AmrLevel::desc_lst[m_index].identicalInterps(scomp,ncomp));

    for (int l = 0; l <= m_amrlevel.level; ++l)
    {
        amrLevels[l].state[m_index].RegisterData(m_mfcd, m_mfid[l]);
    }
    for (int i = 0, N = m_leveldata.boxArray().size(); i < N; ++i)
    {
        //
        // A couple typedefs we'll use in the next code segment.
        //
        typedef std::map<int,Array<Array<Box> > >::value_type IntAABoxMapValType;

        typedef std::map<int,Array<Array<Array<FillBoxId> > > >::value_type IntAAAFBIDMapValType;

        if (m_leveldata.DistributionMap()[i] != MyProc) continue;
        //
        // Insert with a hint since the indices are ordered lowest to highest.
        //
        IntAAAFBIDMapValType v1(i,Array<Array<Array<FillBoxId> > >());

        m_fbid.insert(m_fbid.end(),v1)->second.resize(m_amrlevel.level+1);

        IntAABoxMapValType v2(i,Array<Array<Box> >());

        m_fbox.insert(m_fbox.end(),v2)->second.resize(m_amrlevel.level+1);
        m_cbox.insert(m_cbox.end(),v2)->second.resize(m_amrlevel.level+1);

        m_ba.insert(m_ba.end(),std::map<int,Box>::value_type(i,BoxLib::grow(m_leveldata.boxArray()[i],m_growsize)));
    }

    BoxList        tempUnfillable(boxType);
    BoxList        unfillableThisLevel(boxType);
    Array<Box>     unfilledThisLevel;
    Array<Box>     crse_boxes;
    Array<IntVect> pshifts(27);

    for (std::map<int,Box>::const_iterator it = m_ba.begin(), End = m_ba.end();
         it != End;
         ++it)
    {
        const int  idx = it->first;
        const Box& box = it->second;

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

                for (Array<IntVect>::const_iterator pit = pshifts.begin(),
                         End = pshifts.end();
                     pit != End;
                     ++pit)
                {
                    const IntVect& iv = *pit;

                    Box shbox = box + iv;
                    shbox    &= topPDomain;

                    if (boxType.nodeCentered())
                    {
                        for (int dir = 0; dir < BL_SPACEDIM; dir++)
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
                        BoxList bl = BoxLib::boxDiff(shbox,box);

                        unfilledThisLevel.insert(unfilledThisLevel.end(), bl.begin(), bl.end());
                    }
                }
            }
        }

	// cells outside physical boundaries are not included in unfilledThisLevel

        bool Done = false;

        Array< Array<Box> >&                TheCrseBoxes = m_cbox[idx];
        Array< Array<Box> >&                TheFineBoxes = m_fbox[idx];
        Array< Array< Array<FillBoxId> > >& TheFBIDs     = m_fbid[idx];

        for (int l = m_amrlevel.level; l >= 0 && !Done; --l)
        {
            unfillableThisLevel.clear();

            AmrLevel&       theAmrLevel = amrLevels[l];
            StateData&      theState    = theAmrLevel.state[m_index];
            const Box&      thePDomain  = theState.getDomain();
            const Geometry& theGeom     = theAmrLevel.geom;
            const bool      is_periodic = theGeom.isAnyPeriodic();
            const IntVect&  fine_ratio  = theAmrLevel.fine_ratio;
            Array<Box>&     FineBoxes   = TheFineBoxes[l];
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

            for (Array<Box>::const_iterator fit = FineBoxes.begin(),
                     End = FineBoxes.end();
                 fit != End;
                 ++fit)
            {
                crse_boxes.push_back(*fit);

                if (l != m_amrlevel.level)
                {
                    const Box& cbox = m_map->CoarseBox(*fit,fine_ratio);

		    crse_boxes.back() = cbox;

                    if (is_periodic && !thePDomain.contains(cbox))
                    {
                        theGeom.periodicShift(thePDomain,cbox,pshifts);

                        for (Array<IntVect>::const_iterator pit = pshifts.begin(),
                                 End = pshifts.end();
                             pit != End;
                             ++pit)
                        {
                            const IntVect& iv = *pit;

                            Box shbox = cbox + iv;
                            shbox    &= thePDomain;

                            if (boxType.nodeCentered())
                            {
                                for (int dir = 0; dir < BL_SPACEDIM; dir++)
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

            Array< Array<FillBoxId> >& FBIDs     = TheFBIDs[l];
            Array<Box>&                CrseBoxes = TheCrseBoxes[l];

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
                               int  index,
                               int  scomp,
                               int  ncomp)
{
    BL_PROFILE("FillPatchIterator::Initialize");

    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >= 1);
    BL_ASSERT(0 <= index && index < AmrLevel::desc_lst.size());

    const StateDescriptor& desc = AmrLevel::desc_lst[index];

    m_ncomp = ncomp;
    m_range = desc.sameInterps(scomp,ncomp);

    m_fabs.define(m_leveldata.boxArray(),m_ncomp,boxGrow,Fab_allocate);

    BL_ASSERT(m_leveldata.DistributionMap() == m_fabs.DistributionMap());

    const IndexType& boxType = m_leveldata.boxArray().ixType();
    const int level = m_amrlevel.level;

    for (int i = 0, DComp = 0; i < m_range.size(); i++)
    {
        const int SComp = m_range[i].first;
        const int NComp = m_range[i].second;

	if (level == 0)
	{
	    FillFromLevel0(time, index, SComp, DComp, NComp);
	}
	else
	{
	    if (level == 1 || 
		BoxLib::ProperlyNested(m_amrlevel.crse_ratio,
				       m_amrlevel.parent->blockingFactor(m_amrlevel.level),
				       boxGrow, boxType, desc.interp(SComp)))
	    {
		FillFromTwoLevels(time, index, SComp, DComp, NComp);
	    } else {
		static bool first = true;
		if (first) {
		    first = false;
		    if (ParallelDescriptor::IOProcessor()) {
			int new_blocking_factor = 2*m_amrlevel.parent->blockingFactor(m_amrlevel.level);
			for (int i = 0; i < 10; ++i) {
			    if (BoxLib::ProperlyNested(m_amrlevel.crse_ratio,
						       new_blocking_factor,
						       boxGrow, boxType, desc.interp(SComp))) {
				break;
			    } else {
				new_blocking_factor *= 2;
			    }
			}
			std::cout << "WARNING: Grids are not properly nested.  We might have to use\n"
				  << "         two coarse levels to do fillpatch.  Consider using\n";
			if (new_blocking_factor < 128) {
			    std::cout << "         amr.blocking_factor=" << new_blocking_factor;
			} else {
			    std::cout << "         larger amr.blocking_factor. ";
			}
			std::cout << std::endl;
		    }
		}

		FillPatchIteratorHelper* fph = 0;
		fph = new FillPatchIteratorHelper(m_amrlevel,
						  m_leveldata,
						  boxGrow,
						  time,
						  index,
						  SComp,
						  NComp,
						  desc.interp(SComp));
		
#ifdef CRSEGRNDOMP
#ifdef _OPENMP
#pragma omp parallel
#endif
#endif
		for (MFIter mfi(m_fabs); mfi.isValid(); ++mfi)
		{
		    fph->fill(m_fabs[mfi],DComp,mfi.index());
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
                                             index,
                                             scomp,
                                             0,
                                             ncomp,
                                             time);
}

void
FillPatchIterator::FillFromLevel0 (Real time, int index, int scomp, int dcomp, int ncomp)
{
    BL_ASSERT(m_amrlevel.level == 0);

    StateData& statedata = m_amrlevel.state[index];

    PArray<MultiFab> smf;
    std::vector<Real> stime;
    statedata.getData(smf,stime,time);

    const Geometry& geom = m_amrlevel.geom;

    StateDataPhysBCFunct physbcf(statedata,scomp,geom);

    BoxLib::FillPatchSingleLevel (m_fabs, time, smf, stime, scomp, dcomp, ncomp, geom, physbcf);
}

void
FillPatchIterator::FillFromTwoLevels (Real time, int index, int scomp, int dcomp, int ncomp)
{
    int ilev_fine = m_amrlevel.level;
    int ilev_crse = ilev_fine-1;

    BL_ASSERT(ilev_crse >= 0);

    AmrLevel& fine_level = m_amrlevel;
    AmrLevel& crse_level = m_amrlevel.parent->getLevel(ilev_crse);

    const Geometry& geom_fine = fine_level.geom;
    const Geometry& geom_crse = crse_level.geom;
    
    PArray<MultiFab> smf_crse;
    std::vector<Real> stime_crse;
    StateData& statedata_crse = crse_level.state[index];
    statedata_crse.getData(smf_crse,stime_crse,time);
    StateDataPhysBCFunct physbcf_crse(statedata_crse,scomp,geom_crse);

    PArray<MultiFab> smf_fine;
    std::vector<Real> stime_fine;
    StateData& statedata_fine = fine_level.state[index];
    statedata_fine.getData(smf_fine,stime_fine,time);
    StateDataPhysBCFunct physbcf_fine(statedata_fine,scomp,geom_fine);

    const StateDescriptor& desc = AmrLevel::desc_lst[index];

    BoxLib::FillPatchTwoLevels(m_fabs, time, 
			       smf_crse, stime_crse, 
			       smf_fine, stime_fine,
			       scomp, dcomp, ncomp, 
			       geom_crse, geom_fine,
			       physbcf_crse, physbcf_fine,
			       crse_level.fineRatio(), 
			       desc.interp(scomp), desc.getBCs());
}

static
bool
HasPhysBndry (const Box&      b,
              const Box&      dmn,
              const Geometry& geom)
{
    for (int i = 0; i < BL_SPACEDIM; i++)
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

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        if (!TheGeom.isPeriodic(dir))
        {
            const int lo = ProbDomain.smallEnd(dir) - fab.box().smallEnd(dir);
            const int hi = fab.box().bigEnd(dir)    - ProbDomain.bigEnd(dir);
            if (lo > 0) GrownDomain.growLo(dir,lo);
            if (hi > 0) GrownDomain.growHi(dir,hi);
        }
    }

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
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

    Array< PArray<FArrayBox> >          cfab(m_amrlevel.level+1);
    Array< Array<Box> >&                TheCrseBoxes = m_cbox[idx];
    Array< Array<Box> >&                TheFineBoxes = m_fbox[idx];
    Array< Array< Array<FillBoxId> > >& TheFBIDs     = m_fbid[idx];
    const bool                          extrap       = AmrLevel::desc_lst[m_index].extrap();
    PArray<AmrLevel>&                   amrLevels    = m_amrlevel.parent->getAmrLevels();
    //
    // Build all coarse fabs from which we'll interpolate and
    // fill them with coarse data as best we can.
    //
    for (int l = 0; l <= m_amrlevel.level; l++)
    {
        StateData&                       TheState  = amrLevels[l].state[m_index];
        const Array<Box>&                CrseBoxes = TheCrseBoxes[l];
        PArray<FArrayBox>&               CrseFabs  = cfab[l];
        const Array< Array<FillBoxId> >& FBIDs     = TheFBIDs[l];
        const int                        NC        = CrseBoxes.size();

        CrseFabs.resize(NC,PArrayManage);

        for (int i = 0; i < NC; i++)
        {
            BL_ASSERT(CrseBoxes[i].ok());
            CrseFabs.set(i, new FArrayBox(CrseBoxes[i],m_ncomp));
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
				   CrseFabs[i],
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
        PArray<FArrayBox>& CrseFabs   = cfab[l];
        AmrLevel&          TheLevel   = amrLevels[l];
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
                FArrayBox& dstfab = CrseFabs[i];

                if (ThePDomain.contains(dstfab.box())) continue;

                Array<IntVect> pshifts(27);

                TheLevel.geom.periodicShift(ThePDomain,dstfab.box(),pshifts);

                for (Array<IntVect>::const_iterator pit = pshifts.begin(),
                         End = pshifts.end();
                     pit != End;
                     ++pit)
                {
                    const IntVect& iv = *pit;

                    Box fullsrcbox = dstfab.box() + iv;
                    fullsrcbox    &= ThePDomain;

                    for (int j = 0; j < NC; j++)
                    {
                        const FArrayBox& srcfab = CrseFabs[j];
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
            if ( ! ThePDomain.contains(CrseFabs[i].box()))
            {
                TheLevel.setPhysBoundaryValues(CrseFabs[i],
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
                FixUpPhysCorners(CrseFabs[i],TheLevel,m_index,m_time,m_scomp,0,m_ncomp);
            }
        }
        //
        // Interpolate up to next level.
        //
        AmrLevel&           crseAmrLevel  = amrLevels[l];
        AmrLevel&           fineAmrLevel  = amrLevels[l+1];
        const IntVect&      fine_ratio    = crseAmrLevel.fine_ratio;
        const Array<Box>&   FineBoxes     = TheFineBoxes[l];
        StateData&          fState        = fineAmrLevel.state[m_index];
        const Box&          fDomain       = fState.getDomain();
        PArray<FArrayBox>&  FinerCrseFabs = cfab[l+1];
        const Array<BCRec>& theBCs        = AmrLevel::desc_lst[m_index].getBCs();
        const int           NF            = FineBoxes.size();

        for (int ifine = 0; ifine < NF; ++ifine)
        {
            Array<BCRec> bcr(m_ncomp);
            FArrayBox    finefab(FineBoxes[ifine],m_ncomp);
            FArrayBox    crsefab(m_map->CoarseBox(finefab.box(),fine_ratio),m_ncomp);
            //
            // Fill crsefab from m_cbox via copy on intersect.
            //
            for (int j = 0; j < NC; j++) {
                crsefab.copy(CrseFabs[j]);
	    }
            //
            // Get boundary conditions for the fine patch.
            //
            BoxLib::setBC(finefab.box(),
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
		FinerCrseFabs[j].copy(finefab);
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
    PArray<FArrayBox>& FinestCrseFabs = cfab[m_amrlevel.level];
    //
    // Copy intersect coarse into destination fab.
    //
    for (int i = 0, N = FinestCrseFabs.size(); i < N; ++i) {
        fab.copy(FinestCrseFabs[i],0,dcomp,m_ncomp);
    }

    if (FineGeom.isAnyPeriodic() && !FineDomain.contains(fab.box()))
    {
        Array<IntVect> pshifts(27);

        FineGeom.periodicShift(FineDomain,fab.box(),pshifts);

        for (int i = 0, N = FinestCrseFabs.size(); i < N; i++)
        {
            for (Array<IntVect>::const_iterator pit = pshifts.begin(),
                     End = pshifts.end();
                 pit != End;
                 ++pit)
            {
                const IntVect& iv = *pit;

                fab.shift(iv);

                Box src_dst = FinestCrseFabs[i].box() & fab.box();
                src_dst    &= FineDomain;

                if (src_dst.ok())
                    fab.copy(FinestCrseFabs[i],src_dst,0,src_dst,dcomp,m_ncomp);

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

FillPatchIterator::~FillPatchIterator () {}

void
AmrLevel::FillCoarsePatch (MultiFab& mf,
                           int       dcomp,
                           Real      time,
                           int       index,
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
    BL_ASSERT(0 <= index && index < desc_lst.size());

    int                     DComp   = dcomp;
    const StateDescriptor&  desc    = desc_lst[index];
    const Box&              pdomain = state[index].getDomain();
    const BoxArray&         mf_BA   = mf.boxArray();
    AmrLevel&               clev    = parent->getLevel(level-1);

    std::vector< std::pair<int,int> > ranges  = desc.sameInterps(scomp,ncomp);

    BL_ASSERT(desc.inRange(scomp, ncomp));

    for (int i = 0; i < ranges.size(); i++)
    {
        const int     SComp  = ranges[i].first;
        const int     NComp  = ranges[i].second;
        Interpolater* mapper = desc.interp(SComp);

        BoxArray crseBA(mf_BA.size());
        
        for (int j = 0, N = crseBA.size(); j < N; ++j)
        {
            BL_ASSERT(mf_BA[j].ixType() == desc.getType());

            crseBA.set(j,mapper->CoarseBox(BoxLib::grow(mf_BA[j],nghost),crse_ratio));
        }

	MultiFab crseMF(crseBA,NComp,0);

	if ( level == 1 
	     || BoxLib::ProperlyNested(crse_ratio, parent->blockingFactor(level),
				       nghost, mf_BA.ixType(), mapper) )
	{
	    StateData& statedata = clev.state[index];
	    
	    PArray<MultiFab> smf;
	    std::vector<Real> stime;
	    statedata.getData(smf,stime,time);

	    const Geometry& geom = clev.geom;

	    StateDataPhysBCFunct physbcf(statedata,SComp,geom);

	    BoxLib::FillPatchSingleLevel(crseMF,time,smf,stime,SComp,0,NComp,geom,physbcf);
	}
	else
	{
	    FillPatch(clev,crseMF,0,time,index,SComp,NComp,0);
	}

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	    const Box& dbx = BoxLib::grow(mfi.validbox(),nghost);
	    
	    Array<BCRec> bcr(ncomp);
	    
	    BoxLib::setBC(dbx,pdomain,SComp,0,NComp,desc.getBCs(),bcr);
	    
	    mapper->interp(crseMF[mfi],
			   0,
			   mf[mfi],
			   DComp,
			   NComp,
			   dbx,
			   crse_ratio,
			   clev.geom,
			   geom,
			   bcr,
			   SComp,
			   index);
	}

        DComp += NComp;
    }
}

MultiFab*
AmrLevel::derive (const std::string& name,
                  Real           time,
                  int            ngrow)
{
    BL_ASSERT(ngrow >= 0);

    MultiFab* mf = 0;

    int index, scomp, ncomp;

    if (isStateVariable(name, index, scomp))
    {
        mf = new MultiFab(state[index].boxArray(), 1, ngrow);
        FillPatch(*this,*mf,ngrow,time,index,scomp,1);
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

        MultiFab srcMF(srcBA, rec->numState(), ngrow_src);

        for (int k = 0, dc = 0; k < rec->numRange(); k++, dc += ncomp)
        {
            rec->getRange(k, index, scomp, ncomp);
            FillPatch(*this,srcMF,ngrow_src,time,index,scomp,ncomp,dc);
        }

        mf = new MultiFab(dstBA, rec->numDerive(), ngrow);

#ifdef CRSEGRNDOMP
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
            const RealBox& temp = RealBox(gtbx,geom.CellSize(),geom.ProbLo());
            const Real* xlo     = temp.lo();
            Real        dt      = parent->dtLevel(level);

	    if (rec->derFunc() != static_cast<DeriveFunc>(0)){
		rec->derFunc()(ddat,ARLIM(dlo),ARLIM(dhi),&n_der,
			       cdat,ARLIM(clo),ARLIM(chi),&n_state,
			       lo,hi,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
			       &level,&grid_no);
	    } else if (rec->derFunc3D() != static_cast<DeriveFunc3D>(0)){
		rec->derFunc3D()(ddat,ARLIM_3D(dlo),ARLIM_3D(dhi),&n_der,
				 cdat,ARLIM_3D(clo),ARLIM_3D(chi),&n_state,
				 ARLIM_3D(lo),ARLIM_3D(hi),
				 ARLIM_3D(dom_lo),ARLIM_3D(dom_hi),
				 ZFILL(dx),ZFILL(xlo),
				 &time,&dt,
				 BCREC_3D(bcr),
				 &level,&grid_no);
	    } else {
		BoxLib::Error("AmeLevel::derive: no function available");
	    }
        }
#else
        for (MFIter mfi(srcMF); mfi.isValid(); ++mfi)
        {
            int         grid_no = mfi.index();
            RealBox     gridloc = RealBox(grids[grid_no],geom.CellSize(),geom.ProbLo());
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
		rec->derFunc()(ddat,ARLIM(dlo),ARLIM(dhi),&n_der,
			       cdat,ARLIM(clo),ARLIM(chi),&n_state,
			       dlo,dhi,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
			       &level,&grid_no);
	    } else if (rec->derFunc3D() != static_cast<DeriveFunc3D>(0)){
		rec->derFunc3D()(ddat,ARLIM_3D(dlo),ARLIM_3D(dhi),&n_der,
				 cdat,ARLIM_3D(clo),ARLIM_3D(chi),&n_state,
				 ARLIM_3D(dlo),ARLIM_3D(dhi),
				 ARLIM_3D(dom_lo),ARLIM_3D(dom_hi),
				 ZFILL(dx),ZFILL(xlo),
				 &time,&dt,
				 BCREC_3D(bcr),
				 &level,&grid_no);
	    } else {
		BoxLib::Error("AmeLevel::derive: no function available");
	    }
        }
#endif
    }
    else
    {
        //
        // If we got here, cannot derive given name.
        //
        std::string msg("AmrLevel::derive(MultiFab*): unknown variable: ");
        msg += name;
        BoxLib::Error(msg.c_str());
    }

    return mf;
}

void
AmrLevel::derive (const std::string& name,
                  Real           time,
                  MultiFab&      mf,
                  int            dcomp)
{
    BL_ASSERT(dcomp < mf.nComp());

    const int ngrow = mf.nGrow();

    int index, scomp, ncomp;

    if (isStateVariable(name,index,scomp))
    {
        FillPatch(*this,mf,ngrow,time,index,scomp,1);
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

        MultiFab srcMF(srcBA,rec->numState(),ngrow_src);

        for (int k = 0, dc = 0; k < rec->numRange(); k++, dc += ncomp)
        {
            rec->getRange(k,index,scomp,ncomp);

            FillPatch(*this,srcMF,ngrow_src,time,index,scomp,ncomp,dc);
        }

#ifdef CRSEGRNDOMP
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
		rec->derFunc()(ddat,ARLIM(dlo),ARLIM(dhi),&n_der,
			       cdat,ARLIM(clo),ARLIM(chi),&n_state,
			       lo,hi,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
			       &level,&idx);
	    } else if (rec->derFunc3D() != static_cast<DeriveFunc3D>(0)){
		rec->derFunc3D()(ddat,ARLIM_3D(dlo),ARLIM_3D(dhi),&n_der,
				 cdat,ARLIM_3D(clo),ARLIM_3D(chi),&n_state,
				 ARLIM_3D(lo),ARLIM_3D(hi),
				 ARLIM_3D(dom_lo),ARLIM_3D(dom_hi),
				 ZFILL(dx),ZFILL(xlo),
				 &time,&dt,
				 BCREC_3D(bcr),
				 &level,&idx);
	    } else {
		BoxLib::Error("AmeLevel::derive: no function available");
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
		rec->derFunc()(ddat,ARLIM(dlo),ARLIM(dhi),&n_der,
			       cdat,ARLIM(clo),ARLIM(chi),&n_state,
			       dlo,dhi,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
			       &level,&idx);
	    } else if (rec->derFunc3D() != static_cast<DeriveFunc3D>(0)){
		rec->derFunc3D()(ddat,ARLIM_3D(dlo),ARLIM_3D(dhi),&n_der,
				 cdat,ARLIM_3D(clo),ARLIM_3D(chi),&n_state,
				 ARLIM_3D(dlo),ARLIM_3D(dhi),
				 ARLIM_3D(dom_lo),ARLIM_3D(dom_hi),
				 ZFILL(dx),ZFILL(xlo),
				 &time,&dt,
				 BCREC_3D(bcr),
				 &level,&idx);
	    } else {
		BoxLib::Error("AmeLevel::derive: no function available");
	    }
        }
#endif
    }
    else
    {
        //
        // If we got here, cannot derive given name.
        //
        std::string msg("AmrLevel::derive(MultiFab*): unknown variable: ");
        msg += name;
        BoxLib::Error(msg.c_str());
    }
}

Array<int>
AmrLevel::getBCArray (int State_Type,
                      int gridno,
                      int strt_comp,
                      int ncomp)
{
    Array<int> bc(2*BL_SPACEDIM*ncomp);

    BCRec bcr;

    for (int n = 0; n < ncomp; n++)
    {
        bcr = state[State_Type].getBC(strt_comp+n,gridno);
        const int* b_rec = bcr.vect();
        for (int m = 0; m < 2*BL_SPACEDIM; m++)
            bc[2*BL_SPACEDIM*n + m] = b_rec[m];
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
  
}

AmrLevel::TimeLevel
AmrLevel::which_time (int  indx,
                      Real time) const
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

const BoxArray& AmrLevel::getAreaNotToTag()
{
    return m_AreaNotToTag;
}

const Box& AmrLevel::getAreaToTag()
{
    return m_AreaToTag;
}

void AmrLevel::setAreaNotToTag(BoxArray& ba)
{
    m_AreaNotToTag = ba;
}

void AmrLevel::constructAreaNotToTag()
{
    if (level == 0 || !Amr::useFixedCoarseGrids || Amr::useFixedUpToLevel>level)
        return;

    // We are restricting the tagging on the finest fixed level
    if (Amr::useFixedUpToLevel==level)
    {
        // We use the next coarser level shrunk by one blockingfactor
        //    as the region in which we allow tagging. 
        // Why level-1? Because we always use the full domain at level 0 
        //    and therefore level 0 in initialba is level 1 in the AMR hierarchy, etc.
        const Array<BoxArray>& initialba = parent->getInitialBA();
        Box tagarea(initialba[level-1].minimalBox());
        tagarea.grow(-parent->blockingFactor(level));
        m_AreaToTag = tagarea;

        // We disallow tagging in the remaining part of the domain.
        BoxArray tagba = BoxLib::boxComplement(parent->Geom(level).Domain(),m_AreaToTag);
        m_AreaNotToTag = tagba;

        BoxArray bxa(parent->Geom(level).Domain());
        BL_ASSERT(bxa.contains(m_AreaNotToTag));
    }

    if (Amr::useFixedUpToLevel<level)
    {
        Box tagarea = parent->getLevel(level-1).getAreaToTag();
        tagarea.refine(parent->refRatio(level-1));
        tagarea.grow(-parent->blockingFactor(level));
        m_AreaToTag = tagarea;
        BoxArray tagba = BoxLib::boxComplement(parent->Geom(level).Domain(),m_AreaToTag);
        m_AreaNotToTag = tagba;
    }
}

void
AmrLevel::FillPatch(AmrLevel& amrlevel,
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
AmrLevel::AddProcsToComp(Amr *aptr, int nSidecarProcs, int prevSidecarProcs,
                         int ioProcNumSCS, int ioProcNumAll, int scsMyId,
			 MPI_Comm scsComm)
{
#if BL_USE_MPI
      if(scsMyId != ioProcNumSCS) {
        parent = aptr;
      }

      // ---- ints
      ParallelDescriptor::Bcast(&level, 1, ioProcNumAll, scsComm);

      // ---- IntVects
      Array<int> allIntVects;
      if(scsMyId == ioProcNumSCS) {
        for(int i(0); i < BL_SPACEDIM; ++i)    { allIntVects.push_back(crse_ratio[i]); }
        for(int i(0); i < BL_SPACEDIM; ++i)    { allIntVects.push_back(fine_ratio[i]); }
      }
      BoxLib::BroadcastArray(allIntVects, scsMyId, ioProcNumSCS, scsComm);

      if(scsMyId != ioProcNumSCS) {
        int count(0);
        for(int i(0); i < BL_SPACEDIM; ++i)    { crse_ratio[i] = allIntVects[count++]; }
        for(int i(0); i < BL_SPACEDIM; ++i)    { fine_ratio[i] = allIntVects[count++]; }
      }


      // ---- Boxes
      BoxLib::BroadcastBox(m_AreaToTag, scsMyId, ioProcNumSCS, scsComm);
      
      // ---- Geometry
      Geometry::BroadcastGeometry(geom, ioProcNumSCS, scsComm);
      
      // ---- BoxArrays
      BoxLib::BroadcastBoxArray(grids, scsMyId, ioProcNumSCS, scsComm);
      BoxLib::BroadcastBoxArray(m_AreaNotToTag, scsMyId, ioProcNumSCS, scsComm);


#ifdef USE_PARTICLES
      BoxLib::BroadcastBoxArray(particle_grids, scsMyId, ioProcNumSCS, scsComm);

      int sentinelProc(ParallelDescriptor::MyProcComp());
      BoxLib::BroadcastDistributionMapping(particle_dmap, sentinelProc, scsMyId,
                                           ioProcNumSCS, scsComm, true);
#endif

#ifdef USE_SLABSTAT
      BoxLib::Abort("**** Error in AmrLevel::MSS:  USE_SLABSTAT not implemented");
#endif

      // ---- state
      int stateSize(state.size());
      ParallelDescriptor::Bcast(&stateSize, 1, ioProcNumSCS, scsComm);
      if(scsMyId != ioProcNumSCS) {
        state.resize(stateSize);
      }
      for(int i(0); i < state.size(); ++i) {
        state[i].AddProcsToComp(desc_lst[i], ioProcNumSCS, ioProcNumAll, scsMyId, scsComm);
      }
#endif
}



void
AmrLevel::Check() const
{
    for(int i(0); i < state.size(); ++i) {
      state[i].Check();
    }
}


