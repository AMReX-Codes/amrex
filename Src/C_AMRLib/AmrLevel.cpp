//BL_COPYRIGHT_NOTICE

//
// $Id: AmrLevel.cpp,v 1.49 1999-01-15 21:44:22 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdio>
#include <cstring>
#include <strstream>
#else
#include <stdio.h>
#include <string.h>
#include <strstream.h>
#endif

#include <AmrLevel.H>
#include <Derive.H>
#include <BoxDomain.H>
#include <ParallelDescriptor.H>
#include <RunStats.H>
#include <Utility.H>

DescriptorList AmrLevel::desc_lst;
DeriveList     AmrLevel::derive_lst;

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
{
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

    state.resize(desc_lst.length());

    for (int i = 0; i < state.length(); i++)
    {
        state[i].define(geom.Domain(),
                        grids,
                        desc_lst[i],
                        time,
                        parent->dtLevel(lev));
    }

    finishConstructor();
}

void
AmrLevel::restart (Amr&     papa,
                   istream& is,
		   bool bReadSpecial)
{
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

    if(bReadSpecial) {
      readBoxArray(grids, is, bReadSpecial);
    } else {
      grids.define(is);
    }

    int nstate;
    is >> nstate;
    int ndesc = desc_lst.length();
    assert(nstate == ndesc);

    state.resize(ndesc);
    for (int i = 0; i < ndesc; i++)
    {
        state[i].restart(is, desc_lst[i], papa.theRestartFile(), bReadSpecial);
    }

    finishConstructor();
}

void
AmrLevel::finishConstructor ()
{
    //
    // Set physical locations of grids.
    //
    grid_loc.resize(grids.length());

    for (int i = 0; i < grid_loc.length(); i++)
    {
        grid_loc[i] = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
    }
}

void
AmrLevel::setTimeLevel (Real time,
                        Real dt_old,
                        Real dt_new)
{
    for (int k = 0; k < desc_lst.length(); k++)
    {
        state[k].setTimeLevel(time,dt_old,dt_new);
    }
}

int
AmrLevel::isStateVariable (const aString& name,
                           int&           typ,
                           int&            n)
{
    for (typ = 0; typ < desc_lst.length(); typ++)
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
AmrLevel::countCells ()
{
    long cnt = 0;
    for (int i = 0; i < grids.length(); i++)
    {
        cnt += grids[i].numPts();
    }
    return cnt;
}

void
AmrLevel::checkPoint (const aString& dir,
                      ostream&       os,
                      VisMF::How     how)
{
    int ndesc = desc_lst.length(), i;
    //
    // Build directory to hold the MultiFabs in the StateData at this level.
    // The directory is relative the the directory containing the Header file.
    //
    char buf[64];
    sprintf(buf, "Level_%d", level);
    aString Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    aString FullPath = dir;
    if (!FullPath.isNull() && FullPath[FullPath.length()-1] != '/')
    {
        FullPath += '/';
    }
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!Utility::UtilCreateDirectory(FullPath, 0755))
            Utility::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

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
        aString PathNameInHeader = Level;
        sprintf(buf, "/SD_%d", i);
        PathNameInHeader += buf;
        aString FullPathName = FullPath;
        FullPathName += buf;
        state[i].checkPoint(PathNameInHeader, FullPathName, os, how);
    }
}

AmrLevel::~AmrLevel ()
{
    parent = 0;
}

void
AmrLevel::allocOldData ()
{
    for (int i = 0; i < desc_lst.length(); i++)
    {
        state[i].allocOldData();
    }
}

void
AmrLevel::removeOldData ()
{
    for (int i = 0; i < desc_lst.length(); i++)
    {
        state[i].removeOldData();
    }
}

void
AmrLevel::reset ()
{
    for (int i = 0; i < desc_lst.length(); i++)
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
    else
    {
        BoxLib::Error("get_data: invalid time");
        static MultiFab bogus;
        return bogus;
    }
}

void
AmrLevel::setPhysBoundaryValues (int  state_indx,
                                 int  comp,
                                 int  ncomp,
                                 Real time)
{
    const Real old_time = state[state_indx].prevTime();
    const Real new_time = state[state_indx].curTime();
    const Real eps = 0.001*(new_time - old_time);

    int do_new;
    if (time > old_time-eps && time < old_time+eps)
    {
        do_new = 0;
    }
    else if (time > new_time-eps && time < new_time+eps)
    {
        do_new = 1;
    }
    else
    {
        BoxLib::Error("setPhysBndryValues: invalid time");
    }

    state[state_indx].FillBoundary(geom.CellSize(),
                                   geom.ProbDomain(),
                                   comp,
                                   ncomp,
                                   do_new);
}

FillPatchIterator::FillPatchIterator (AmrLevel& amrlevel,
                                      MultiFab& leveldata)
    :
    MultiFabIterator(leveldata),
    m_amrlevel(amrlevel),
    m_leveldata(leveldata),
    m_mfid(m_amrlevel.level+1),
    m_map(m_amrlevel.level+1),
    m_cfab(m_amrlevel.level+1),
    m_finebox(m_leveldata.boxArray().length()),
    m_crsebox(m_leveldata.boxArray().length()),
    m_fbid(m_leveldata.boxArray().length()),
    m_ba(m_leveldata.boxArray().length()),
    m_init(false)
{}

FillPatchIterator::FillPatchIterator (AmrLevel&     amrlevel,
                                      MultiFab&     leveldata,
                                      int           boxGrow,
                                      Real          time,
                                      int           state_index,
                                      int           src_comp,
                                      int           ncomp,
                                      Interpolater* mapper)
    :
    MultiFabIterator(leveldata),
    m_amrlevel(amrlevel),
    m_leveldata(leveldata),
    m_mfid(m_amrlevel.level+1),
    m_map(m_amrlevel.level+1),
    m_cfab(m_amrlevel.level+1),
    m_finebox(m_leveldata.boxArray().length()),
    m_crsebox(m_leveldata.boxArray().length()),
    m_fbid(m_leveldata.boxArray().length()),
    m_ba(m_leveldata.boxArray().length()),
    m_time(time),
    m_growsize(boxGrow),
    m_stateindex(state_index),
    m_scomp(src_comp),
    m_ncomp(ncomp),
    m_init(false)
{
    Initialize(boxGrow,time,state_index,src_comp,ncomp,mapper);
}
 
//
// Used in a couple RunStat calls in FillPatchIterator.
//
static const aString RunstatString("fill_patch");

void
FillPatchIterator::Initialize (int           boxGrow,
                               Real          time,
                               int           state_index,
                               int           src_comp,
                               int           ncomp,
                               Interpolater* mapper)
{
    RunStats stats(RunstatString, m_amrlevel.Level());

    stats.start();

    m_time       = time;
    m_growsize   = boxGrow;
    m_stateindex = state_index;
    m_scomp      = src_comp;
    m_ncomp      = ncomp;

    const int         MyProc     = ParallelDescriptor::MyProc();
    PArray<AmrLevel>& amrLevels  = m_amrlevel.parent->getAmrLevels();
    const AmrLevel&   topLevel   = amrLevels[m_amrlevel.level];
    const Box&        topPDomain = topLevel.state[m_stateindex].getDomain();
    const IndexType   boxType    = m_leveldata.boxArray()[0].ixType();

    for (int l = 0; l <= m_amrlevel.level; ++l)
    {
        amrLevels[l].state[m_stateindex].RegisterData(m_mfcd, m_mfid[l]);

        if ((m_map[l] = mapper) == 0)
            m_map[l] = amrLevels[l].desc_lst[m_stateindex].interp();
    }
    for (int i = 0; i < m_ba.length(); ++i)
    {
        if (m_leveldata.DistributionMap()[i] == MyProc)
        {
            m_ba.set(i, m_leveldata.boxArray()[i]);
            m_fbid[i].resize(m_amrlevel.level + 1);
            m_finebox[i].resize(m_amrlevel.level + 1);
            m_crsebox[i].resize(m_amrlevel.level + 1);
        }
    }
    m_ba.grow(m_growsize);  // These are the ones we want to fillpatch.

    BoxList unfillableThisLevel(boxType), tempUnfillable(boxType);
    vector<Box> unfilledThisLevel, crse_boxes;

    if (topLevel.geom.isAnyPeriodic())
        m_pshifts.resize(27);

    for (int ibox = 0; ibox < m_ba.length(); ++ibox)
    {
        if (m_leveldata.DistributionMap()[ibox] != MyProc)
            continue;

        unfilledThisLevel.clear();
        unfilledThisLevel.push_back(m_ba[ibox]);

        if (!topPDomain.contains(m_ba[ibox]))
        {
            unfilledThisLevel.back() &= topPDomain;

            if (topLevel.geom.isAnyPeriodic())
            {
                //
                // May need to add additional unique pieces of valid region
                // in order to do periodic copies into ghost cells.
                //
                topLevel.geom.periodicShift(topPDomain,m_ba[ibox],m_pshifts);

                for (int iiv = 0; iiv < m_pshifts.length(); iiv++)
                {
                    Box shbox  = m_ba[ibox] + m_pshifts[iiv];
                    shbox     &= topPDomain;

                    if (boxType.nodeCentered())
                    {
                        for (int dir = 0; dir < BL_SPACEDIM; dir++)
                        {
                            if (m_pshifts[iiv][dir] > 0)
                                shbox.growHi(dir,-1);
                            else if (m_pshifts[iiv][dir] < 0)
                                shbox.growLo(dir,-1);
                        }
                    }

                    if (shbox.ok())
                    {
                        BoxList bl = ::boxDiff(shbox,m_ba[ibox]);
                        for (BoxListIterator bli(bl); bli; ++bli)
                            unfilledThisLevel.push_back(bli());
                    }
                }
            }
        }

        bool Done = false;

        for (int l = m_amrlevel.level; l >= 0 && !Done; --l)
        {
            unfillableThisLevel.clear();

            StateData&      theState      = amrLevels[l].state[m_stateindex];
            const Box&      thePDomain    = theState.getDomain();
            const Geometry& theGeom       = amrLevels[l].geom;
            const bool      is_periodic   = theGeom.isAnyPeriodic();
            const IntVect&  fine_ratio    = amrLevels[l].fine_ratio;
            //
            // These are the boxes on this level contained in thePDomain
            // that need to be filled in order to directly fill at the
            // highest level or to interpolate up to the next higher level.
            //
            m_finebox[ibox][l].resize(unfilledThisLevel.size());

            for (int i = 0; i < unfilledThisLevel.size(); i++)
                m_finebox[ibox][l][i] = unfilledThisLevel[i];
            //
            // Now build coarse boxes needed to interpolate to fine.
            //
            // If we're periodic and we're not at the finest level, we may
            // need to get some additional data at this level in order to
            // properly fill the CoarseBox()d versions of the fineboxes.
            //
            crse_boxes.clear();

            const Array<Box>& FineBoxes = m_finebox[ibox][l];

            for (int i = 0; i < FineBoxes.length(); i++)
            {
                crse_boxes.push_back(FineBoxes[i]);

                if (l != m_amrlevel.level)
                {
                    Box& cbox = crse_boxes.back();

                    cbox = m_map[l]->CoarseBox(FineBoxes[i],fine_ratio);

                    if (is_periodic && !thePDomain.contains(cbox))
                    {
                        theGeom.periodicShift(thePDomain,cbox,m_pshifts);

                        for (int iiv = 0; iiv < m_pshifts.length(); iiv++)
                        {
                            Box shbox = cbox + m_pshifts[iiv];
                            shbox    &= thePDomain;

                            if (boxType.nodeCentered())
                            {
                                for (int dir = 0; dir < BL_SPACEDIM; dir++)
                                {
                                    if (m_pshifts[iiv][dir] > 0)
                                        shbox.growHi(dir,-1);
                                    else if (m_pshifts[iiv][dir] < 0)
                                        shbox.growLo(dir,-1);
                                }
                            }

                            if (shbox.ok())
                                crse_boxes.push_back(shbox);
                        }
                    }
                }
            }

            m_crsebox[ibox][l].resize(crse_boxes.size());

            m_fbid[ibox][l].resize(crse_boxes.size());
            //
            // Now attempt to get as much coarse data as possible.
            //
            Array<Box>& CrseBoxes = m_crsebox[ibox][l];

            for (int i = 0; i < CrseBoxes.length(); i++)
            {
                assert(tempUnfillable.isEmpty());

                CrseBoxes[i] = crse_boxes[i];

                assert(CrseBoxes[i].intersects(thePDomain));

                theState.linInterpAddBox(m_mfcd,
                                         m_mfid[l],
                                         &tempUnfillable,
                                         m_fbid[ibox][l][i],
                                         CrseBoxes[i],
                                         m_time,
                                         m_scomp,
                                         0,
                                         m_ncomp);

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

                for (BoxListIterator bli(unfillableThisLevel); bli; ++bli)
                    unfilledThisLevel.push_back(bli());
            }
        }
    }

    m_mfcd.CollectData();

    m_init = true;

    stats.end();
}

bool
FillPatchIterator::isValid ()
{
    assert(m_init);

    if (!MultiFabIterator::isValid())
        return false;

    RunStats stats(RunstatString, m_amrlevel.Level());

    stats.start();

    Array<BCRec>      bcr(m_ncomp);
    PArray<AmrLevel>& amrLevels  = m_amrlevel.parent->getAmrLevels();
    //
    // The ultimate destination.
    //
    m_fab.resize(m_ba[index()], m_ncomp);
    //
    // Set to special value we'll later check to ensure we've filled the FAB.
    //
    m_fab.setVal(2.e30);
    //
    // Build all coarse fabs from which we'll interpolate and
    // fill them with coarse data as best we can.
    //
    for (int l = 0; l <= m_amrlevel.level; l++)
    {
        StateData&         TheState = amrLevels[l].state[m_stateindex];
        PArray<FArrayBox>& CrseFabs = m_cfab[l];

        m_cfab[l].resize(m_crsebox[currentIndex][l].length(),PArrayManage);

        for (int i = 0; i < CrseFabs.length(); i++)
        {
            const Box& cbox = m_crsebox[currentIndex][l][i];

            assert(cbox.ok());
            //
            // Set to special value we'll later check
            // to ensure we've filled the FABs at the coarse level.
            //
            CrseFabs.set(i, new FArrayBox(cbox,m_ncomp));

            CrseFabs[i].setVal(3.e30);

            TheState.linInterpFillFab(m_mfcd,
                                      m_mfid[l],
                                      m_fbid[currentIndex][l][i],
                                      CrseFabs[i],
                                      m_time,
                                      0,
                                      0,
                                      m_ncomp);
        }
    }
    //
    // Now work from the bottom up interpolating to next higher level.
    //
    FArrayBox fine_fab, crse_fab;

    for (int l = 0; l < m_amrlevel.level; l++)
    {
        StateData&         TheState = amrLevels[l].state[m_stateindex];
        PArray<FArrayBox>& CrseFabs = m_cfab[l];
        const Geometry&    TheGeom  = amrLevels[l].geom;
        //
        // Set non-periodic BCs in coarse data -- what we interpolate with.
        //
        for (int i = 0; i < CrseFabs.length(); i++)
        {
            if (!TheState.getDomain().contains(CrseFabs[i].box()))
            {
                TheState.FillBoundary(CrseFabs[i],
                                      m_time,
                                      TheGeom.CellSize(),
                                      TheGeom.ProbDomain(),
                                      0,
                                      m_scomp,
                                      m_ncomp);
            }
        }

        if (TheGeom.isAnyPeriodic())
        {
            //
            // Fill CrseFabs with periodic data in preparation for interp().
            //
            const Box& thePDomain  = TheState.getDomain();

            for (int i = 0; i < CrseFabs.length(); i++)
            {
                FArrayBox& dstfab = CrseFabs[i];

                if (!thePDomain.contains(dstfab.box()))
                {
                    TheGeom.periodicShift(thePDomain,dstfab.box(),m_pshifts);

                    for (int iiv = 0; iiv < m_pshifts.length(); iiv++)
                    {
                        Box fullsrcbox = dstfab.box() + m_pshifts[iiv];
                        fullsrcbox    &= thePDomain;

                        for (int j = 0; j < CrseFabs.length(); j++)
                        {
                            FArrayBox& srcfab = CrseFabs[j];

                            if (fullsrcbox.intersects(srcfab.box()))
                            {
                                Box srcbox = fullsrcbox & srcfab.box();
                                Box dstbox = srcbox - m_pshifts[iiv];

                                dstfab.copy(srcfab,
                                            srcbox,
                                            0,
                                            dstbox,
                                            0,
                                            m_ncomp);
                            }
                        }
                    }
                }
            }
        }
        //
        // Interpolate up to next level.
        //
        const IntVect&     fine_ratio    = amrLevels[l].fine_ratio;
        const Array<Box>&  FineBoxes     = m_finebox[currentIndex][l];
        const StateDescriptor& fDesc     = amrLevels[l+1].desc_lst[m_stateindex];
        StateData&         fState        = amrLevels[l+1].state[m_stateindex];
        const Box&         fDomain       = fState.getDomain();
        PArray<FArrayBox>& FinerCrseFabs = m_cfab[l+1];

        for (int i = 0; i < FineBoxes.length(); i++)
        {
            fine_fab.resize(FineBoxes[i],m_ncomp);

            Box crse_box = m_map[l]->CoarseBox(fine_fab.box(),fine_ratio);

            crse_fab.resize(crse_box,m_ncomp);
            //
            // Fill crse_fab from m_crsebox via copy on intersect.
            //
            for (int j = 0; j < CrseFabs.length(); j++)
            {
                crse_fab.copy(CrseFabs[j]);
            }
            //
            // Get boundary conditions for the fine patch.
            //
            setBC(fine_fab.box(),
                  fDomain,
                  m_scomp,
                  0,
                  m_ncomp,
                  fDesc.getBCs(),
                  bcr);
            //
            // The coarse FAB had better be completely filled with "good" data.
            //
            assert(crse_fab.norm(0) < 3.e30);
            //
            // Interpolate up to fine patch.
            //
            m_map[l]->interp(crse_fab,
                             0,
                             fine_fab,
                             0,
                             m_ncomp,
                             fine_fab.box(),
                             fine_ratio,
                             amrLevels[l].geom,
                             amrLevels[l+1].geom,
                             bcr);
            //
            // Copy intersect fine_fab into next level m_crseboxes.
            //
            for (int j = 0; j < FinerCrseFabs.length(); j++)
            {
                FinerCrseFabs[j].copy(fine_fab);
            }
        }
        //
        // No longer need coarse data at this level.
        //
        CrseFabs.clear();
    }
    //
    // Now for the finest level stuff.
    //
    StateData&         FineState      = m_amrlevel.state[m_stateindex];
    const Box&         FineDomain     = FineState.getDomain();
    PArray<FArrayBox>& FinestCrseFabs = m_cfab[m_amrlevel.level];
    //
    // Copy intersect coarse into destination fab.
    //
    for (int i = 0; i < FinestCrseFabs.length(); i++)
    {
        assert(FineState.getDomain().contains(FinestCrseFabs[i].box()));

        m_fab.copy(FinestCrseFabs[i]);
    }

    if (m_amrlevel.geom.isAnyPeriodic() && !FineDomain.contains(m_fab.box()))
    {
        m_amrlevel.geom.periodicShift(FineDomain,m_fab.box(),m_pshifts);

        for (int iiv = 0; iiv < m_pshifts.length(); iiv++)
        {
            m_fab.shift(m_pshifts[iiv]);

            for (int i = 0; i < FinestCrseFabs.length(); i++)
            {
                Box src_dst_box = FinestCrseFabs[i].box() & m_fab.box();
                src_dst_box    &= FineDomain;

                if (src_dst_box.ok())
                {
                    m_fab.copy(FinestCrseFabs[i],
                               src_dst_box,
                               0,
                               src_dst_box,
                               0,
                               m_ncomp);
                }
            }

            m_fab.shift(-m_pshifts[iiv]);
        }
    }
    //
    // No longer need coarse data at finest level.
    //
    FinestCrseFabs.clear();
    //
    // Final set of non-periodic BCs.
    //
    if (!FineState.getDomain().contains(m_fab.box()))
    {
        FineState.FillBoundary(m_fab,
                               m_time,
                               m_amrlevel.geom.CellSize(),
                               m_amrlevel.geom.ProbDomain(),
                               0,
                               m_scomp,
                               m_ncomp);
    }
    //
    // The final FAB had better be completely filled with "good" data.
    //
    assert(m_fab.norm(0) < 2.e30);

    stats.end();

    return true;
}

FillPatchIterator::~FillPatchIterator () {}

void
AmrLevel::FillCoarsePatch (MultiFab&     mfdest,
                           int           dcomp,
                           Real          time,
                           int           state_indx,
                           int           scomp,
                           int           ncomp,
                           Interpolater* mapper)
{
    //
    // Must fill this region on crse level and interpolate.
    //
    assert(level != 0);
    assert(ncomp <= (mfdest.nComp()-dcomp));
    assert((0 <= state_indx) && (state_indx < desc_lst.length()));

    const StateDescriptor& desc        = desc_lst[state_indx];
    Interpolater*          map         = mapper == 0 ? desc.interp() : mapper;
    const RealBox&         prob_domain = geom.ProbDomain();
    const Box&             p_domain    = state[state_indx].getDomain();
    AmrLevel&              crse_lev    = parent->getLevel(level-1);

    assert(desc.inRange(scomp, ncomp));
    //
    // Build a properly coarsened boxarray.
    //
    BoxArray crse_regBoxArray(mfdest.boxArray().length());

    for (int ibox = 0; ibox < crse_regBoxArray.length(); ++ibox)
    {
        Box dbox = mfdest.fabbox(ibox) & p_domain;

        assert(dbox.ixType() == desc.getType());
        //
        // Coarsen unfilled region and widen by interpolater stencil width.
        //
        crse_regBoxArray.set(ibox,map->CoarseBox(dbox,crse_ratio));
    }

    Array<BCRec> bcr(ncomp);

    MultiFab mf(crse_regBoxArray, ncomp, 0, Fab_noallocate);

    FillPatchIterator fpi(crse_lev,mf,0,time,state_indx,scomp,ncomp,mapper);

    for ( ; fpi.isValid(); ++fpi)
    {
        DependentMultiFabIterator mfdest_mfi(fpi, mfdest);

        assert(mfdest_mfi.fabbox() == mfdest_mfi().box());

        Box dbox = mfdest_mfi().box() & p_domain;
        //
        // Get bndry conditions for this patch.
        //
        setBC(dbox,p_domain,scomp,0,ncomp,desc.getBCs(),bcr);
        //
        // Interpolate up to fine patch.
        //
        map->interp(fpi(),
                    0,
                    mfdest_mfi(),
                    dcomp,
                    ncomp,
                    dbox,
                    crse_ratio,
                    crse_lev.geom,geom,bcr);

        if (!p_domain.contains(mfdest_mfi().box()))
        {
            state[state_indx].FillBoundary(mfdest_mfi(),
                                           time,
                                           geom.CellSize(),
                                           prob_domain,
                                           dcomp,
                                           scomp,
                                           ncomp);
        }
    }
}

MultiFab*
AmrLevel::derive (const aString& name,
                  Real           time,
                  int            ngrow)
{
    assert(ngrow >= 0);

    MultiFab* mf = 0;

    int state_indx, src_comp, num_comp;

    if (isStateVariable(name, state_indx, src_comp))
    {
        mf = new MultiFab(state[state_indx].boxArray(), 1, ngrow);

        FillPatchIterator fpi(*this, get_new_data(state_indx), ngrow,
                              time, state_indx, src_comp, 1);

        for ( ; fpi.isValid(); ++fpi)
        {
            assert((*mf)[fpi.index()].box() == fpi().box());

            (*mf)[fpi.index()].copy(fpi());
        }
    }
    else if (const DeriveRec* rec = derive_lst.get(name))
    {
        rec->getRange(0, state_indx, src_comp, num_comp);

        BoxArray srcBA(state[state_indx].boxArray());
        BoxArray dstBA(state[state_indx].boxArray());

        srcBA.convert(rec->boxMap());
        dstBA.convert(rec->deriveType());

        MultiFab srcMF(srcBA, rec->numState(), ngrow);

        for (int k = 0, dc = 0; k < rec->numRange(); k++, dc += num_comp)
        {
            rec->getRange(k, state_indx, src_comp, num_comp);

            FillPatchIterator fpi(*this, srcMF, ngrow, time, state_indx,
                                  src_comp, num_comp, rec->interp());

            for ( ; fpi.isValid(); ++fpi)
            {
                DependentMultiFabIterator dmfi(fpi, srcMF);
                dmfi().copy(fpi(), 0, dc, num_comp);
            }
        }

        mf = new MultiFab(dstBA, rec->numDerive(), ngrow);

        for (MultiFabIterator mfi(srcMF); mfi.isValid(); ++mfi)
        {
            int grid_no       = mfi.index();
            Real* ddat        = (*mf)[grid_no].dataPtr();
            const int* dlo    = (*mf)[grid_no].loVect();
            const int* dhi    = (*mf)[grid_no].hiVect();
            int n_der         = rec->numDerive();
            Real* cdat        = mfi().dataPtr();
            const int* clo    = mfi().loVect();
            const int* chi    = mfi().hiVect();
            int n_state       = rec->numState();
            const int* dom_lo = state[state_indx].getDomain().loVect();
            const int* dom_hi = state[state_indx].getDomain().hiVect();
            const Real* dx    = geom.CellSize();
            const int* bcr    = rec->getBC();
            const Real* xlo   = grid_loc[grid_no].lo();
            Real dt           = parent->dtLevel(level);

            rec->derFunc()(ddat,ARLIM(dlo),ARLIM(dhi),&n_der,
                           cdat,ARLIM(clo),ARLIM(chi),&n_state,
                           dlo,dhi,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
                           &level,&grid_no);
        }
    }
    else
    {
        //
        // If we got here, cannot derive given name.
        //
        aString msg("AmrLevel::derive(MultiFab*): unknown variable: ");
        msg += name;
        BoxLib::Error(msg.c_str());
    }

    return mf;
}

Array<int>
AmrLevel::getBCArray (int State_Type,
                      int gridno,
                      int strt_comp,
                      int num_comp)
{
    Array<int> bc(2*BL_SPACEDIM*num_comp);

    for (int n = 0; n < num_comp; n++)
    {
        const int* b_rec = state[State_Type].getBC(strt_comp+n,gridno).vect();
        for (int m = 0; m < 2*BL_SPACEDIM; m++)
            bc[2*BL_SPACEDIM*n + m] = b_rec[m];
    }

    return bc;
}
