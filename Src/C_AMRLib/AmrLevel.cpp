//BL_COPYRIGHT_NOTICE

//
// $Id: AmrLevel.cpp,v 1.44 1998-10-07 21:18:19 vince Exp $
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
    m_ratio(m_amrlevel.level+1),
    m_map(m_amrlevel.level+1),
    m_finebox(m_leveldata.boxArray().length()),
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
    m_ratio(m_amrlevel.level+1),
    m_map(m_amrlevel.level+1),
    m_finebox(m_leveldata.boxArray().length()),
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

inline
void
Intersect (vector<Box>& boxes,
           const Box&   box)
{
    vector<Box> tmp;

    for (int i = 0; i < boxes.size(); i++)
    {
        if (boxes[i].intersects(box))
        {
            tmp.push_back(boxes[i] & box);
        }
    }

    tmp.swap(boxes);
}

inline
void
Coarsen (vector<Box>&   boxes,
         const IntVect& ratio)
{
    for (int i = 0; i < boxes.size(); i++)
    {
        boxes[i].coarsen(ratio);
    }
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
    //
    // This function sets up and performs the communication pattern for
    // filling fabs of size m_leveldata[i].box().grow(boxGrow) from amrlevel's
    // state data, possibly interpolating between the new and old data
    // the fill is done from this level and, if necessary, coarser levels.
    //
    m_time       = time;
    m_growsize   = boxGrow;
    m_stateindex = state_index;
    m_scomp      = src_comp;
    m_ncomp      = ncomp;

    Array<IntVect> pshifts(27);

    const int         MyProc    = ParallelDescriptor::MyProc();
    PArray<AmrLevel>& amrLevels = m_amrlevel.parent->getAmrLevels();

    m_ratio[m_amrlevel.level] = IntVect::TheUnitVector();
    for (int level = m_amrlevel.level - 1; level >= 0; --level)
    {
        m_ratio[level] = m_ratio[level+1] * amrLevels[level+1].crse_ratio;
    }
    for (int level = 0; level <= m_amrlevel.level; ++level)
    {
        StateData& state            = amrLevels[level].state[m_stateindex];
        const StateDescriptor& desc = amrLevels[level].desc_lst[m_stateindex];
        state.RegisterData(m_mfcd, m_mfid[level]);
        m_map[level] = mapper;
        if (m_map[level] == 0)
            m_map[level] = desc.interp();
    }
    for (int i = 0; i < m_ba.length(); ++i)
    {
        if (m_leveldata.DistributionMap()[i] == MyProc)
        {
            m_ba.set(i, m_leveldata.boxArray()[i]);
            m_fbid[i].resize(m_amrlevel.level + 1);
            m_finebox[i].resize(m_amrlevel.level + 1);
        }
    }
    m_ba.grow(m_growsize);  // These are the ones we want to fillpatch.

    const IndexType boxType(m_leveldata.boxArray()[0].ixType());

    vector<Box> unfilledThisLevel, tempUnfilledBoxes;

    BoxList unfillableThisLevel(boxType);
    //
    // Do this for all local (grown) fab boxes.
    //
    for (int ibox = 0; ibox < m_ba.length(); ++ibox)
    {
        if (m_leveldata.DistributionMap()[ibox] != MyProc)
            continue;

        unfilledThisLevel.clear();

        unfilledThisLevel.push_back(m_ba[ibox]);
        //
        // Find the boxes that can be filled on each level.
        // These are all defined at their level of refinement.
        //
        bool needsFilling = true;

        for (int level = m_amrlevel.level; level >= 0 && needsFilling; --level)
        {
            unfillableThisLevel.clear();

            StateData& currentState   = amrLevels[level].state[m_stateindex];
            const Box& currentPDomain = currentState.getDomain();
            bool is_periodic          = amrLevels[level].geom.isAnyPeriodic();

            if (is_periodic)
            {
                tempUnfilledBoxes.clear();

                for (int i = 0; i < unfilledThisLevel.size(); i++)
                {
                    assert(unfilledThisLevel[i].ok());

                    if (!currentPDomain.contains(unfilledThisLevel[i]))
                    {
                        amrLevels[level].geom.periodicShift(currentPDomain,
                                                            unfilledThisLevel[i],
                                                            pshifts);

                        for (int iiv = 0; iiv < pshifts.length(); iiv++)
                        {
                            Box shbox(unfilledThisLevel[i]);
                            shbox.shift(pshifts[iiv]);
                            shbox &= currentPDomain;
                            assert(shbox.ok());
                            tempUnfilledBoxes.push_back(shbox);
                        }
                    }
                }
                for (int i = 0; i < tempUnfilledBoxes.size(); i++)
                {
                    unfilledThisLevel.push_back(tempUnfilledBoxes[i]);
                }
            }
            m_fbid[ibox][level].resize(unfilledThisLevel.size());
            m_finebox[ibox][level].resize(unfilledThisLevel.size());

            int iBLI = 0;

            for (int i = 0; i < unfilledThisLevel.size(); i++)
            {
                assert(unfilledThisLevel[i].ok());

                if (unfilledThisLevel[i].intersects(currentPDomain))
                {
                    Box fineDestBox = unfilledThisLevel[i] & currentPDomain;

                    fineDestBox.refine(m_ratio[level]);

                    if (fineDestBox.intersects(m_ba[ibox]))
                    {
                        fineDestBox &= m_ba[ibox];
                    }
                    else
                    {
                        assert(is_periodic);
                    }
                    Box crse_box = fineDestBox;

                    if (level != m_amrlevel.level)
                    {
                        crse_box = m_map[level]->CoarseBox(fineDestBox,
                                                           m_ratio[level]);
                    }
                    m_finebox[ibox][level][iBLI] = fineDestBox;

                    BoxList tempUnfillableBoxes(boxType);

                    currentState.linInterpAddBox(m_mfcd,
                                                 m_mfid[level],
                                                 &tempUnfillableBoxes,
                                                 m_fbid[ibox][level][iBLI++],
                                                 crse_box,
                                                 m_time,
                                                 m_scomp,
                                                 0,
                                                 m_ncomp);

                    unfillableThisLevel.join(tempUnfillableBoxes);
                }
            }
            unfilledThisLevel.clear();
            unfillableThisLevel.intersect(currentPDomain);

            if (unfillableThisLevel.isEmpty())
            {
                needsFilling = false;
            }
            else
            {
                //
                // Populate `unfilledThisLevel' from `unfillableThisLevel'.
                //
                for (BoxListIterator bli(unfillableThisLevel); bli; ++bli)
                {
                    unfilledThisLevel.push_back(bli());
                }
                Box coarseLocalMFBox = ::coarsen(m_ba[ibox],m_ratio[level]);

                if (!is_periodic)
                {
                    Intersect(unfilledThisLevel,coarseLocalMFBox);
                }
                Coarsen(unfilledThisLevel,amrLevels[level].crse_ratio);

                if (level == 0)
                {
                    Intersect(unfilledThisLevel,currentPDomain);

                    if (!unfilledThisLevel.empty())
                    {
                        vector<Box> unfilledInside = unfilledThisLevel;
                        Intersect(unfilledInside,coarseLocalMFBox);
                        assert(unfilledInside.empty());
                    }
                }
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

    Array<BCRec>   bcr(m_ncomp);
    Array<IntVect> pshifts(27);

    const int              MyProc    = ParallelDescriptor::MyProc();
    PArray<AmrLevel>&      amrLevels = m_amrlevel.parent->getAmrLevels();
    const Box&             destBox   = m_ba[index()];
    const StateDescriptor& fDesc     = m_amrlevel.desc_lst[m_stateindex];
    StateData&             fState    = m_amrlevel.state[m_stateindex];
    const Box&             fDomain   = fState.getDomain();

    m_fab.resize(destBox, m_ncomp);
    m_fab.setVal(1.e30);

    FArrayBox tempCoarseDestFab, tempCurrentFillPatchedFab;

    for (int level = 0; level <= m_amrlevel.level; ++level)
    {
        const bool is_periodic  = amrLevels[level].geom.isAnyPeriodic();
        StateData& currentState = amrLevels[level].state[m_stateindex];

        for (int iBox = 0; iBox < m_fbid[currentIndex][level].length(); ++iBox)
        {
            const Box& tmpCrseBox = m_fbid[currentIndex][level][iBox][0].box();

            tempCoarseDestFab.resize(tmpCrseBox, m_ncomp);
            tempCoarseDestFab.setVal(1.e30);

            currentState.linInterpFillFab(m_mfcd,
                                          m_mfid[level],
                                          m_fbid[currentIndex][level][iBox],
                                          tempCoarseDestFab,
                                          m_time,
                                          0,
                                          0,
                                          m_ncomp);
            //
            // Set non-periodic BCs in tempCoarseDestFab.
            //
            if (!currentState.getDomain().contains(tempCoarseDestFab.box()))
            {
                currentState.FillBoundary(tempCoarseDestFab,
                                          m_time,
                                          amrLevels[level].geom.CellSize(),
                                          amrLevels[level].geom.ProbDomain(),
                                          0,
                                          0,
                                          m_ncomp);
            }
            Box iSectDest               = m_finebox[currentIndex][level][iBox];
            const BoxArray& filledBoxes = m_fbid[currentIndex][level][iBox][0].FilledBoxes();
            BoxArray fboxes             = filledBoxes;
            FArrayBox* cpFromFab        = 0;
            const BoxArray* cpFromBoxes = 0;

            if (iSectDest.ok())
            {
                if (level != m_amrlevel.level)
                {
                    //
                    // Get boundary conditions for this patch.
                    //
                    setBC(iSectDest,
                          fDomain,
                          m_scomp,
                          0,
                          m_ncomp,
                          fDesc.getBCs(),
                          bcr);

                    fboxes.refine(m_ratio[level]);
                    //
                    // Interpolate up to fine patch.
                    //
                    tempCurrentFillPatchedFab.resize(iSectDest,m_ncomp);

                    m_map[level]->interp(tempCoarseDestFab,
                                         0,
                                         tempCurrentFillPatchedFab,
                                         0,
                                         m_ncomp,
                                         iSectDest,
                                         m_ratio[level],
                                         amrLevels[level].geom,
                                         amrLevels[m_amrlevel.level].geom,
                                         bcr);

                    cpFromFab   = &tempCurrentFillPatchedFab;
                    cpFromBoxes = &fboxes;
                }
                else
                {
                    cpFromFab   = &tempCoarseDestFab;
                    cpFromBoxes = &filledBoxes;
                }

                for (int i = 0; i < cpFromBoxes->length(); ++i)
                {
                    Box srcdestBox = (*cpFromBoxes)[i];
                    srcdestBox    &= m_fab.box();
                    srcdestBox    &= iSectDest;

                    if (srcdestBox.ok())
                    {
                        m_fab.copy(*cpFromFab,
                                   srcdestBox,
                                   0,
                                   srcdestBox,
                                   0,
                                   m_ncomp);
                    }
                }

                if (is_periodic && !fDomain.contains(m_fab.box()))
                {
                    m_amrlevel.geom.periodicShift(fDomain,m_fab.box(),pshifts);

                    for (int iiv = 0; iiv < pshifts.length(); iiv++)
                    {
                        m_fab.shift(pshifts[iiv]);

                        for (int i = 0; i < cpFromBoxes->length(); ++i)
                        {
                            Box srcdestBox = (*cpFromBoxes)[i];
                            srcdestBox    &= m_fab.box();
                            srcdestBox    &= cpFromFab->box();

                            if (srcdestBox.ok())
                            {
                                m_fab.copy(*cpFromFab,
                                           srcdestBox,
                                           0,
                                           srcdestBox,
                                           0,
                                           m_ncomp);
                            }
                        }

                        m_fab.shift(-pshifts[iiv]);
                    }
                }
            }
        }
    }
    //
    // Do non-periodic BCs on the finest level.
    //
    if (!fDomain.contains(destBox))
    {
        fState.FillBoundary(m_fab,
                            m_time,
                            m_amrlevel.geom.CellSize(),
                            m_amrlevel.geom.ProbDomain(),
                            0,
                            m_scomp,
                            m_ncomp);
    }

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

    const StateDescriptor& desc = desc_lst[state_indx];
    Interpolater* map           = (mapper == 0) ? desc.interp() : mapper;
    const RealBox& prob_domain  = geom.ProbDomain();
    const Box& p_domain         = state[state_indx].getDomain();
    AmrLevel& crse_lev          = parent->getLevel(level-1);

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
