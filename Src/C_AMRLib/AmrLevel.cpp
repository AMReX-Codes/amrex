//BL_COPYRIGHT_NOTICE

//
// $Id: AmrLevel.cpp,v 1.27 1998-04-01 18:21:51 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <strstream>
#else
#include <strstream.h>
#endif

#ifdef BL_USE_NEW_HFILES
#include <cstdio>
#include <cstring>
#else
#include <stdio.h>
#include <string.h>
#endif

#include <AmrLevel.H>
#include <Derive.H>
#include <BoxDomain.H>
#include <ParallelDescriptor.H>
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
    level = lev;
    parent = &papa;

    fine_ratio = IntVect::TheUnitVector(); fine_ratio.scale(-1);
    crse_ratio = IntVect::TheUnitVector(); crse_ratio.scale(-1);

    if (level > 0)
        crse_ratio = parent->refRatio(level-1);
    if (level < parent->maxLevel())
        fine_ratio = parent->refRatio(level);

    Real dt = parent->dtLevel(lev);
    int ndesc = desc_lst.length();
    state.resize(ndesc);
    const Box& domain = geom.Domain();
    for (int i = 0; i < ndesc; i++)
        state[i].define(domain,grids,desc_lst[i],time, dt);

    finishConstructor();
}

void
AmrLevel::restart (Amr&     papa,
                   istream& is)
{
    parent = &papa;

    is >> level;
    is >> geom;

    fine_ratio = IntVect::TheUnitVector(); fine_ratio.scale(-1);
    crse_ratio = IntVect::TheUnitVector(); crse_ratio.scale(-1);

    if (level > 0)
        crse_ratio = parent->refRatio(level-1);
    if (level < parent->maxLevel())
        fine_ratio = parent->refRatio(level);

    grids.define(is);

    int nstate;
    is >> nstate;
    int ndesc = desc_lst.length();
    assert(nstate == ndesc);

    state.resize(ndesc);
    for (int i = 0; i < ndesc; i++)
    {
        state[i].restart(is, desc_lst[i], papa.theRestartFile());
    }

    finishConstructor();
}

void
AmrLevel::finishConstructor()
{
    //
    // Set physical locations of grids.
    //
    int num_grids = grids.length();
    grid_loc.resize(num_grids);
    const Real* prob_lo = geom.ProbLo();
    const Real* dx = geom.CellSize();
    for (int i = 0; i < num_grids; i++)
        grid_loc[i] = RealBox(grids[i],dx,prob_lo);
}

void
AmrLevel::setTimeLevel (Real time,
                        Real dt_old,
                        Real dt_new)
{
    int ndesc = desc_lst.length();
    for (int k = 0; k < ndesc; k++)
        state[k].setTimeLevel(time,dt_old,dt_new);
}

int
AmrLevel::isStateVariable (const aString& name,
                           int&           typ,
                           int&            n)
{
    int ndesc = desc_lst.length();
    for (typ = 0; typ < ndesc; typ++)
    {
        const StateDescriptor &desc = desc_lst[typ];
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
        cnt += grids[i].numPts();
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
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!Utility::CreateDirectory(FullPath, 0755))
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
    int ndesc = desc_lst.length();
    for (int i = 0; i < ndesc; i++)
        state[i].allocOldData();
}

void
AmrLevel::removeOldData()
{
    int ndesc = desc_lst.length();
    for (int i = 0; i < ndesc; i++)
        state[i].removeOldData();
}

void
AmrLevel::reset ()
{
    int ndesc = desc_lst.length();
    for (int i = 0; i < ndesc; i++)
        state[i].reset();
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
        BoxLib::Error("setPhysBndryValues: invalid time");

    state[state_indx].FillBoundary(geom.CellSize(),geom.ProbDomain(),comp,
                                   ncomp, do_new);
}

FillPatchIterator::FillPatchIterator (AmrLevel& amrlevel,
                                      MultiFab& leveldata)
    :
    MultiFabIterator(leveldata),
    amrLevel(amrlevel),
    levelData(leveldata),
    multiFabCopyDesc(true),
    bIsInitialized(false)
{}

#ifndef BL_NEWFPMINBOX
FillPatchIterator::FillPatchIterator (AmrLevel&     amrlevel,
                                      MultiFab&     leveldata,
                                      int           boxGrow,
                                      int           dest_comp,
                                      Real          time,
                                      int           state_index,
                                      int           src_comp,
                                      int           ncomp,
                                      Interpolater* mapper)
    :
    MultiFabIterator(leveldata),
    amrLevel(amrlevel),
    levelData(leveldata),
    growSize(boxGrow),
    stateIndex(state_index),
    srcComp(src_comp),
    destComp(dest_comp),
    nComp(ncomp),
    interpTime(time),
    multiFabCopyDesc(true),
    bIsInitialized(false)
{
    Initialize(boxGrow, dest_comp, time, state_index, src_comp, ncomp, mapper);
}

void
FillPatchIterator::Initialize (int           boxGrow,
                               int           dest_comp,
                               Real          time,
                               int           state_index,
                               int           src_comp,
                               int           ncomp,
                               Interpolater* mapper)
{
    //
    // This function sets up and performs the communication pattern for
    // filling fabs of size levelData[i].box().grow(boxGrow) from amrlevel's
    // state data, possibly interpolating between the new and old data
    // the fill is done from this level and, if necessary, coarser levels.
    //
    growSize   = boxGrow;
    stateIndex = state_index;
    srcComp    = src_comp;
    destComp   = dest_comp;
    nComp      = ncomp;
    interpTime = time;

    int myproc = ParallelDescriptor::MyProc();
    int currentLevel;
    PArray<AmrLevel> &amrLevels = amrLevel.parent->getAmrLevels();
    cumulativeRefRatios.resize(amrLevel.level + 1);
    map.resize(amrLevel.level + 1);

    cumulativeRefRatios[amrLevel.level] = IntVect::TheUnitVector();
    for (currentLevel = amrLevel.level - 1; currentLevel >= 0; --currentLevel)
    {
        cumulativeRefRatios[currentLevel] = cumulativeRefRatios[currentLevel + 1] *
            amrLevels[currentLevel + 1].crse_ratio;
    }

    stateDataMFId.resize(amrLevel.level + 1);
    for (currentLevel = 0; currentLevel <= amrLevel.level; ++currentLevel)
    {
        StateData& currentState = amrLevels[currentLevel].state[stateIndex];
        const StateDescriptor& currentDesc =
            amrLevels[currentLevel].desc_lst[stateIndex];
        currentState.RegisterData(multiFabCopyDesc, stateDataMFId[currentLevel]);
        map[currentLevel] = mapper;
        if (map[currentLevel] == 0)
            map[currentLevel] = currentDesc.interp();
    }

    localMFBoxes.resize(levelData.boxArray().length());
    fillBoxId.resize(levelData.boxArray().length());
    savedFineBox.resize(levelData.boxArray().length());
    for (int iLocal = 0; iLocal < localMFBoxes.length(); ++iLocal)
    {
        if (levelData.DistributionMap().ProcessorMap()[iLocal] == myproc)
        {
            //
            // Local.
            //
            localMFBoxes.set(iLocal, levelData.boxArray()[iLocal]);
            fillBoxId[iLocal].resize(amrLevel.level + 1);
            savedFineBox[iLocal].resize(amrLevel.level + 1);
        }
    }
    localMFBoxes.grow(growSize);  // These are the ones we want to fillpatch.

    IndexType boxType(levelData.boxArray()[0].ixType());
    BoxList unfilledBoxesOnThisLevel(boxType);
    BoxList unfillableBoxesOnThisLevel(boxType);
    //
    // Do this for all local (grown) fab boxes.
    //
    for (int ibox = 0; ibox < localMFBoxes.length(); ++ibox)
    {
        if (levelData.DistributionMap().ProcessorMap()[ibox] != myproc)
            continue;  // Not local.

        unfilledBoxesOnThisLevel.clear();
        assert(unfilledBoxesOnThisLevel.ixType() == boxType);
        assert(unfilledBoxesOnThisLevel.ixType() == localMFBoxes[ibox].ixType());
        unfilledBoxesOnThisLevel.add(localMFBoxes[ibox]);
        //
        // Find the boxes that can be filled on each level--these are all
        // defined at their level of refinement.
        //
        bool needsFilling = true;
        for (currentLevel = amrLevel.level; currentLevel >= 0 && needsFilling;
            --currentLevel)
        {
            unfillableBoxesOnThisLevel.clear();

            StateData& currentState   = amrLevels[currentLevel].state[stateIndex];
            const Box& currentPDomain = currentState.getDomain();

// vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv periodic
            bool is_periodic = amrLevels[currentLevel].geom.isAnyPeriodic();

            if (is_periodic)
            {
                BoxList tempUnfilledBoxes(boxType);
                for (BoxListIterator pbli(unfilledBoxesOnThisLevel); pbli; ++pbli) {
                    assert(pbli().ok());
                    const Box& dbox = pbli();
                    if (!currentPDomain.contains(dbox))
                    {
                        Array<IntVect> pshifts(27);
                        amrLevels[currentLevel].geom.periodicShift(currentPDomain,dbox,pshifts);
                        for (int iiv = 0; iiv < pshifts.length(); iiv++)
                        {
                            const IntVect& iv = pshifts[iiv];
                            Box shbox(dbox);
                            D_TERM(shbox.shift(0,iv[0]);,
                                   shbox.shift(1,iv[1]);,
                                   shbox.shift(2,iv[2]);)
                            shbox &= currentPDomain;
                            if (shbox.ok())
                                tempUnfilledBoxes.add(shbox);
                        }
                    }
                }
                unfilledBoxesOnThisLevel.join(tempUnfilledBoxes);
            }
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ end periodic

            fillBoxId[ibox][currentLevel].resize(unfilledBoxesOnThisLevel.length());
            savedFineBox[ibox][currentLevel].resize(unfilledBoxesOnThisLevel.length());

            int currentBLI = 0;
            for (BoxListIterator bli(unfilledBoxesOnThisLevel); bli; ++bli)
            {
                assert(bli().ok());
                Box coarseDestBox(bli());
                Box fineTruncDestBox(coarseDestBox & currentPDomain);
                if (fineTruncDestBox.ok())
                {
                    fineTruncDestBox.refine(cumulativeRefRatios[currentLevel]);
                    Box tempCoarseBox;
                    if (currentLevel == amrLevel.level)
                    {
                        tempCoarseBox = fineTruncDestBox;
                    }
                    else
                    {
                        tempCoarseBox = map[currentLevel]->CoarseBox(fineTruncDestBox,
                                                                     cumulativeRefRatios[currentLevel]);
                    }

                    savedFineBox[ibox][currentLevel][currentBLI] = fineTruncDestBox;
                    if (!is_periodic)
                    {
                        assert(localMFBoxes[ibox].intersects(fineTruncDestBox));
                    }

                    BoxList tempUnfillableBoxes(boxType);
                    currentState.linInterpAddBox(multiFabCopyDesc,
                                                 stateDataMFId[currentLevel],
                                                 tempUnfillableBoxes,
                                                 fillBoxId[ibox][currentLevel][currentBLI],
                                                 tempCoarseBox,
                                                 interpTime, srcComp, destComp, nComp);

                    unfillableBoxesOnThisLevel.join(tempUnfillableBoxes);
                    ++currentBLI;
                }
            }

            unfilledBoxesOnThisLevel.clear();
            unfilledBoxesOnThisLevel =
                unfillableBoxesOnThisLevel.intersect(currentPDomain);

            if (unfilledBoxesOnThisLevel.isEmpty())
            {
                needsFilling = false;
            }
            else
            {
                Box coarseLocalMFBox(localMFBoxes[ibox]);
                coarseLocalMFBox.coarsen(cumulativeRefRatios[currentLevel]);
                if (!is_periodic)
                {
                    unfilledBoxesOnThisLevel.intersect(coarseLocalMFBox);
                }
                unfilledBoxesOnThisLevel.coarsen(amrLevels[currentLevel].crse_ratio);

                if (currentLevel == 0)
                {
                    BoxList unfilledInside =
                        unfilledBoxesOnThisLevel.intersect(currentPDomain);
                    if (!unfilledInside.isEmpty())
                    {
                        unfilledInside.intersect(coarseLocalMFBox);
                        assert(unfilledInside.isEmpty());
                    }
                }
            }
        }
    }

    multiFabCopyDesc.CollectData();

    bIsInitialized = true;
}

bool
FillPatchIterator::isValid (bool bDoSync)
{
    assert(bIsInitialized);
    //
    // If the currentIndex is valid,
    // this function will fill the currentFillPatchedFab from state
    // so it is ready to be used if requested by operator()
    // the fill is done from this level and, if necessary, coarser levels
    // with values from the FillPatchIterator constructor.
    //
    if (!MultiFabIterator::isValid(bDoSync))
        //
        // This does a sync if not valid.
        //
        return false; // (if bDoSync == true)

    int myproc = ParallelDescriptor::MyProc();
    PArray<AmrLevel> &amrLevels = amrLevel.parent->getAmrLevels();

    Box destBox(validbox());
    destBox.grow(growSize);

    currentFillPatchedFab.resize(destBox, nComp);
    currentFillPatchedFab.setVal(1.e30);

    int currentLevel;
    for (currentLevel = 0; currentLevel <= amrLevel.level; ++currentLevel)
    {
        bool is_periodic = amrLevels[currentLevel].geom.isAnyPeriodic();
        StateData& currentState = amrLevels[currentLevel].state[stateIndex];
        const Box& currentPDomain = currentState.getDomain();

        for (int currentBox = 0;
            currentBox < fillBoxId[currentIndex][currentLevel].length();
            ++currentBox)
        {
            Box tempCoarseBox(fillBoxId[currentIndex][currentLevel][currentBox][0].box());
            FArrayBox tempCoarseDestFab(tempCoarseBox, nComp);
            tempCoarseDestFab.setVal(1.e30);
            currentState.linInterpFillFab(multiFabCopyDesc,
                                          stateDataMFId[currentLevel],
                                          fillBoxId[currentIndex][currentLevel][currentBox],
                                          tempCoarseDestFab, interpTime,
                                          0, destComp, nComp);

            Box intersectDestBox(savedFineBox[currentIndex][currentLevel][currentBox]);
            if (!is_periodic)
                intersectDestBox &= currentFillPatchedFab.box();

            const BoxArray& filledBoxes =
                fillBoxId[currentIndex][currentLevel][currentBox][0].FilledBoxes();
            BoxArray fboxes(filledBoxes);
            FArrayBox* copyFromThisFab;
            const BoxArray* copyFromTheseBoxes;
            FArrayBox tempCurrentFillPatchedFab;

            if (intersectDestBox.ok())
            {
                if (currentLevel != amrLevel.level)
                {
                    //
                    // Get boundary conditions for this patch.
                    //
                    Array<BCRec> bcCoarse(nComp);
                    const StateDescriptor &desc = amrLevels[currentLevel].desc_lst[stateIndex];
                    setBC(intersectDestBox, currentPDomain, srcComp, 0, nComp,
                          desc.getBCs(), bcCoarse);
                    fboxes.refine(cumulativeRefRatios[currentLevel]);
                    //
                    // Interpolate up to fine patch.
                    //
                    tempCurrentFillPatchedFab.resize(intersectDestBox, nComp);
                    map[currentLevel]->interp(tempCoarseDestFab, 0,
                                              tempCurrentFillPatchedFab,
                                              destComp, nComp,
                                              intersectDestBox,
                                              cumulativeRefRatios[currentLevel],
                                              amrLevels[currentLevel].geom,
                                              amrLevels[amrLevel.level].geom,
                                              bcCoarse);
                    copyFromThisFab = &tempCurrentFillPatchedFab;
                    copyFromTheseBoxes = &fboxes;
                }
                else
                {
                    copyFromThisFab = &tempCoarseDestFab;
                    copyFromTheseBoxes = &filledBoxes;
                }
                int iFillBox = 0;
                for ( ; iFillBox < copyFromTheseBoxes->length(); ++iFillBox)
                {
                    Box srcdestBox((*copyFromTheseBoxes)[iFillBox]);
                    srcdestBox &= currentFillPatchedFab.box();
                    srcdestBox &= intersectDestBox;
                    if (srcdestBox.ok())
                    {
                        currentFillPatchedFab.copy(*copyFromThisFab,
                                                   srcdestBox, 0, srcdestBox,
                                                   destComp, nComp);
                    }
                }

                if (is_periodic)
                {
                    StateData& fineState   = amrLevels[amrLevel.level].state[stateIndex];
                    const Box& finePDomain = fineState.getDomain();

                    if (!finePDomain.contains(currentFillPatchedFab.box()))
                    {
                        Array<IntVect> pshifts(27);
                        const Box& dbox = currentFillPatchedFab.box();
                        amrLevels[amrLevel.level].geom.periodicShift(finePDomain, dbox, pshifts);
                        for (int iiv = 0; iiv < pshifts.length(); iiv++)
                        {
                            IntVect iv = pshifts[iiv];
                            currentFillPatchedFab.shift(iv);
                            for (int iFillBox = 0;
                                 iFillBox < copyFromTheseBoxes->length();
                                ++iFillBox)
                            {
                                Box srcdestBox((*copyFromTheseBoxes)[iFillBox]);
                                srcdestBox &= currentFillPatchedFab.box();
                                srcdestBox &= copyFromThisFab->box();
                                if (srcdestBox.ok())
                                    currentFillPatchedFab.copy(*copyFromThisFab,
                                                               srcdestBox, 0,
                                                               srcdestBox,
                                                               destComp, nComp);
                            }
                            currentFillPatchedFab.shift(-iv);
                        }
                    }
                }
            }
        }
    }
    //
    // Do non-periodic BCs on the finest level.
    //
    StateData& currentState = amrLevel.state[stateIndex];
    const Box& p_domain = amrLevel.state[stateIndex].getDomain();
    if (!p_domain.contains(destBox))
    {
        const Real* dx = amrLevel.geom.CellSize();
        const RealBox& realProbDomain = amrLevel.geom.ProbDomain();
        currentState.FillBoundary(currentFillPatchedFab, interpTime, dx,
                                  realProbDomain, destComp, srcComp, nComp);
    }

    return true;
}

#else /*BL_NEWFPMINBOX*/

FillPatchIterator::FillPatchIterator (AmrLevel&     amrlevel,
                                      MultiFab&     levelData,
                                      int           boxGrow,
                                      int           dest_comp,
                                      Real          time,
                                      int           state_index,
                                      int           src_comp,
                                      int           ncomp,
                                      Interpolater* mapper)
    :
    MultiFabIterator(levelData),
    amrLevel(amrlevel),
    growSize(boxGrow),
    stateIndex(state_index),
    srcComp(src_comp),
    destComp(dest_comp),
    nComp(ncomp),
    interpTime(time),
    multiFabCopyDesc(true)
{
    //
    // This function sets up and performs the communication pattern for
    // filling fabs of size levelData[i].box().grow(boxGrow) from amrlevel's
    // state data, possibly interpolating between the new and old data
    // the fill is done from this level and, if necessary, coarser levels.
    //
    int myproc = ParallelDescriptor::MyProc();
    int currentLevel;
    PArray<AmrLevel> &amrLevels = amrLevel.parent->getAmrLevels();
    map.resize(amrLevel.level + 1);
    stateDataMFId.resize(amrLevel.level + 1);
    for (currentLevel = 0; currentLevel <= amrLevel.level; ++currentLevel)
    {
        StateData &currentState = amrLevels[currentLevel].state[stateIndex];
        const StateDescriptor &currentDesc =
            amrLevels[currentLevel].desc_lst[stateIndex];
        currentState.RegisterData(multiFabCopyDesc, stateDataMFId[currentLevel]);
        map[currentLevel] = mapper;
        if (map[currentLevel] == 0)
            map[currentLevel] = currentDesc.interp();
    }

    localMFBoxes.resize(levelData.boxArray().length());
    fillBoxId.resize(levelData.boxArray().length());
    savedFineBox.resize(levelData.boxArray().length());
    for (int iLocal = 0; iLocal < localMFBoxes.length(); ++iLocal)
    {
        if (levelData.DistributionMap().ProcessorMap()[iLocal] == myproc)
        {
            //
            // Local.
            //
            localMFBoxes.set(iLocal, levelData.boxArray()[iLocal]);
            fillBoxId[iLocal].resize(amrLevel.level + 1);
            savedFineBox[iLocal].resize(amrLevel.level + 1);
        }
    }
    localMFBoxes.grow(growSize);  // These are the ones we want to fillpatch.

    Box unfilledBoxOnThisLevel;
    BoxList unfillableBoxesOnThisLevel;
    //
    // Do this for all local (grown) fab boxes.
    //
    for (int ibox = 0; ibox < localMFBoxes.length(); ++ibox)
    {
        if (levelData.DistributionMap().ProcessorMap()[ibox] != myproc)
            continue;  // Not local.

        unfilledBoxOnThisLevel = localMFBoxes[ibox] &
            amrLevels[amrLevel.level].state[stateIndex].getDomain();
        assert(unfilledBoxOnThisLevel.ok());
        bool needsFilling = true;

        for (currentLevel = amrLevel.level; currentLevel >= 0 && needsFilling;
            --currentLevel)
        {
            int refRatio = amrLevels[amrLevel.level].crse_ratio;
            refRatio = 2;

            StateData &currentState   = amrLevels[currentLevel].state[stateIndex];
            const Box &currentPDomain = currentState.getDomain();
            unfillableBoxesOnThisLevel.clear();

            fillBoxId[ibox][currentLevel].resize(1);
            savedFineBox[ibox][currentLevel].resize(1);

            int currentBLI = 0;
            savedFineBox[ibox][currentLevel][currentBLI] = unfilledBoxOnThisLevel;

            Box tempCoarseBox;
            if (currentLevel == amrLevel.level)
            {
                tempCoarseBox = unfilledBoxOnThisLevel;
            }
            else
            {
                tempCoarseBox = map[currentLevel]->CoarseBox(unfilledBoxOnThisLevel,
                                                             refRatio);
            }

            currentState.linInterpAddBox(multiFabCopyDesc,
                                         stateDataMFId[currentLevel],
                                         unfillableBoxesOnThisLevel,
                                         fillBoxId[ibox][currentLevel][currentBLI],
                                         tempCoarseBox,
                                         interpTime, srcComp, destComp, nComp);

            unfillableBoxesOnThisLevel.intersect(currentPDomain);
            unfilledBoxOnThisLevel = unfillableBoxesOnThisLevel.minimalBox();
            unfilledBoxOnThisLevel &= currentPDomain;

            if (unfilledBoxOnThisLevel.ok())
            {
            }
            else
            {
                needsFilling = false;
            }
        }
    }

  multiFabCopyDesc.CollectData();
}

bool
FillPatchIterator::isValid (bool bDoSync)
{
    assert(bIsInitialized);
    //
    // if the currentIndex is valid,
    // this function will fill the currentFillPatchedFab from state
    // so it is ready to be used if requested by operator()
    // the fill is done from this level and, if necessary, coarser levels
    // with values from the FillPatchIterator constructor.
    //
    if (!MultiFabIterator::isValid(bDoSync))
        //
        // This does a sync if not valid.
        //
        return false; // (if bDoSync == true)

    int myproc = ParallelDescriptor::MyProc();
    PArray<AmrLevel> &amrLevels = amrLevel.parent->getAmrLevels();

    Box destBox(box());
    destBox.grow(growSize);

    currentFillPatchedFab.resize(destBox, nComp);
    currentFillPatchedFab.setVal(1.e30);

    FArrayBox tempCoarseDestFab, tempFineDestFab;
    FArrayBox *coarseDestFabPtr = 0, *fineDestFabPtr = 0;
    int currentLevel;
    int coarsestFillLevel = amrLevel.level;
    for (currentLevel = 0; currentLevel < amrLevel.level; ++currentLevel)
    {
        if (fillBoxId[currentIndex][currentLevel].length() > 0)
        {
            coarsestFillLevel = currentLevel;
            break;
        }
    }
    assert(coarsestFillLevel >= 0 && coarsestFillLevel <= amrLevel.level);

    for (currentLevel = coarsestFillLevel; currentLevel < amrLevel.level;
        ++currentLevel)
    {
        if (fillBoxId[currentIndex][currentLevel].length() == 0)
            continue;

        assert(fillBoxId[currentIndex][currentLevel].length() == 1);

        int ivRefRatio = 2;
        int currentBox = 0;
        StateData &currentState = amrLevels[currentLevel].state[stateIndex];
        Box tempCoarseBox(fillBoxId[currentIndex][currentLevel][currentBox][0].box());

        if (currentLevel == coarsestFillLevel)
        {
            assert(tempCoarseBox.ok());
            tempCoarseDestFab.resize(tempCoarseBox, nComp);
            tempCoarseDestFab.setVal(1.e30);
            coarseDestFabPtr = &tempCoarseDestFab;
        }
        else
        {
            assert(fineDestFabPtr != 0);
            coarseDestFabPtr = fineDestFabPtr;
        }
        assert(coarseDestFabPtr != 0);

        currentState.linInterpFillFab(multiFabCopyDesc,
                                      stateDataMFId[currentLevel],
                                      fillBoxId[currentIndex][currentLevel][currentBox],
                                      *coarseDestFabPtr,
                                      interpTime, srcComp, destComp, nComp);

        const Real *dx = amrLevels[currentLevel].geom.CellSize();
        const RealBox &realProbDomain = amrLevels[currentLevel].geom.ProbDomain();

        currentState.FillBoundary(*coarseDestFabPtr, interpTime, dx,
                                  realProbDomain, destComp, srcComp, nComp);

        const Box &currentPDomain = currentState.getDomain();

        Box intersectDestBox(savedFineBox[currentIndex][currentLevel][currentBox]);

        assert(intersectDestBox.ok());
        Array<BCRec> bcCoarse(nComp);
        const StateDescriptor &desc = amrLevels[currentLevel].desc_lst[stateIndex];
        setBC(intersectDestBox, currentPDomain, srcComp, 0, nComp,
              desc.getBCs(), bcCoarse);

        // interpolate up to fine patch
        const BoxArray &filledBoxes =
            fillBoxId[currentIndex][currentLevel][currentBox][0].FilledBoxes();
        BoxArray fboxes(filledBoxes);
        fboxes.refine(ivRefRatio);
        assert(fboxes.length() == 1);
        int iFillBox = 0;
        Box srcdestBox(fboxes[iFillBox]);
        srcdestBox &= currentFillPatchedFab.box();
        srcdestBox &= intersectDestBox;
        if ((currentLevel + 1) == amrLevel.level)
        {
            fineDestFabPtr = &currentFillPatchedFab;
        }
        else
        {
            tempFineDestFab.resize(fillBoxId[currentIndex][currentLevel+1][currentBox][0
                ].box(), nComp);
            fineDestFabPtr = &tempFineDestFab;
        }

        map[currentLevel]->interp(*coarseDestFabPtr, 0, *fineDestFabPtr,
                                  destComp, nComp, intersectDestBox,
                                  ivRefRatio,
                                  amrLevels[currentLevel].geom,
                                  amrLevels[currentLevel + 1].geom,
                                  bcCoarse);
    }

    currentLevel = amrLevel.level;
    int currentBox = 0;
    StateData &currentState = amrLevel.state[stateIndex];
    currentState.linInterpFillFab(multiFabCopyDesc,
                                  stateDataMFId[currentLevel],
                                  fillBoxId[currentIndex][currentLevel][currentBox],
                                  currentFillPatchedFab,
                                  interpTime, srcComp, destComp, nComp);
    //
    // Do non-periodic BCs on the finest level.
    //
    const Box &p_domain = amrLevel.state[stateIndex].getDomain();
    if (!p_domain.contains(destBox))
    {
        const Real *dx = amrLevel.geom.CellSize();
        const RealBox &realProbDomain = amrLevel.geom.ProbDomain();

        currentState.FillBoundary(currentFillPatchedFab, interpTime, dx,
                                  realProbDomain, destComp, srcComp, nComp);
    }

    return true;
}
#endif /*!BL_NEWFPMINBOX*/

FillPatchIterator::~FillPatchIterator () {}

void
AmrLevel::FillCoarsePatch (MultiFab&     mfdest,
                           int           dest_comp,
                           Real          time,
                           int           state_indx,
                           int           src_comp,
                           int           ncomp,
                           Interpolater* mapper)
{
    //
    // Must fill this region on crse level and interpolate.
    //
    assert(level != 0);
    assert(ncomp <= (mfdest.nComp() - dest_comp));
    assert((0 <= state_indx) && (state_indx < desc_lst.length()));

    const StateDescriptor& desc = desc_lst[state_indx];
    assert(desc.inRange(src_comp, ncomp));
    Interpolater* map = mapper;
    if (map == 0)
        map = desc.interp();

    const RealBox& prob_domain = geom.ProbDomain();
    const Box& p_domain = state[state_indx].getDomain();
    AmrLevel& crse_lev = parent->getLevel(level-1);
    //
    // Build a properly coarsened boxarray.
    //
    BoxArray mfdestBoxArray(mfdest.boxArray());
    BoxArray crse_regBoxArray(mfdestBoxArray.length());
    for (int ibox = 0; ibox < mfdestBoxArray.length(); ++ibox)
    {
        Box dbox(mfdest.fabbox(ibox));
        assert(dbox.ixType() == desc.getType());
        dbox &= p_domain;
        //
        // Coarsen unfilled region and widen by interpolater stencil width.
        //
        Box crse_reg(map->CoarseBox(dbox,crse_ratio));
        crse_regBoxArray.set(ibox, crse_reg);
    }

    int boxGrow = 0;
    MultiFab mf_crse_reg(crse_regBoxArray, ncomp, boxGrow);

    for (FillPatchIterator fpi(crse_lev, mf_crse_reg, boxGrow, dest_comp,
                               time, state_indx, src_comp, ncomp, mapper);
         fpi.isValid();
         ++fpi)
    {
        DependentMultiFabIterator mfdest_mfi(fpi, mfdest);
        assert(mfdest_mfi.fabbox() == mfdest_mfi().box());
        Box dbox(mfdest_mfi().box());
        dbox &= p_domain;

        FArrayBox& crse = fpi();
        FArrayBox& dest = mfdest_mfi();
        //
        // Get bndry conditions for this patch.
        //
        Array<BCRec> bc_crse(ncomp);
        setBC(dbox,p_domain,src_comp,0,ncomp,desc.getBCs(),bc_crse);
        //
        // Interpolate up to fine patch.
        //
        const Geometry& crse_geom = crse_lev.geom;
        map->interp(crse,0,dest,dest_comp,ncomp,dbox,
                    crse_ratio,crse_geom,geom,bc_crse);

        if (!p_domain.contains(dbox)) 
        {
            const Real *dx = geom.CellSize();
            state[state_indx].FillBoundary(dest,time,dx,prob_domain,
                                           dest_comp,src_comp,ncomp);
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
                              0, time, state_indx, src_comp, 1);

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

            FillPatchIterator fpi(*this, srcMF, ngrow, 0, time, state_indx,
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
