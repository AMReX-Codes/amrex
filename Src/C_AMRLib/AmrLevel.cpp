//BL_COPYRIGHT_NOTICE

//
// $Id: AmrLevel.cpp,v 1.19 1998-02-18 21:35:32 vince Exp $
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
#ifdef BL_PARALLEL_IO
        state[i].restart(is, desc_lst[i], papa.theRestartFile());
#else
        state[i].restart(is, desc_lst[i]);
#endif
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

#ifdef BL_PARALLEL_IO
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
    ParallelDescriptor::Synchronize();

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

#else

void
AmrLevel::checkPoint (ostream& os)
{
    int ndesc = desc_lst.length(), i;

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
        state[i].checkPoint(os);
}
#endif /*BL_PARALLEL_IO*/

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

#if (NEWFPMINBOX == 0)
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

  //cout << "_in FillPatchIterator::Init():  CopyDescriptor stats: " << endl;
  //multiFabCopyDesc.PrintStats();

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

#endif /*NEWFPMINBOX==0*/

#if (NEWFPMINBOX == 1)
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

        //cout << "dest.box() = " << localMFBoxes[ibox] << '\n';
        unfilledBoxOnThisLevel = localMFBoxes[ibox] &
            amrLevels[amrLevel.level].state[stateIndex].getDomain();
        assert(unfilledBoxOnThisLevel.ok());
        bool needsFilling = true;

        for (currentLevel = amrLevel.level; currentLevel >= 0 && needsFilling;
            --currentLevel)
        {
            //cout << '\n' << "currentLevel = " << currentLevel << '\n';
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
                //cout << "currentLevel = " << currentLevel << '\n';
                //cout << "mapped->CoarseBox = " << tempCoarseBox << "  from box " << unfilledBoxOnThi
                //sLevel << '\n';
            }
            //cout << "about to linInterp with tempCoarseBox = " << tempCoarseBox << '\n';

            currentState.linInterpAddBox(multiFabCopyDesc,
                                         stateDataMFId[currentLevel],
                                         unfillableBoxesOnThisLevel,
                                         fillBoxId[ibox][currentLevel][currentBLI],
                                         tempCoarseBox,
                                         interpTime, srcComp, destComp, nComp);

            //cout << "      fillBoxId.box[]  = " << fillBoxId[ibox][currentLevel][currentBLI][0].
            //box()  << '\n';

            unfillableBoxesOnThisLevel.intersect(currentPDomain);
            unfilledBoxOnThisLevel = unfillableBoxesOnThisLevel.minimalBox();
            unfilledBoxOnThisLevel &= currentPDomain;

            if (unfilledBoxOnThisLevel.ok())
            {
                //cout << "unfilled box on level " << currentLevel << " = " << unfilledBoxOnThisLevel
                     //<< '\n';
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
    //cout << "coarsestFillLevel = " << coarsestFillLevel << '\n';
    assert(coarsestFillLevel >= 0 && coarsestFillLevel <= amrLevel.level);

    for (currentLevel = coarsestFillLevel; currentLevel < amrLevel.level;
        ++currentLevel)
    {
        if (fillBoxId[currentIndex][currentLevel].length() == 0)
            continue;

        assert(fillBoxId[currentIndex][currentLevel].length() == 1);
        //cout << "currentLevel     = " << currentLevel     << '\n';

        int ivRefRatio = 2;
        int currentBox = 0;
        StateData &currentState = amrLevels[currentLevel].state[stateIndex];
        Box tempCoarseBox(fillBoxId[currentIndex][currentLevel][currentBox][0].box());

        if (currentLevel == coarsestFillLevel)
        {
            assert(tempCoarseBox.ok());
            //cout << "Resizing coarse fab to " << tempCoarseBox << '\n';
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

        //cout << "linInterp on coarse.box() = " << coarseDestFabPtr->box() << '\n';
        currentState.linInterpFillFab(multiFabCopyDesc,
                                      stateDataMFId[currentLevel],
                                      fillBoxId[currentIndex][currentLevel][currentBox],
                                      *coarseDestFabPtr,
                                      interpTime, srcComp, destComp, nComp);

        const Real *dx = amrLevels[currentLevel].geom.CellSize();
        const RealBox &realProbDomain = amrLevels[currentLevel].geom.ProbDomain();

        //cout << "FillBoundary on coarse.box() = " << coarseDestFabPtr->box()
             //<< "  outside boundary " << realProbDomain << '\n';
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

        //cout << "map->interp coarse.box() = " << coarseDestFabPtr->box()
             //<< " to finer box = " << fineDestFabPtr->box()
             //<< " on int_region = " << intersectDestBox << '\n';
        //cout << "crse0_geom = " << amrLevels[currentLevel].geom << '\n';
        //cout << "crse_geom  = " << amrLevels[currentLevel + 1].geom  << '\n';
        //for(int ibc=0; ibc < nComp; ++ibc) {
            //cout << "bc_crse[" << ibc << "] = " << bcCoarse[ibc] << '\n';
        //}
        map[currentLevel]->interp(*coarseDestFabPtr, 0, *fineDestFabPtr,
                                  destComp, nComp, intersectDestBox,
                                  ivRefRatio,
                                  amrLevels[currentLevel].geom,
                                  amrLevels[currentLevel + 1].geom,
                                  bcCoarse);
    }

    currentLevel = amrLevel.level;
    int currentBox = 0;
    //cout << "linInterp on currentFillPatched.box() = "
         //<< currentFillPatchedFab.box() << '\n';
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

        //cout << "FillBoundary on dest.box() = "
             //<< currentFillPatchedFab.box() << " outside boundary "
             //<< realProbDomain << '\n';
        currentState.FillBoundary(currentFillPatchedFab, interpTime, dx,
                                  realProbDomain,
                                  destComp, srcComp, nComp);
    }

    return true;
}
#endif /*NEWFPMINBOX==1*/

FillPatchIterator::~FillPatchIterator () {}

#if (USEUNRAVELEDFILLPATCH == 1)
void
AmrLevel::FillPatch (FArrayBox&    dest,
                     int           dest_comp,
                     Real          time,
                     int           stateIndex,
                     int           src_comp,
                     int           ncomp,
                     Interpolater* mapper)
{
cout << '\n';
cout << "_in old FillPatch" << '\n';
cout << "dest.box() = " << dest.box() << '\n';
    dest.setVal(1.e30);
    const Box& dbox = dest.box();
    Box truncdbox(dbox);
    const StateDescriptor& desc = desc_lst[stateIndex];
    const RealBox& prob_domain = geom.ProbDomain();
    const Box& p_domain = state[stateIndex].getDomain();

    bool inside = p_domain.contains(dbox); // Does grid intersect boundary.
    if (!inside )
        truncdbox &= p_domain;

    Box unfilled_region;
    BoxDomain fd(dbox.ixType());
    fd.add(truncdbox);

    Box enclosing_box = fd.minimalBox();
    const BoxArray &grds = state[stateIndex].boxArray();
    int i;
    for (i = 0; i < grds.length(); i++)
    {
        const Box& gbox = grds[i];
        if (enclosing_box.intersects(gbox))
            fd.rmBox(gbox);
    }
    unfilled_region = fd.minimalBox();
    fd.clear();

    FArrayBox crse;
    if (unfilled_region.ok())
    {
      AmrLevel &crse_lev = parent->getLevel(level-1);
      const Geometry& crse_geom = crse_lev.geom;
      //
      // Must fill on this region on crse level and interpolate.
      //
      assert(level != 0);

      Interpolater *map = mapper;
      if (map == 0)
          map = desc.interp();
      //
      // Intersect unfilled_region with the domain, this is necessary
      // because the unfilled_region may contain stuff carried by
      // periodic BC to be outside the domain of dest.
      //
      Box int_region = unfilled_region & dbox;
      if (int_region.ok())
      {
cout << "unfilled box on level " << level << " = " << int_region << '\n';
        //
        // Coarsen unfilled region and widen if necessary.
        //
        Box crse_reg(map->CoarseBox(int_region,crse_ratio));
cout << "mapped->CoarseBox = " << crse_reg << '\n';
cout << "Resizing crse fab to " << crse_reg << '\n';
        crse.resize(crse_reg,ncomp);         // alloc patch for crse level
        crse.setVal(1.e30);
        Box crse_dbox(crse.box());
        Box crse_truncdbox(crse_dbox);
        const RealBox &crse_prob_domain = crse_geom.ProbDomain();
        const Box &crse_p_domain = crse_lev.state[stateIndex].getDomain();
        int crse_inside = crse_p_domain.contains(crse_dbox);

        if (!inside)
            crse_truncdbox &= crse_p_domain;

        Box crse_unfilled_region;
        BoxDomain crse_fd(crse_dbox.ixType());
        crse_fd.add(crse_truncdbox);

        Box crse_enclosing_box = crse_fd.minimalBox();
        const BoxArray &crse_grds = crse_lev.state[stateIndex].boxArray();
        for (i = 0; i < crse_grds.length(); i++)
        {
            const Box& crse_gbox = crse_grds[i];
            if (crse_enclosing_box.intersects(crse_gbox))
                crse_fd.rmBox(crse_gbox);
        }
        crse_unfilled_region = crse_fd.minimalBox();
        crse_fd.clear();

        FArrayBox crse0;
        if (crse_unfilled_region.ok())
        {
          AmrLevel &crse0_lev = parent->getLevel(level-2);  // should be 0 here
          const Geometry& crse0_geom = crse0_lev.geom;
          //
          // Must fill on this region on crse0 level and interpolate.
          //
          Box crse_int_region = crse_unfilled_region & crse_dbox;
          if (crse_int_region.ok())
          {
cout << "unfilled box on level " << level-1 << " = " << crse_int_region << '\n';
            //
            // Coarsen unfilled region and widen if necessary.
            //
            int crse0_ratio = 2;  // for testing
            Box crse0_reg(map->CoarseBox(crse_int_region,crse0_ratio));
cout << "mapped->CoarseBox = " << crse0_reg << '\n';
cout << "Resizing crse0 fab to " << crse0_reg << '\n';
            crse0.resize(crse0_reg,ncomp);         // alloc patch for crse0 level

            crse0.setVal(1.e30);
            Box crse0_dbox(crse0.box());
            Box crse0_truncdbox(crse0_dbox);

            const RealBox &crse0_prob_domain = crse0_geom.ProbDomain();
            const Box &crse0_p_domain = crse0_lev.state[stateIndex].getDomain();
            int crse0_inside = crse0_p_domain.contains(crse0_dbox);

            if (!crse0_inside)
                crse0_truncdbox &= crse0_p_domain;
cout << "linInterp on crse0.box() = " << crse0.box() << '\n';
            crse0_lev.state[stateIndex].linInterp(crse0,crse0.box(),time,
                                                 src_comp,0,ncomp);
            if (!crse0_inside)
            {
              const Real* crse0_dx = crse0_geom.CellSize();
cout << "FillBoundary on crse0.box() = " << crse0.box() << "  outside boundary " <<
crse0_prob_domain << '\n';
              crse0_lev.state[stateIndex].FillBoundary(crse0, time, crse0_dx,
                                      crse0_prob_domain, 0, src_comp, ncomp);
            }
            Array<BCRec> bc_crse(ncomp); // get bndry conditions for this patch
            setBC(crse_int_region,crse_p_domain,src_comp,0,
                     ncomp,desc.getBCs(),bc_crse);
            //
            // Interpolate up to fine patch.
            //
cout << "map->interp crse0.box() = " << crse0.box() << " to crse.box() = " << crse.b
ox() << " on crse_int_region = " << crse_int_region << '\n';
cout << "crse0_geom = " << crse0_geom << '\n';
cout << "crse_geom  = " << crse_geom  << '\n';
for (int ibc=0; ibc < ncomp; ++ibc) {
  cout << "bc_crse[" << ibc << "] = " << bc_crse[ibc] << '\n';
}
            map->interp(crse0,0,crse,dest_comp,ncomp,crse_int_region,
                        crse_ratio,crse0_geom,geom,bc_crse);
          }
        }

cout << "linInterp on crse.box() = " << crse.box() << '\n';
        crse_lev.state[stateIndex].linInterp(crse,crse.box(),time,
                                             src_comp,dest_comp,ncomp);
        if (!inside)
        {
            //
            // Do non-periodic BC's on this level.
            //
          const Real* crse_dx = crse_geom.CellSize();
cout << "FillBoundary on crse.box() = " << crse.box() << "  outside boundary " << cr
se_prob_domain << '\n';
          crse_lev.state[stateIndex].FillBoundary(crse,time,crse_dx,
                                                  crse_prob_domain,
                                                  dest_comp,src_comp,ncomp);
        }
        Array<BCRec> bc_crse(ncomp); // Get bndry conditions for this patch.
        setBC(int_region,p_domain,src_comp,0,ncomp,desc.getBCs(),bc_crse);
        //
        // Interpolate up to fine patch.
        //
cout << "map->interp crse.box() = " << crse.box() << " to dest.box() = " << dest.box
() << " on int_region = " << int_region << '\n';
cout << "crse_geom  = " << crse_geom  << '\n';
cout << "geom       = " << geom       << '\n';
for(int ibc=0; ibc < ncomp; ++ibc) {
  cout << "bc_crse[" << ibc << "] = " << bc_crse[ibc] << '\n';
}
        map->interp(crse,0,dest,dest_comp,ncomp,int_region,
                    crse_ratio,crse_geom,geom,bc_crse);
      }
    }

  state[stateIndex].linInterp(dest,dest.box(),time,src_comp,dest_comp,ncomp);

  if (!inside)
  {
      //
      // Do non-periodic BC's on this level.
      //
    const Real* dx = geom.CellSize();
cout << "FillBoundary on dest.box() = " << dest.box() << "  outside boundary " << pr
ob_domain << '\n';
    state[stateIndex].FillBoundary(dest,time,dx,prob_domain,
                                   dest_comp,src_comp,ncomp);
  }
cout << "_out old FillPatch" << '\n';
cout << '\n';
}
#endif /*USEUNRAVELEDFILLPATCH==1*/

#if (USEUNRAVELEDFILLPATCH == 0)
//
// old filPatch  (recursive)
//
void
AmrLevel::FillPatch (FArrayBox&    dest,
                     int           dest_comp,
                     Real          time,
                     int           state_indx,
                     int           src_comp,
                     int           ncomp,
                     Interpolater* mapper)
{
cout << '\n';
cout << "]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]] _in old FillPatch" << '\n';
cout << "currentLevel     = " << this->Level() << '\n';
    dest.setVal(1.e30);
    const Box& dbox = dest.box();
cout << "box to FillPatch = " << dest.box() << '\n';
    Box truncdbox(dbox);
    int nv = dest.nComp();
    int ndesc = desc_lst.length();
    const StateDescriptor& desc = desc_lst[state_indx];
    const RealBox& prob_domain = geom.ProbDomain();
    const Box& p_domain = state[state_indx].getDomain();

    bool inside = p_domain.contains(dbox);
    if (!inside )
        truncdbox &= p_domain;

    Box unfilled_region;
    BoxDomain fd(dbox.ixType());
    fd.add(truncdbox);

    Box enclosing_box = fd.minimalBox();
    //
    // Step 2: take away stuff that can be gotten from this level of refinement
    //
    const BoxArray &grds = state[state_indx].boxArray();
    for (int i = 0; i < grds.length(); i++)
    {
        if (enclosing_box.intersects(grds[i]))
            fd.rmBox(grds[i]);
    }
    unfilled_region = fd.minimalBox();
    fd.clear();

    FArrayBox crse;
    if (unfilled_region.ok())
    {
      AmrLevel &crse_lev = parent->getLevel(level-1);
      const Geometry& crse_geom = crse_lev.geom;
      //
      // Must fill on this region on crse level and interpolate.
      //
      assert(level != 0);

      Interpolater *map = mapper;
      if (map == 0)
          map = desc.interp();

      Box int_region = unfilled_region & dbox;
      if (int_region.ok())
      {
          //
          // Coarsen unfilled region and widen if necessary.
          //
          Box crse_reg(map->CoarseBox(int_region,crse_ratio));
          //
          // alloc patch for crse level.
          //
          crse.resize(crse_reg,ncomp);
          //
          // Fill patch at lower level.
          //
          crse_lev.FillPatch(crse,0,time,state_indx,src_comp,ncomp,mapper);
          //
          // Get bndry conditions for this patch.
          //
          Array<BCRec> bc_crse(ncomp);
          setBC(int_region,p_domain,src_comp,0,ncomp,
                desc.getBCs(),bc_crse);

cout << "]]----------- about to map->interp:" << '\n';
cout << "]]----------- crse.box()       = " << crse.box() << '\n';
cout << "]]----------- dest.box()       = " << dest.box() << '\n';
cout << "]]----------- intersectDestBox = " << int_region << '\n';
        //
        // Interpolate up to fine patch.
        //
        map->interp(crse,0,dest,dest_comp,ncomp,int_region,
                    crse_ratio,crse_geom,geom,bc_crse);
      }
    }

cout << "]]=========== about to linInterp:" << '\n';
cout << "]]=========== tempCoarseDestFab.box() = " << dest.box() << '\n';
    //
    // Copy from data on this level.
    //
    state[state_indx].linInterp(dest,dest.box(),time,src_comp,dest_comp,ncomp);
    //
    // Do non-periodic BC's on this level.
    //
    if (!inside)
    {
        const Real* dx = geom.CellSize();
        state[state_indx].FillBoundary(dest,time,dx,prob_domain,
                                   dest_comp,src_comp,ncomp);
    }
cout << "]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]] _out old FillPatch" << '\n';
cout << '\n';
}
#endif /*USEUNRAVELEDFILLPATCH==0*/

void
AmrLevel::FillCoarsePatch (FArrayBox&   dest,
                           int          dest_comp,
                           Real         time,
                           int          state_indx,
                           int          src_comp,
                           int          ncomp,
                          Interpolater* mapper)
{
    //
    // Must fill this region on crse level and interpolate.
    //
    if (level == 0)
        BoxLib::Error("crsePatch: called at level 0");

    Box dbox = dest.box();

    int nv = dest.nComp();
    int ndesc = desc_lst.length();

    assert((0<=state_indx) && (state_indx<ndesc));
    const StateDescriptor& desc = desc_lst[state_indx];

    assert(dbox.ixType() == desc.getType());
    assert(desc.inRange(src_comp, ncomp));
    assert(ncomp <= (nv-dest_comp));

    const RealBox& prob_domain = geom.ProbDomain();

    const Box& p_domain = state[state_indx].getDomain();
    //
    // Intersect with problem domain at this level.
    //
    dbox &= p_domain;

    Interpolater* map = mapper;
    if (map == 0)
        map = desc.interp();
    //
    // Coarsen unfilled region and widen by interpolater stencil width.
    //
    Box crse_reg(map->CoarseBox(dbox,crse_ratio));
    //
    // Alloc patch for crse level.
    //
    FArrayBox crse(crse_reg,ncomp);
    //
    // Fill patch at lower level.
    //
    AmrLevel &crse_lev = parent->getLevel(level-1);
    crse_lev.FillPatch(crse,0,time,state_indx,src_comp,ncomp,mapper);
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
        const Real* dx = geom.CellSize();
        state[state_indx].FillBoundary(dest,time,dx,prob_domain,
                                   dest_comp,src_comp,ncomp);
    }
}

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
                  Real           time)
{
    BoxLib::Error("AmrLevel::derive(MultiFab*) not implemented");
    //
    // TODO -- implement this!!!
    //
    // The below code returns PArray<FArrayBox*> which isn't OK for parallel.
    //
#if 0
    int state_indx, src_comp;

    if (isStateVariable(name,state_indx,src_comp))
    {
        const StateDescriptor& desc = desc_lst[state_indx];
        int nc = desc.nComp();
        const BoxArray& grds = state[state_indx].boxArray();
        PArray<FArrayBox> *df = new PArray<FArrayBox>(grids.length(),PArrayManage);
        for (int i = 0; i < grds.length(); i++)
        {
            FArrayBox* dest = new FArrayBox(grds[i],1);
            state[state_indx].linInterp(*dest,grds[i],time,src_comp,0,1);
            df->set(i,dest);
        }
        return df;
    }
    //
    // Can quantity be derived?
    //
    const DeriveRec* d;
    if (d = derive_lst.get(name))
    {
        PArray<FArrayBox> *df = new PArray<FArrayBox>(grids.length(),PArrayManage);

        const Real* dx = geom.CellSize();
        int state_indx, src_comp, num_comp;
        d->getRange(0,state_indx,src_comp,num_comp);
        const BoxArray& grds = state[state_indx].boxArray();
        int n_state = d->numState();
        int n_der = d->numDerive();
        int nsr = d->numRange();
        IndexType der_typ = d->deriveType();
        //
        // Can do special fill
        int i;
        for (i = 0; i < grds.length(); i++)
        {
            //
            // Build destination FAB and install.
            //
            Box dbox(grids[i]);
            dbox.convert(der_typ);
            FArrayBox *dest = new FArrayBox(dbox,n_der);
            df->set(i,dest);
            //
            // Build src fab and fill with component state data.
            //
            Box sbox(d->boxMap()(dbox));
            FArrayBox src(sbox,n_state);
            int dc = 0;
            int k;
            for (k = 0; k < nsr; k++)
            {
                d->getRange(k,state_indx,src_comp,num_comp); 
                const StateDescriptor &desc = desc_lst[state_indx];
                if (grds[i].contains(sbox))
                {
                    //
                    // Can do copy.
                    //
                    state[state_indx].linInterp(src,sbox,time,
                                                src_comp,dc,num_comp);
                }
                else
                {
                    //
                    // Must filpatch.
                    //
                    FillPatch(src,dc,time,state_indx,src_comp,num_comp);
                }
                dc += num_comp;
            }
            //
            // Call deriving function.
            //
            Real *ddat = dest->dataPtr();
            const int* dlo = dest->loVect();
            const int* dhi = dest->hiVect();
            Real *cdat = src.dataPtr();
            const int* clo = src.loVect();
            const int* chi = src.hiVect();
            const int* dom_lo = state[state_indx].getDomain().loVect();
            const int* dom_hi = state[state_indx].getDomain().hiVect();
            const int* bcr = d->getBC();
            const Real* xlo = grid_loc[i].lo();
            d->derFunc()(ddat,ARLIM(dlo),ARLIM(dhi),&n_der,
                         cdat,ARLIM(clo),ARLIM(chi),&n_state,
                         dlo,dhi,dom_lo,dom_hi,dx,xlo,&time,bcr,
                         &level,&i);
        }
        return df;
    }
    //
    // If we got here, cannot derive given name.
    //
    aString msg("AmrLevel::derive(MultiFab*): unknown variable: ");
    msg += name;
    BoxLib::Error(msg.c_str());
#endif
    //
    // Just to keep the compiler happy
    //
    return 0;
}

FArrayBox*
AmrLevel::derive (const Box&     b,
                  const aString& name,
                  Real           time)
{
    BoxLib::Warning("AmrLevel::derive(FAB*) not implemented in parallel");

    if (ParallelDescriptor::NProcs() > 1)
        BoxLib::Error();

    int state_indx, src_comp;
    //
    // Is it a state variable?
    //
    if (isStateVariable(name,state_indx,src_comp))
    {
        FArrayBox *dest = new FArrayBox(b,1);
        FillPatch(*dest,0,time,state_indx,src_comp,1);
        return dest;
    }

    const DeriveRec* d;
    if (d = derive_lst.get(name))
    {
        FArrayBox *dest = new FArrayBox(b,d->numDerive());
        FillDerive(*dest,b,name,time);
        return dest;
    }
    //
    // If we got here, cannot derive given name.
    //
    aString msg("AmrLevel::derive(FAB): unknown variable: ");
    msg += name;
    BoxLib::Error(msg.c_str());
    //
    // Just to keep the compiler happy.
    //
    return 0;
}

void
AmrLevel::FillDerive (FArrayBox&     dest,
                      const Box&     subbox,
                      const aString& name,
                      Real           time)
{
    const DeriveRec* d = derive_lst.get(name);
    assert (d != 0);
    //
    // Get domain and BoxArray.
    //
    IndexType der_typ = d->deriveType();
    int n_der = d->numDerive();
    int cell_centered = der_typ.cellCentered();

    Box dom(geom.Domain());
    if (!cell_centered)
        dom.convert(der_typ);
    //
    // Only fill portion on interior of domain.
    //
    Box dbox(subbox);
    dbox &= dom;
    //
    // create a FIDIL domain to keep track of what can be filled at this level.
    //
    Box unfilled_region;
    BoxDomain fd(dbox.ixType());
    fd.add(dbox);
    int i;
    for (i = 0; i < grids.length(); i++)
    {
        Box gbox(grids[i]);
        if (!cell_centered)
            gbox.convert(der_typ);
        if (dbox.intersects(gbox))
            fd.rmBox(gbox);
    }
    unfilled_region = fd.minimalBox();
    fd.clear();

    if (unfilled_region.ok())
    {
        //
        // Must fill on this region on crse level and interpolate.
        //
        if (level == 0)
            BoxLib::Error("FillPatch: unfilled region at level 0");
        //
        // Coarsen unfilled region and widen if necessary.
        //
        Interpolater *mapper = d->interp();
        Box crse_reg(mapper->CoarseBox(unfilled_region,crse_ratio));
        //
        // Alloc patch for crse level.
        //
        FArrayBox crse(crse_reg,n_der);
        //
        // Fill patch at lower level.
        //
        AmrLevel &crse_lev = parent->getLevel(level-1);
        crse_lev.FillDerive(crse,crse_reg,name,time);
        //
        // Get bncry conditions for this patch.
        //
        Array<BCRec> bc_crse;
        //
        // Interpolate up to fine patch.
        //
        const Geometry& crse_geom = crse_lev.geom;
        mapper->interp(crse,0,dest,0,n_der,unfilled_region,
                       crse_ratio,crse_geom,geom,bc_crse);
    }
    //
    // Walk through grids, deriving on intersect.
    //
    int n_state = d->numState();
    int nsr = d->numRange();
    const Real* dx = geom.CellSize();
    Real dt = parent->dtLevel(level);
    for (i = 0; i < grids.length(); i++)
    {
        Box g(grids[i]);
        if (!cell_centered)
            g.convert(der_typ);
        g &= dbox;
        if (g.ok())
        {
            Box sbox(d->boxMap()(g));
            FArrayBox src(sbox,n_state);
            int dc = 0;
            int state_indx, sc, nc;
            int k;
            for (k = 0; k < nsr; k++)
            {
                d->getRange(k,state_indx,sc,nc);
                const Box& sg = state[state_indx].boxArray()[i];
                if (sg.contains(sbox))
                {
                    //
                    // Can do copy.
                    //
                    state[state_indx].linInterp(src,sbox,time,sc,dc,nc);
                }
                else
                {
                    //
                    // Must filpatch.
                    //
                    FillPatch(src,dc,time,state_indx,sc,nc);
                }
                dc += nc;
            }
            //
            // Now derive.
            //
            Real *ddat = dest.dataPtr();
            const int* dlo = dest.loVect();
            const int* dhi = dest.hiVect();
            const int* lo = g.loVect();
            const int* hi = g.hiVect();
            Real *cdat = src.dataPtr();
            const int* clo = src.loVect();
            const int* chi = src.hiVect();
            const int* dom_lo = dom.loVect();
            const int* dom_hi = dom.hiVect();
            const int* bcr = d->getBC();
            Real xlo[BL_SPACEDIM];
            geom.LoNode(dest.smallEnd(),xlo);
            d->derFunc()(ddat,ARLIM(dlo),ARLIM(dhi),&n_der,
                         cdat,ARLIM(clo),ARLIM(chi),&n_state,
                         lo,hi,dom_lo,dom_hi,dx,xlo,&time,&dt,bcr,
                         &level,&i);
        }
    }
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

void
AmrLevel::probe (ostream&       os,
                 const IntVect& iv,
                 int            rad,
                 Real           time,
                 int            state_indx,
                 int            src_comp,
                 int            num_comp)
{
    Box bx(iv,iv);
    bx.grow(rad);
    FArrayBox fab(bx,num_comp);
    FillPatch(fab,0,time,state_indx,src_comp,num_comp);
    os << "----------------------------------------------------\n"
       << " >>>> PROBE of:";
    const StateDescriptor &desc = desc_lst[state_indx];
    desc.dumpNames(os,src_comp,num_comp);
    os << '\n';
    IntVect lo = bx.smallEnd();
    IntVect hi = bx.bigEnd();
    IntVect point;
    char buf[80];
    for (point=lo; point <= hi; bx.next(point))
    {
        os << point;
        int k;
        for (k = 0; k < num_comp; k++)
        {
            sprintf(buf,"    %20.14f",fab(point,k));
            os << buf;
        }
        os << '\n';
    }
    os << '\n';
    os << "----------------------------------------------------\n";
}
