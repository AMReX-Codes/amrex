#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBGEOMDEBUGOUT_H_
#define _EBGEOMDEBUGOUT_H_

#include "AMReX_IFData.H"
#include "AMReX_CutCellMoments.H"


void dumpIFData1(const IFData<1>* a_ifData);
void dumpIFData2(const IFData<2>* a_ifData);
void dumpIFData3(const IFData<3>* a_ifData);

void dumpCCM1(const CutCellMoments<1>* a_ccm);
void dumpCCM2(const CutCellMoments<2>* a_ccm);
void dumpCCM3(const CutCellMoments<3>* a_ccm);


#endif
