#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMReX_EBGeomDebugOut.H"



void dumpIFData1(const IFData<1>* a_ifData)
{
  a_ifData->print(pout());
}

void dumpIFData2(const IFData<2>* a_ifData)
{
  a_ifData->print(pout());
}

void dumpIFData3(const IFData<3>* a_ifData)
{
  a_ifData->print(pout());
}

void dumpCCM1(const CutCellMoments<1>* a_ccm)
{
  a_ccm->print(pout());
}

void dumpCCM2(const CutCellMoments<2>* a_ccm)
{
  a_ccm->print(pout());
}

void dumpCCM3(const CutCellMoments<3>* a_ccm)
{
  a_ccm->print(pout());
}

