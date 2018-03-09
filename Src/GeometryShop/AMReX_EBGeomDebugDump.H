#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBGEOMDEBUGDUMP_H_
#define _EBGEOMDEBUGDUMP_H_

#include "EBGeomDebugOut.H"

#include "NamespaceHeader.H"

void EBGeomDontYouDareCallMe()
{
  dumpIFData1(NULL);
  dumpIFData2(NULL);
  dumpIFData3(NULL);

  dumpCCM1(NULL);
  dumpCCM2(NULL);
  dumpCCM3(NULL);
}

#include "NamespaceFooter.H"

#endif
