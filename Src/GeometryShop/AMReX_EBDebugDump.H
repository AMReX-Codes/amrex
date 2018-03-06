#ifndef _EBDebugDump___H_
#define _EBDebugDump___H_

#include "AMReX_EBDebugOut.H"

namespace amrex
{

  void dontCallThisFunctionEver()
  {

    dumpFAB(NULL);

    dumpDBL(NULL);

    dumpBFR(NULL);

    dumpBFI(NULL);
 
    dumpBA(NULL);

    dumpIVS(NULL);

    dumpBox(NULL);

    dumpEBFAB(NULL);

    dumpIVFAB(NULL);
    
    dumpEBFace(NULL);


    dumpVVoFs( NULL);

    dumpVFaces(NULL);

    dumpEBLevel(NULL);

    dumpVoFSten(NULL);

    dumpFaceSten(NULL);
  }
  
}

#endif

