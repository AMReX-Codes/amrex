// -------------------------------------------------------------
// ArrayView.H
// -------------------------------------------------------------
#ifndef _ARRAYVIEW_H_
#define _ARRAYVIEW_H_

#include "FArrayBox.H"
#include "MultiFab.H"
#ifdef BL_ARRAYVIEW_TAGBOX
#include "TagBox.H"
#endif
#include "DatasetClient.H"

// -------------------------------------------------------------
void NeverCalled() {
  ArrayView(NULL);
  ArrayViewFab(NULL);
  ArrayViewFabFormatLabel(NULL, "%7.5f", "FabLabel");

  ArrayViewMultiFabElement(NULL, 0);
  ArrayViewMultiFabElementFormatLabel(NULL, 0, "%7.5f", "MultiFabLabel");

  ArrayViewMultiFab(NULL);
  ArrayViewMultiFabFormatLabel(NULL, "%7.5f", "MultiFabLabel");

  ArrayViewReal(NULL, NULL, NULL);
  ArrayViewRealFormatLabel(NULL, NULL, NULL, "%7.5f", "RealLabel");

  ArrayViewRealNVar(NULL, 0, NULL, NULL);
  ArrayViewRealNVarFormatLabel(NULL, 0, NULL, NULL, "%7.5f", "RealLabel");

#if (BL_SPACEDIM == 2)
  ArrayViewRealDims(NULL, 0, 0, 1, 1);
  ArrayViewRealDimsFormatLabel(NULL, 0, 0, 1, 1, "%7.5f", "RealLabel");

  ArrayViewRealNVarDims(NULL, 1, 0, 0, 1, 1);
  ArrayViewRealNVarDimsFormatLabel(NULL, 1, 0, 0, 1, 1, "%7.5f", "RealLabel");
#else
  ArrayViewRealDims(NULL, 0, 0, 0, 1, 1, 1);
  ArrayViewRealDimsFormatLabel(NULL, 0, 0, 0, 1, 1, 1, "%7.5f", "RealLabel");
  ArrayViewRealNVarDims(NULL, 1, 0, 0, 0, 1, 1, 1);
  ArrayViewRealNVarDimsFormatLabel(NULL, 1, 0, 0, 0, 1, 1, 1, "%7.5f", "RealLabel");
#endif

  ArrayViewRealPtrArrayNVarDims(NULL, 1, NULL, NULL, "%7.5f", "RealLabel");

#ifdef BL_ARRAYVIEW_TAGBOX
  ArrayViewTagBox(NULL);
  ArrayViewTagBoxArray(NULL);
#endif

}
// -------------------------------------------------------------
// -------------------------------------------------------------
#endif
