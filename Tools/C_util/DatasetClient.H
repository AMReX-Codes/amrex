// -------------------------------------------------------------
// DatasetClient.H
// -------------------------------------------------------------
#ifndef _DATASETCLIENT_H_
#define _DATASETCLIENT_H_

#include "REAL.H"

class FArrayBox;
class MultiFab;

#ifdef BL_ARRAYVIEW_TAGBOX
class TagBox;
class TagBoxArray;
#endif

extern "C" {

// -------------------------------------------------------------
// "format" is the printf format for numerical display ( ex:  "%5.3f" ).
// "label"  is the first part of the window title ( ex: "myfabname" ).
// "element" is a single element of a MultiFab (which is a fab).
// "nvar" is the number of variables in a Real array.
// "lodim" and "hidim" are array dimensions of size [BL_SPACEDIM]
// -------------------------------------------------------------
// -------------------------------------------------------------
// pointer to fab interface

  bool ArrayView(FArrayBox *fab);
  bool ArrayViewFab(FArrayBox *fab);
  bool ArrayViewFabFormatLabel(FArrayBox *fab, const char *format,
                               const char *label);

// -------------------------------------------------------------
// pointer to MultiFab interface

  // view a single fab from the multifab
  bool ArrayViewMultiFabElement(MultiFab *multifab, int element);
  bool ArrayViewMultiFabElementFormatLabel(MultiFab *multifab, int element,
                            const char *format,  const char *label);

  // view the entire multifab
  bool ArrayViewMultiFab(MultiFab *multifab);
  bool ArrayViewMultiFabFormatLabel(MultiFab *multifab, const char *format,
                                    const char *label);

// -------------------------------------------------------------
// pointer to real interface

  bool ArrayViewReal(Real *data, const int *lodim, const int *hidim);
  bool ArrayViewRealFormatLabel(Real *data, const int *lodim, const int *hidim,
                           const char *format, const char *label);
  bool ArrayViewRealNVar(Real *data, int nvar, const int *lodim, const int *hidim);
  bool ArrayViewRealNVarFormatLabel(Real *data, int nvar,
                               const int *lodim, const int *hidim,
                               const char *format, const char *label);

#if (BL_SPACEDIM == 2)
  bool ArrayViewRealDims(Real *data, int xlo, int ylo, int xhi, int yhi);
  bool ArrayViewRealDimsFormatLabel(Real *data,
                         int xlo, int ylo, int xhi, int yhi,
                         const char *format, const char *label);

  bool ArrayViewRealNVarDims(Real *data, int nvar,
                         int xlo, int ylo, int xhi, int yhi);
  bool ArrayViewRealNVarDimsFormatLabel(Real *data, int nvar,
                         int xlo, int ylo, int xhi, int yhi,
                         const char *format, const char *label);
#else
  bool ArrayViewRealDims(Real *data, int xlo, int ylo, int zlo,
                         int xhi, int yhi, int zhi);
  bool ArrayViewRealDimsFormatLabel(Real *data, int xlo, int ylo, int zlo,
                         int xhi, int yhi, int zhi,
                         const char *format, const char *label);
  bool ArrayViewRealNVarDims(Real *data, int nvar, int xlo, int ylo, int zlo,
                             int xhi, int yhi, int zhi);
  bool ArrayViewRealNVarDimsFormatLabel(Real *data, int nvar,
                             int xlo, int ylo, int zlo,
                             int xhi, int yhi, int zhi,
                             const char *format, const char *label);
#endif


  bool ArrayViewRealPtrArrayNVarDims(Real *data[], int nvar,  // size nvar
                   const int *lodim, const int *hidim,        // size BL_SPACEDIM
                   const char *format, const char *label);

#ifdef BL_ARRAYVIEW_TAGBOX
// -------------------------------------------------------------
// pointer to TagBox interface

  bool ArrayViewTagBox(TagBox *tagbox);

// -------------------------------------------------------------
// pointer to TagBoxArray interface

  bool ArrayViewTagBoxArray(TagBoxArray *tagboxarray);
#endif

};
// -------------------------------------------------------------
// -------------------------------------------------------------
#endif
