
#ifndef _DATASETCLIENT_H_
#define _DATASETCLIENT_H_

//
// $Id: DatasetClient.H,v 1.2 2000-10-02 20:48:41 lijewski Exp $
//

extern "C"
{
    //
    // "format" is the printf format for numerical display ( ex:  "%5.3f" ).
    // "label"  is the first part of the window title ( ex: "myfabname" ).
    // "element" is a single element of a MultiFab (which is a fab).
    // "nvar" is the number of variables in a Real array.
    // "lodim" and "hidim" are array dimensions of size [BL_SPACEDIM]
    //
    bool ArrayView (FArrayBox* fab);
    bool ArrayViewFab (FArrayBox* fab);
    bool ArrayViewFabFormatLabel (FArrayBox*  fab,
                                  const char* format,
                                  const char* label);
    //
    // View a single fab from the multifab.
    //
    bool ArrayViewMultiFabElement (MultiFab* multifab,
                                   int       element);
    bool ArrayViewMultiFabElementFormatLabel (MultiFab*   multifab,
                                              int         element,
                                              const char* format,
                                              const char* label);
    //
    // View the entire multifab.
    //
    bool ArrayViewMultiFab (MultiFab* multifab);
    bool ArrayViewMultiFabFormatLabel (MultiFab*   multifab,
                                       const char* format,
                                       const char* label);

    bool ArrayViewReal (Real*      data,
                        const int* lodim,
                        const int* hidim);
    bool ArrayViewRealFormatLabel (Real*       data,
                                   const int*  lodim,
                                   const int*  hidim,
                                   const char* format,
                                   const char* label);
    bool ArrayViewRealNVar (Real*      data,
                            int        nvar,
                            const int* lodim,
                            const int* hidim);
    bool ArrayViewRealNVarFormatLabel (Real*       data,
                                       int         nvar,
                                       const int*  lodim,
                                       const int*  hidim,
                                       const char* format,
                                       const char* label);
#if (BL_SPACEDIM == 2)
    bool ArrayViewRealDims (Real* data,
                            int   xlo,
                            int   ylo,
                            int   xhi,
                            int   yhi);
    bool ArrayViewRealDimsFormatLabel (Real*       data,
                                       int         xlo,
                                       int         ylo,
                                       int         xhi,
                                       int         yhi,
                                       const char* format,
                                       const char* label);

    bool ArrayViewRealNVarDims (Real* data,
                                int   nvar,
                                int   xlo,
                                int   ylo,
                                int   xhi,
                                int   yhi);
    bool ArrayViewRealNVarDimsFormatLabel (Real*       data,
                                           int         nvar,
                                           int         xlo,
                                           int         ylo,
                                           int         xhi,
                                           int         yhi,
                                           const char* format,
                                           const char* label);
#else
    bool ArrayViewRealDims (Real* data,
                            int   xlo,
                            int   ylo,
                            int   zlo,
                            int   xhi,
                            int   yhi,
                            int   zhi);
    bool ArrayViewRealDimsFormatLabel (Real*       data,
                                       int         xlo,
                                       int         ylo,
                                       int         zlo,
                                       int         xhi,
                                       int         yhi,
                                       int         zhi,
                                       const char* format,
                                       const char* label);
    bool ArrayViewRealNVarDims (Real* data,
                                int   nvar,
                                int   xlo,
                                int   ylo,
                                int   zlo,
                                int   xhi,
                                int   yhi,
                                int   zhi);
    bool ArrayViewRealNVarDimsFormatLabel (Real*       data,
                                           int         nvar,
                                           int         xlo,
                                           int         ylo,
                                           int         zlo,
                                           int         xhi,
                                           int         yhi,
                                           int         zhi,
                                           const char* format,
                                           const char* label);
#endif

    bool ArrayViewRealPtrArrayNVarDims (Real*       data[],
                                        int         nvar,
                                        const int*  lodim,
                                        const int*  hidim,
                                        const char* format,
                                        const char* label);
    //
    // Stuff for TagBoxs and TagBoxArrays as well.
    //
    bool ArrayViewTagBox (TagBox* tb);
    bool ArrayViewTagBoxArray (TagBoxArray* tba);
};

//
// This class is a helper class for the ArrayView stuff.
//
// We define it's member functions in the same module as the above defined
// functions that we want accessible by the debugger.  We then define an
// object of this class in any module that includes this file.  Since the
// compiler MUST construct and destruct that object, it'll have to link in
// the module that contains those member functions, which just happens to
// contain the definitions for the above functions.
//

class ArrayViewHelperClass
{
public:
    ArrayViewHelperClass ();
    ~ArrayViewHelperClass ();
};

static ArrayViewHelperClass An_Array_View_Helper_Class_Object;

#endif /*_DATASETCLIENT_H_*/
