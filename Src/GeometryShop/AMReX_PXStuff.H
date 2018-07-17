#ifndef AMREX_PXSTUFF_H_
#define AMREX_PXSTUFF_H_

namespace amrex
{
  enum PXE_Boolean { PXE_False, PXE_True, PXE_BooleanLast};


#ifndef unlikely
#define unlikely(expr) __builtin_expect((expr), 0)
#endif

#ifndef likely
#define likely(expr) __builtin_expect((expr), 1)
#endif


//#ifdef __cplusplus
//#define RESTRICT __restrict__
//#else
//#define RESTRICT restrict
//#endif

  extern void PXErrorReport( const char *file, int line, const char *call, int ierr);


#define PXError(v,X)             { int _ier; _ier = (X); if (unlikely(_ier != PX_NO_ERROR)){ PXErrorReport(__FILE__, __LINE__, #X, _ier); } (v) = _ier;}
#define PXErrorBreak(X)        do{ int _ier; _ier = (X); if (unlikely(_ier != PX_NO_ERROR)){ PXErrorReport(__FILE__, __LINE__, #X, _ier);    exit(-1);} } while (0)
#define PXErrorReturn(X)       do{ int _ier; _ier = (X); if (unlikely(_ier != PX_NO_ERROR)){ PXErrorReport(__FILE__, __LINE__, #X, _ier); return _ier;} } while (0)
#define PXErrorReturnCode(X,Y) do{ int _ier; _ier = (X); if (unlikely(_ier != PX_NO_ERROR)){ PXErrorReport(__FILE__, __LINE__, #X, _ier); return ( Y);} } while (0)
#define PXErrorReturnSilent(X) do{ int _ier; _ier = (X); if (unlikely(_ier != PX_NO_ERROR)){ return _ier;} } while (0)


#define PX_NO_ERROR 0
#define PX_MEMORY_ERROR -1
#define PX_BAD_INPUT -2
#define PX_NON_PHYSICAL -3
#define PX_READWRITE_ERROR -4
#define PX_GRID_ERROR -5
#define PX_SEARCH_NOT_FOUND -6
#define PX_NO_UPDATE -7
#define PX_PARALLEL_ERROR -8
#define PX_CODE_FLOW_ERROR -9
#define PX_SYSTEM_ERROR -10
#define PX_DYNAMIC_LIBRARY_ERROR -11
#define PX_NOT_CONVERGED -12
#define PX_VIZ_ERROR -13
#define PX_LAPACK_ERROR -14
#define PX_HARD_EXIT -15
#define PX_CGNS_ERROR -16

  extern const char *PXE_ErrorCodeName[];

#define AssertDblWithin( truth, actual, delta)                          \
  if (Abs((truth) - (actual)) > (delta)) PXErrorReturn(PX_CODE_FLOW_ERROR);

#define PXAssertDblEquals( truth, actual )                              \
  CuAssertDblEquals( private_test_struct, truth, actual, 1.0e-15)

#define PXAssertDblRelErrWithin( truth, actual, delta)                  \
  CuAssertDblRelErrWithin( private_test_struct, truth, actual, delta )

#define PXAssertIntWithin( truth, actual, delta )                       \
  if ( abs(truth - actual) > delta){                                    \
    sprintf(private_test_message,"expected <%d> but was <%d>", truth, actual ); \
    CuFail( private_test_struct, private_test_message); }

#define AssertIntEquals( truth, actual )                        \
  if ((truth) != (actual)) PXErrorReturn(PX_CODE_FLOW_ERROR);

#define PXAssertPtrNotNull( ptr )                       \
  CuAssertPtrNotNull( private_test_struct, ptr )

#define PXAssertPtrEquals( truth, actual )                      \
  CuAssertPtrEquals( private_test_struct, truth, actual )

#define PXAssertStrEquals( truth, actual )                      \
  CuAssertStrEquals( private_test_struct, truth, actual )

#define AssertDblGreater( min, actual )                         \
  if ( !((actual) > (min)) ) PXErrorReturn(PX_CODE_FLOW_ERROR);

#define PXAssertIsNotNaN( value )                               \
  if ( !(value==value) ){                                       \
    sprintf(private_test_message,                               \
            "NaN encountered.\n");                              \
    CuFail( private_test_struct, private_test_message ); }

#define ErrorAssert(X) do{ int _ier; _ier = (X); if (_ier != PX_NO_ERROR){ printf("Error %d (%s) has o\
ccured.\n File : %s  Line : %d\n Call : %s\n", _ier,PXE_ErrorCodeName[-_ier], __FILE__, __LINE__, #X); f \
                                                                                                         flush(stdout); AssertIntEquals(PX_NO_ERROR, _ier);} } while (0)


}

#endif
