//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
#ifdef _OPENMP
#include <omp.h>
#endif
//------------------------------------------------------------------------------------------------------------------------------
#include "timers.h"
#include "defines.h"
#include "level.h"
#include "operators.h"
//------------------------------------------------------------------------------------------------------------------------------
#define MyPragma(a) _Pragma(#a)
//------------------------------------------------------------------------------------------------------------------------------
#if (_OPENMP>=201107) // OpenMP 3.1 supports max reductions...
  // XL C/C++ 12.01.0000.0009 sets _OPENMP to 201107, but does not support the max clause within a _Pragma().  
  // This issue was fixed by XL C/C++ 12.01.0000.0011
  // If you do not have this version of XL C/C++ and run into this bug, uncomment these macros...
  //#warning not threading norm() calculations due to issue with XL/C, _Pragma, and reduction(max:bmax)
  //#define PRAGMA_THREAD_ACROSS_BLOCKS(    level,b,nb     )    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1)                     )
  //#define PRAGMA_THREAD_ACROSS_BLOCKS_SUM(level,b,nb,bsum)    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1) reduction(  +:bsum) )
  //#define PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,b,nb,bmax)    
  #define PRAGMA_THREAD_ACROSS_BLOCKS(    level,b,nb     )    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1)                     )
  #define PRAGMA_THREAD_ACROSS_BLOCKS_SUM(level,b,nb,bsum)    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1) reduction(  +:bsum) )
  #define PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,b,nb,bmax)    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1) reduction(max:bmax) )
#elif _OPENMP // older OpenMP versions don't support the max reduction clause
  #warning Threading max reductions requires OpenMP 3.1 (July 2011).  Please upgrade your compiler.                                                           
  #define PRAGMA_THREAD_ACROSS_BLOCKS(    level,b,nb     )    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1)                     )
  #define PRAGMA_THREAD_ACROSS_BLOCKS_SUM(level,b,nb,bsum)    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1) reduction(  +:bsum) )
  #define PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,b,nb,bmax)    
#else // flat MPI should not define any threading...
  #define PRAGMA_THREAD_ACROSS_BLOCKS(    level,b,nb     )    
  #define PRAGMA_THREAD_ACROSS_BLOCKS_SUM(level,b,nb,bsum)    
  #define PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,b,nb,bmax)    
#endif
//------------------------------------------------------------------------------------------------------------------------------
void apply_BCs(level_type * level, int x_id, int shape){apply_BCs_p2(level,x_id,shape);} // 27pt uses cell centered, not cell averaged
//void apply_BCs(level_type * level, int x_id, int shape){apply_BCs_v2(level,x_id,shape);}
//------------------------------------------------------------------------------------------------------------------------------
#define STENCIL_COEF0 (-4.2666666666666666666)  // -128.0/30.0;
#define STENCIL_COEF1 ( 0.4666666666666666666)  //   14.0/30.0;
#define STENCIL_COEF2 ( 0.1000000000000000000)  //    3.0/30.0;
#define STENCIL_COEF3 ( 0.0333333333333333333)  //    1.0/30.0;
//------------------------------------------------------------------------------------------------------------------------------
#ifdef STENCIL_VARIABLE_COEFFICIENT
  #error This implementation does not support variable-coefficient operators
#endif
#ifdef STENCIL_FUSE_BC
  #error This implementation does not support fusion of the boundary conditions with the operator
#endif
//------------------------------------------------------------------------------------------------------------------------------
#define Dinv_ijk() Dinv[ijk]        // simply retrieve it rather than recalculating it
//------------------------------------------------------------------------------------------------------------------------------
#define apply_op_ijk(x)				\
(						\
  a*x[ijk] - b*h2inv*(				\
    STENCIL_COEF3*(x[ijk-kStride-jStride-1] +	\
                   x[ijk-kStride-jStride+1] +	\
                   x[ijk-kStride+jStride-1] +	\
                   x[ijk-kStride+jStride+1] +	\
                   x[ijk+kStride-jStride-1] +	\
                   x[ijk+kStride-jStride+1] +	\
                   x[ijk+kStride+jStride-1] +	\
                   x[ijk+kStride+jStride+1] ) +	\
    STENCIL_COEF2*(x[ijk-kStride-jStride  ] +	\
                   x[ijk-kStride        -1] +	\
                   x[ijk-kStride        +1] +	\
                   x[ijk-kStride+jStride  ] +	\
                   x[ijk        -jStride-1] +	\
                   x[ijk        -jStride+1] +	\
                   x[ijk        +jStride-1] +	\
                   x[ijk        +jStride+1] +	\
                   x[ijk+kStride-jStride  ] +	\
                   x[ijk+kStride        -1] +	\
                   x[ijk+kStride        +1] +	\
                   x[ijk+kStride+jStride  ] ) +	\
    STENCIL_COEF1*(x[ijk-kStride          ] +	\
                   x[ijk        -jStride  ] +	\
                   x[ijk                -1] +	\
                   x[ijk                +1] +	\
                   x[ijk        +jStride  ] +	\
                   x[ijk+kStride          ] ) +	\
    STENCIL_COEF0*(x[ijk                  ] )	\
  )						\
)
//------------------------------------------------------------------------------------------------------------------------------
int stencil_get_radius(){return(1);} // 27pt = dense 3^3
int stencil_get_shape(){return(STENCIL_SHAPE_BOX);} // needs faces, edges, and corners
//------------------------------------------------------------------------------------------------------------------------------
void rebuild_operator(level_type * level, level_type *fromLevel, double a, double b){
  // form restriction of alpha[], beta_*[] coefficients from fromLevel
  if(fromLevel != NULL){
    restriction(level,VECTOR_ALPHA ,fromLevel,VECTOR_ALPHA ,RESTRICT_CELL  );
    restriction(level,VECTOR_BETA_I,fromLevel,VECTOR_BETA_I,RESTRICT_FACE_I);
    restriction(level,VECTOR_BETA_J,fromLevel,VECTOR_BETA_J,RESTRICT_FACE_J);
    restriction(level,VECTOR_BETA_K,fromLevel,VECTOR_BETA_K,RESTRICT_FACE_K);
  } // else case assumes alpha/beta have been set

  // exchange alpha/beta/...  (must be done before calculating Dinv)
  exchange_boundary(level,VECTOR_ALPHA ,STENCIL_SHAPE_BOX); // safe
  exchange_boundary(level,VECTOR_BETA_I,STENCIL_SHAPE_BOX);
  exchange_boundary(level,VECTOR_BETA_J,STENCIL_SHAPE_BOX);
  exchange_boundary(level,VECTOR_BETA_K,STENCIL_SHAPE_BOX);

  // black box rebuild of D^{-1}, l1^{-1}, dominant eigenvalue, ...
  rebuild_operator_blackbox(level,a,b,2);

  // exchange Dinv/L1inv/...
  exchange_boundary(level,VECTOR_DINV ,STENCIL_SHAPE_BOX); // safe
  exchange_boundary(level,VECTOR_L1INV,STENCIL_SHAPE_BOX);
}


//------------------------------------------------------------------------------------------------------------------------------
#ifdef  USE_GSRB
#warning GSRB is not recommended for the 27pt operator
#define GSRB_OOP
#define NUM_SMOOTHS      2 // RBRB
#include "operators/gsrb.c"
#elif   USE_CHEBY
#define NUM_SMOOTHS      1
#define CHEBYSHEV_DEGREE 4 // i.e. one degree-4 polynomial smoother
#include "operators/chebyshev.c"
#elif   USE_JACOBI
#define NUM_SMOOTHS      6
#include "operators/jacobi.c"
#elif   USE_L1JACOBI
#define NUM_SMOOTHS      6
#include "operators/jacobi.c"
#elif   USE_SYMGS
#define NUM_SMOOTHS      2 // FBFB
#include "operators/symgs.c"
#else
#error You must compile with either -DUSE_GSRB, -DUSE_CHEBY, -DUSE_JACOBI, -DUSE_L1JACOBI, or -DUSE_SYMGS
#endif
#include "operators/residual.c"
#include "operators/apply_op.c"
#include "operators/rebuild.c"
//------------------------------------------------------------------------------------------------------------------------------
#include "operators/blockCopy.c"
#include "operators/misc.c"
#include "operators/exchange_boundary.c"
#include "operators/boundary_fd.c" // 27pt uses cell centered, not cell averaged
//#include "operators/boundary_fv.c"
#include "operators/restriction.c"
#include "operators/interpolation_p2.c"
//#include "operators/interpolation_v2.c"
//------------------------------------------------------------------------------------------------------------------------------
void interpolation_vcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){interpolation_p2(level_f,id_f,prescale_f,level_c,id_c);} // 27pt uses cell centered, not cell averaged
void interpolation_fcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){interpolation_p2(level_f,id_f,prescale_f,level_c,id_c);}
//void interpolation_vcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){interpolation_v2(level_f,id_f,prescale_f,level_c,id_c);}
//void interpolation_fcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){interpolation_v2(level_f,id_f,prescale_f,level_c,id_c);}
//------------------------------------------------------------------------------------------------------------------------------
#include "operators/problem.p6.c"
//------------------------------------------------------------------------------------------------------------------------------
