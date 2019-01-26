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
#warning operators.old.c represents an older, lower performance, less performance portable approach to smoothers/residual calculation.  It is strongly suggested you use the default operator file.
//------------------------------------------------------------------------------------------------------------------------------
void apply_BCs(level_type * level, int x_id, int shape){apply_BCs_p1(level,x_id,shape);}
//------------------------------------------------------------------------------------------------------------------------------
#define Dinv_ijk() Dinv[ijk]        // simply retrieve it rather than recalculating it
//------------------------------------------------------------------------------------------------------------------------------
#ifdef STENCIL_VARIABLE_COEFFICIENT
#ifdef USE_HELMHOLTZ // variable coefficient Helmholtz...
  #define apply_op_ijk(x)                               \
  (                                                     \
    a*alpha[ijk]*x[ijk]                                 \
   -b*h2inv*(                                           \
      + beta_i[ijk+1      ]*( x[ijk+1      ] - x[ijk] ) \
      + beta_i[ijk        ]*( x[ijk-1      ] - x[ijk] ) \
      + beta_j[ijk+jStride]*( x[ijk+jStride] - x[ijk] ) \
      + beta_j[ijk        ]*( x[ijk-jStride] - x[ijk] ) \
      + beta_k[ijk+kStride]*( x[ijk+kStride] - x[ijk] ) \
      + beta_k[ijk        ]*( x[ijk-kStride] - x[ijk] ) \
    )                                                   \
  )
#else // variable coefficient Poisson...
  #define apply_op_ijk(x)                               \
  (                                                     \
    -b*h2inv*(                                          \
      + beta_i[ijk+1      ]*( x[ijk+1      ] - x[ijk] ) \
      + beta_i[ijk        ]*( x[ijk-1      ] - x[ijk] ) \
      + beta_j[ijk+jStride]*( x[ijk+jStride] - x[ijk] ) \
      + beta_j[ijk        ]*( x[ijk-jStride] - x[ijk] ) \
      + beta_k[ijk+kStride]*( x[ijk+kStride] - x[ijk] ) \
      + beta_k[ijk        ]*( x[ijk-kStride] - x[ijk] ) \
    )                                                   \
  )
#endif
#else  // constant coefficient case...  
  #define apply_op_ijk(x)            \
  (                                \
    a*x[ijk] - b*h2inv*(           \
      + x[ijk+1      ]             \
      + x[ijk-1      ]             \
      + x[ijk+jStride]             \
      + x[ijk-jStride]             \
      + x[ijk+kStride]             \
      + x[ijk-kStride]             \
      - x[ijk        ]*6.0         \
    )                              \
  )
#endif // variable/constant coefficient


//------------------------------------------------------------------------------------------------------------------------------
int stencil_get_radius(){return(1);} // 7pt reaches out 1 point
int stencil_get_shape(){return(STENCIL_SHAPE_STAR);} // needs just faces
//------------------------------------------------------------------------------------------------------------------------------
void rebuild_operator(level_type * level, level_type *fromLevel, double a, double b){
  if(level->my_rank==0){fprintf(stdout,"  rebuilding operator for level...  h=%e  ",level->h);}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // form restriction of alpha[], beta_*[] coefficients from fromLevel
  if(fromLevel != NULL){
    restriction(level,VECTOR_ALPHA ,fromLevel,VECTOR_ALPHA ,RESTRICT_CELL  );
    restriction(level,VECTOR_BETA_I,fromLevel,VECTOR_BETA_I,RESTRICT_FACE_I);
    restriction(level,VECTOR_BETA_J,fromLevel,VECTOR_BETA_J,RESTRICT_FACE_J);
    restriction(level,VECTOR_BETA_K,fromLevel,VECTOR_BETA_K,RESTRICT_FACE_K);
  } // else case assumes alpha/beta have been set


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // exchange alpha/beta/...  (must be done before calculating Dinv)
  exchange_boundary(level,VECTOR_ALPHA ,STENCIL_SHAPE_BOX); // safe
  exchange_boundary(level,VECTOR_BETA_I,STENCIL_SHAPE_BOX);
  exchange_boundary(level,VECTOR_BETA_J,STENCIL_SHAPE_BOX);
  exchange_boundary(level,VECTOR_BETA_K,STENCIL_SHAPE_BOX);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // calculate Dinv, L1inv, and estimate the dominant Eigenvalue
  double _timeStart = getTime();
  int block;

  double dominant_eigenvalue = -1e9;

  PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,block,level->num_my_blocks,dominant_eigenvalue)
  for(block=0;block<level->num_my_blocks;block++){
    const int box = level->my_blocks[block].read.box;
    const int ilo = level->my_blocks[block].read.i;
    const int jlo = level->my_blocks[block].read.j;
    const int klo = level->my_blocks[block].read.k;
    const int ihi = level->my_blocks[block].dim.i + ilo;
    const int jhi = level->my_blocks[block].dim.j + jlo;
    const int khi = level->my_blocks[block].dim.k + klo;
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    double h2inv = 1.0/(level->h*level->h);
    double * __restrict__ alpha  = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_i = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_j = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_k = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride);
    double * __restrict__   Dinv = level->my_boxes[box].vectors[VECTOR_DINV  ] + ghosts*(1+jStride+kStride);
    double * __restrict__  L1inv = level->my_boxes[box].vectors[VECTOR_L1INV ] + ghosts*(1+jStride+kStride);
    double block_eigenvalue = -1e9;

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){ 
      int ijk = i + j*jStride + k*kStride;

      // used for quick linear approximation to zero dirichlet BC
      double ilo_is_valid =1.0;
      double ihi_is_valid =1.0;
      double jlo_is_valid =1.0;
      double jhi_is_valid =1.0;
      double klo_is_valid =1.0;
      double khi_is_valid =1.0;
      if(level->boundary_condition.type != BC_PERIODIC){
         if(level->my_boxes[box].low.i+i-1 <             0)ilo_is_valid = 0.0;
         if(level->my_boxes[box].low.j+j-1 <             0)jlo_is_valid = 0.0;
         if(level->my_boxes[box].low.k+k-1 <             0)klo_is_valid = 0.0;
         if(level->my_boxes[box].low.i+i+1 >= level->dim.i)ihi_is_valid = 0.0;
         if(level->my_boxes[box].low.j+j+1 >= level->dim.j)jhi_is_valid = 0.0;
         if(level->my_boxes[box].low.k+k+1 >= level->dim.k)khi_is_valid = 0.0;
       }

      #ifdef STENCIL_VARIABLE_COEFFICIENT
      // radius of Gershgorin disc is the sum of the absolute values of the off-diagonal elements...
      double sumAbsAij = fabs(b*h2inv) * (
                           fabs( beta_i[ijk        ]*ilo_is_valid )+
                           fabs( beta_j[ijk        ]*jlo_is_valid )+
                           fabs( beta_k[ijk        ]*klo_is_valid )+
                           fabs( beta_i[ijk+1      ]*ihi_is_valid )+
                           fabs( beta_j[ijk+jStride]*jhi_is_valid )+
                           fabs( beta_k[ijk+kStride]*khi_is_valid )
                         );

      // center of Gershgorin disc is the diagonal element...
      double    Aii = a*alpha[ijk] - b*h2inv*(
                        beta_i[ijk        ]*( ilo_is_valid-2.0 )+
                        beta_j[ijk        ]*( jlo_is_valid-2.0 )+
                        beta_k[ijk        ]*( klo_is_valid-2.0 )+
                        beta_i[ijk+1      ]*( ihi_is_valid-2.0 )+
                        beta_j[ijk+jStride]*( jhi_is_valid-2.0 )+
                        beta_k[ijk+kStride]*( khi_is_valid-2.0 ) 
                      );
      #else // Constant coefficient versions with fused BC's...
      // radius of Gershgorin disc is the sum of the absolute values of the off-diagonal elements...
      double sumAbsAij = fabs(b*h2inv) * (
                           ilo_is_valid +
                           jlo_is_valid +
                           klo_is_valid +
                           ihi_is_valid +
                           jhi_is_valid +
                           khi_is_valid 
                         );

      // center of Gershgorin disc is the diagonal element...
      double    Aii = a - b*h2inv*(
                         ilo_is_valid +
                         jlo_is_valid +
                         klo_is_valid +
                         ihi_is_valid +
                         jhi_is_valid +
                         khi_is_valid - 12.0
                      );
      #endif

      // calculate Dinv = D^{-1}, L1inv = ( D+D^{L1} )^{-1}, and the dominant eigenvalue...
                             Dinv[ijk] = 1.0/Aii;				// inverse of the diagonal Aii
                          //L1inv[ijk] = 1.0/(Aii+sumAbsAij);			// inverse of the L1 row norm... L1inv = ( D+D^{L1} )^{-1}
      if(Aii>=1.5*sumAbsAij)L1inv[ijk] = 1.0/(Aii              ); 		// as suggested by eq 6.5 in Baker et al, "Multigrid smoothers for ultra-parallel computing: additional theory and discussion"...
                       else L1inv[ijk] = 1.0/(Aii+0.5*sumAbsAij);		// 
      double Di = (Aii + sumAbsAij)/Aii;if(Di>block_eigenvalue)block_eigenvalue=Di;	// upper limit to Gershgorin disc == bound on dominant eigenvalue
    }}}
    if(block_eigenvalue>dominant_eigenvalue){dominant_eigenvalue = block_eigenvalue;}
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Reduce the local estimates dominant eigenvalue to a global estimate
  #ifdef USE_MPI
  double _timeStartAllReduce = getTime();
  double send = dominant_eigenvalue;
  MPI_Allreduce(&send,&dominant_eigenvalue,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  double _timeEndAllReduce = getTime();
  level->timers.collectives   += (double)(_timeEndAllReduce-_timeStartAllReduce);
  #endif
  if(level->my_rank==0){fprintf(stdout,"eigenvalue_max<%e\n",dominant_eigenvalue);}
  level->dominant_eigenvalue_of_DinvA = dominant_eigenvalue;


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // exchange Dinv/L1inv/...
  exchange_boundary(level,VECTOR_DINV ,STENCIL_SHAPE_BOX); // safe
  exchange_boundary(level,VECTOR_L1INV,STENCIL_SHAPE_BOX);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}


//------------------------------------------------------------------------------------------------------------------------------
#include "operators.old/iterators.c"
#ifdef  USE_GSRB
#define NUM_SMOOTHS      2 // RBRB
#include "operators.old/gsrb.c"
#elif   USE_CHEBY
#define NUM_SMOOTHS      1
#define CHEBYSHEV_DEGREE 4 // i.e. one degree-4 polynomial smoother
#include "operators.old/chebyshev.c"
#elif   USE_JACOBI
#define NUM_SMOOTHS      6
#include "operators.old/jacobi.c"
#elif   USE_L1JACOBI
#define NUM_SMOOTHS      6
#include "operators.old/jacobi.c"
#elif   USE_SYMGS
#define NUM_SMOOTHS      1
#include "operators.old/symgs.c"
#else
#error You must compile with either -DUSE_GSRB, -DUSE_CHEBY, -DUSE_JACOBI, -DUSE_L1JACOBI, or -DUSE_SYMGS
#endif
#include "operators.old/residual.c"
#include "operators.old/apply_op.c"
//------------------------------------------------------------------------------------------------------------------------------
#include "operators/blockCopy.c"
#include "operators/misc.c"
#include "operators/exchange_boundary.c"
#include "operators/boundary_fd.c"
#include "operators/restriction.c"
#include "operators/interpolation_p0.c"
#include "operators/interpolation_p1.c"
//------------------------------------------------------------------------------------------------------------------------------
void interpolation_vcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){interpolation_p0(level_f,id_f,prescale_f,level_c,id_c);}
void interpolation_fcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){interpolation_p1(level_f,id_f,prescale_f,level_c,id_c);}
//------------------------------------------------------------------------------------------------------------------------------
#include "operators/problem.p6.c"
//------------------------------------------------------------------------------------------------------------------------------
