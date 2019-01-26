//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int x_id, int rhs_id, double a, double b){
  if(NUM_SMOOTHS&1){
    fprintf(stderr,"error - NUM_SMOOTHS must be even...\n");
    exit(0);
  }

  #ifdef USE_L1JACOBI
  double weight = 1.0;
  #else
  double weight = 2.0/3.0;
  #endif
 
  int box,s;
  for(s=0;s<NUM_SMOOTHS;s++){
    // exchange ghost zone data... Jacobi ping pongs between x_id and VECTOR_TEMP
    if((s&1)==0){exchange_boundary(level,       x_id,stencil_get_shape());apply_BCs(level,       x_id,stencil_get_shape());}
            else{exchange_boundary(level,VECTOR_TEMP,stencil_get_shape());apply_BCs(level,VECTOR_TEMP,stencil_get_shape());}

    // apply the smoother... Jacobi ping pongs between x_id and VECTOR_TEMP
    double _timeStart = getTime();
    const int  ghosts = level->box_ghosts;
    const int jStride = level->box_jStride;
    const int kStride = level->box_kStride;
    const int     dim = level->box_dim;
    const double h2inv = 1.0/(level->h*level->h);

    PRAGMA_THREAD_ACROSS_BOXES(level,box)
    for(box=0;box<level->num_my_boxes;box++){
      int i,j,k;
      const double * __restrict__ rhs    = level->my_boxes[box].vectors[       rhs_id] + ghosts*(1+jStride+kStride);
      const double * __restrict__ alpha  = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_i = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_j = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_k = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride);
      #ifdef USE_L1JACOBI
      const double * __restrict__ lambda = level->my_boxes[box].vectors[VECTOR_L1INV ] + ghosts*(1+jStride+kStride);
      #else
      const double * __restrict__ lambda = level->my_boxes[box].vectors[VECTOR_DINV  ] + ghosts*(1+jStride+kStride);
      #endif
        const double * __restrict__ x_n;
              double * __restrict__ x_np1;
                      if((s&1)==0){x_n   = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride);
                                   x_np1 = level->my_boxes[box].vectors[VECTOR_TEMP  ] + ghosts*(1+jStride+kStride);}
                              else{x_n   = level->my_boxes[box].vectors[VECTOR_TEMP  ] + ghosts*(1+jStride+kStride);
                                   x_np1 = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride);}
      PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
      for(k=0;k<dim;k++){
      for(j=0;j<dim;j++){
      for(i=0;i<dim;i++){
        int ijk = i + j*jStride + k*kStride;
        double Ax_n = apply_op_ijk(x_n);
        x_np1[ijk] = x_n[ijk] + weight*lambda[ijk]*(rhs[ijk]-Ax_n);
      }}}
    } // box-loop
    level->timers.smooth += (double)(getTime()-_timeStart);
  } // s-loop
}

//------------------------------------------------------------------------------------------------------------------------------
