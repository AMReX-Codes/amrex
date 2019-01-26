//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#if   defined(GSRB_FP)
  #warning Overriding default GSRB implementation and using pre-computed 1.0/0.0 FP array for Red-Black to facilitate vectorization...
#elif defined(GSRB_STRIDE2)
  #if defined(GSRB_OOP)
  #warning Overriding default GSRB implementation and using out-of-place and stride-2 accesses to minimize the number of flops
  #else
  #warning Overriding default GSRB implementation and using stride-2 accesses to minimize the number of flops
  #endif
#elif defined(GSRB_BRANCH)
  #if defined(GSRB_OOP)
  #warning Overriding default GSRB implementation and using out-of-place implementation with an if-then-else on loop indices...
  #else
  #warning Overriding default GSRB implementation and using if-then-else on loop indices...
  #endif
#else
#define GSRB_STRIDE2 // default implementation
#endif
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int phi_id, int rhs_id, double a, double b){
  int box,s;
  for(s=0;s<2*NUM_SMOOTHS;s++){ // there are two sweeps per GSRB smooth

    // exchange the ghost zone...
    #ifdef GSRB_OOP // out-of-place GSRB ping pongs between x and VECTOR_TEMP
    if((s&1)==0){exchange_boundary(level,     phi_id,stencil_get_shape());apply_BCs(level,     phi_id,stencil_get_shape());}
            else{exchange_boundary(level,VECTOR_TEMP,stencil_get_shape());apply_BCs(level,VECTOR_TEMP,stencil_get_shape());}
    #else // in-place GSRB only operates on x
                 exchange_boundary(level,     phi_id,stencil_get_shape());apply_BCs(level,     phi_id,stencil_get_shape());
    #endif


    // apply the smoother...
    double _timeStart = getTime();
    const int  ghosts = level->box_ghosts;
    const int jStride = level->box_jStride;
    const int kStride = level->box_kStride;
    const int     dim = level->box_dim;
    const double h2inv = 1.0/(level->h*level->h);

    PRAGMA_THREAD_ACROSS_BOXES(level,box)
    for(box=0;box<level->num_my_boxes;box++){
      int i,j,k;
      const int color000 = (level->my_boxes[box].low.i^level->my_boxes[box].low.j^level->my_boxes[box].low.k^s)&1;  // is element 000 red or black on *THIS* sweep

      const double * __restrict__ rhs      = level->my_boxes[box].vectors[       rhs_id] + ghosts*(1+jStride+kStride);
      const double * __restrict__ alpha    = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_i   = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_j   = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_k   = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride);
      const double * __restrict__ Dinv     = level->my_boxes[box].vectors[VECTOR_DINV  ] + ghosts*(1+jStride+kStride);
      #ifdef GSRB_OOP
      const double * __restrict__ phi;
            double * __restrict__ phi_new;
                     if((s&1)==0){phi      = level->my_boxes[box].vectors[       phi_id] + ghosts*(1+jStride+kStride);
                                  phi_new  = level->my_boxes[box].vectors[VECTOR_TEMP  ] + ghosts*(1+jStride+kStride);}
                             else{phi      = level->my_boxes[box].vectors[VECTOR_TEMP  ] + ghosts*(1+jStride+kStride);
                                  phi_new  = level->my_boxes[box].vectors[       phi_id] + ghosts*(1+jStride+kStride);}
      #else
      const double * __restrict__ phi      = level->my_boxes[box].vectors[       phi_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
            double * __restrict__ phi_new  = level->my_boxes[box].vectors[       phi_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
      #endif
          

      #if defined(GSRB_FP)
      PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
      for(k=0;k<dim;k++){
      for(j=0;j<dim;j++){
      const double * __restrict__ RedBlack = level->RedBlack_FP + ghosts*(1+jStride) + kStride*((k^color000)&0x1);
      for(i=0;i<dim;i++){
        int ij  = i + j*jStride;
        int ijk = i + j*jStride + k*kStride;
        double Ax     = apply_op_ijk(phi);
        double lambda =     Dinv_ijk();
        phi_new[ijk] = phi[ijk] + RedBlack[ij]*lambda*(rhs[ijk]-Ax);
        //phi_new[ijk] = ((i^j^k^color000)&1) ? phi[ijk] : phi[ijk] + lambda*(rhs[ijk]-Ax);
      }}}


      #elif defined(GSRB_STRIDE2)
      PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
      for(k=0;k<dim;k++){
      for(j=0;j<dim;j++){
        #ifdef GSRB_OOP
        // out-of-place must copy old value...
        for(i=0;i<dim;i++){
          int ijk = i + j*jStride + k*kStride; 
          phi_new[ijk] = phi[ijk];
        }
        #endif
        for(i=((j^k^color000)&1);i<dim;i+=2){ // stride-2 GSRB
          int ijk = i + j*jStride + k*kStride; 
          double Ax     = apply_op_ijk(phi);
          double lambda =     Dinv_ijk();
          phi_new[ijk] = phi[ijk] + lambda*(rhs[ijk]-Ax);
        }
      }}


      #elif defined(GSRB_BRANCH)
      PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
      for(k=0;k<dim;k++){
      for(j=0;j<dim;j++){
      for(i=0;i<dim;i++){
        int ijk = i + j*jStride + k*kStride;
        if((i^j^k^color000^1)&1){ // looks very clean when [0] is i,j,k=0,0,0 
          double Ax     = apply_op_ijk(phi);
          double lambda =     Dinv_ijk();
          phi_new[ijk] = phi[ijk] + lambda*(rhs[ijk]-Ax);
        #ifdef GSRB_OOP
        }else{
          phi_new[ijk] = phi[ijk]; // copy old value when sweep color != cell color
        #endif
        }
      }}}


      #else
      #error no GSRB implementation was specified
      #endif


    } // boxes
    level->timers.smooth += (double)(getTime()-_timeStart);
  } // s-loop
}


//------------------------------------------------------------------------------------------------------------------------------
