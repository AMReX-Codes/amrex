//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// Based on Yousef Saad's Iterative Methods for Sparse Linear Algebra, Algorithm 12.1, page 399
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int x_id, int rhs_id, double a, double b){
  if( (level->dominant_eigenvalue_of_DinvA<=0.0) && (level->my_rank==0) )printf("dominant_eigenvalue_of_DinvA <= 0.0 !\n");
  if((CHEBYSHEV_DEGREE*NUM_SMOOTHS)&1){
    printf("error... CHEBYSHEV_DEGREE*NUM_SMOOTHS must be even for the chebyshev smoother...\n");
    exit(0);
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int box,s;
  int ghosts = level->box_ghosts;
  int communicationAvoiding = ghosts > stencil_get_radius(); 


  // compute the Chebyshev coefficients...
  double beta     = 1.000*level->dominant_eigenvalue_of_DinvA;
//double alpha    = 0.300000*beta;
//double alpha    = 0.250000*beta;
//double alpha    = 0.166666*beta;
  double alpha    = 0.125000*beta;
  double theta    = 0.5*(beta+alpha);		// center of the spectral ellipse
  double delta    = 0.5*(beta-alpha);		// major axis?
  double sigma = theta/delta;
  double rho_n = 1/sigma;			// rho_0
  double chebyshev_c1[CHEBYSHEV_DEGREE];	// + c1*(x_n-x_nm1) == rho_n*rho_nm1
  double chebyshev_c2[CHEBYSHEV_DEGREE];	// + c2*(b-Ax_n)
  chebyshev_c1[0] = 0.0;
  chebyshev_c2[0] = 1/theta;
  for(s=1;s<CHEBYSHEV_DEGREE;s++){
    double rho_nm1 = rho_n;
    rho_n = 1.0/(2.0*sigma - rho_nm1);
    chebyshev_c1[s] = rho_n*rho_nm1;
    chebyshev_c2[s] = rho_n*2.0/delta;
  }


  // if communication-avoiding, need updated RHS for stencils in ghost zones
  if(communicationAvoiding)exchange_boundary(level,rhs_id,0); 

  for(s=0;s<CHEBYSHEV_DEGREE*NUM_SMOOTHS;s+=ghosts){
    // Chebyshev ping pongs between x_id and VECTOR_TEMP
    if((s&1)==0){exchange_boundary(level,       x_id,stencil_is_star_shaped() && !communicationAvoiding);apply_BCs(level,       x_id);}
            else{exchange_boundary(level,VECTOR_TEMP,stencil_is_star_shaped() && !communicationAvoiding);apply_BCs(level,VECTOR_TEMP);}
    
    // now do ghosts communication-avoiding smooths on each box...
    uint64_t _timeStart = CycleTime();
    PRAGMA_THREAD_ACROSS_BOXES(level,box)
    for(box=0;box<level->num_my_boxes;box++){
      int i,j,k,ss;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int     dim = level->my_boxes[box].dim;
      const double h2inv = 1.0/(level->h*level->h);
      const double * __restrict__ rhs      = level->my_boxes[box].vectors[       rhs_id] + ghosts*(1+jStride+kStride);
      const double * __restrict__ alpha    = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_i   = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_j   = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_k   = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride);
      const double * __restrict__ Dinv     = level->my_boxes[box].vectors[VECTOR_DINV  ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ valid    = level->my_boxes[box].vectors[VECTOR_VALID ] + ghosts*(1+jStride+kStride); // cell is inside the domain

      int ghostsToOperateOn=ghosts-1;
      for(ss=s;ss<s+ghosts;ss++,ghostsToOperateOn--){
              double * __restrict__ x_np1;
        const double * __restrict__ x_n;
        const double * __restrict__ x_nm1;
              if((ss&1)==0){x_n    = level->my_boxes[box].vectors[       x_id] + ghosts*(1+jStride+kStride);
                            x_nm1  = level->my_boxes[box].vectors[VECTOR_TEMP] + ghosts*(1+jStride+kStride); 
                            x_np1  = level->my_boxes[box].vectors[VECTOR_TEMP] + ghosts*(1+jStride+kStride);}
                       else{x_n    = level->my_boxes[box].vectors[VECTOR_TEMP] + ghosts*(1+jStride+kStride);
                            x_nm1  = level->my_boxes[box].vectors[       x_id] + ghosts*(1+jStride+kStride); 
                            x_np1  = level->my_boxes[box].vectors[       x_id] + ghosts*(1+jStride+kStride);}
        const double c1 = chebyshev_c1[ss%CHEBYSHEV_DEGREE]; // limit polynomial to degree CHEBYSHEV_DEGREE.
        const double c2 = chebyshev_c2[ss%CHEBYSHEV_DEGREE]; // limit polynomial to degree CHEBYSHEV_DEGREE.
        PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
        for(k=0-ghostsToOperateOn;k<dim+ghostsToOperateOn;k++){
        for(j=0-ghostsToOperateOn;j<dim+ghostsToOperateOn;j++){
        for(i=0-ghostsToOperateOn;i<dim+ghostsToOperateOn;i++){
          int ijk = i + j*jStride + k*kStride;
          // According to Saad... but his was missing a Dinv[ijk] == D^{-1} !!!
          //  x_{n+1} = x_{n} + rho_{n} [ rho_{n-1}(x_{n} - x_{n-1}) + (2/delta)(b-Ax_{n}) ]
          //  x_temp[ijk] = x_n[ijk] + c1*(x_n[ijk]-x_temp[ijk]) + c2*Dinv[ijk]*(rhs[ijk]-Ax_n);
          double Ax_n   = apply_op_ijk(x_n);
          double lambda =     Dinv_ijk();
          x_np1[ijk] = x_n[ijk] + c1*(x_n[ijk]-x_nm1[ijk]) + c2*lambda*(rhs[ijk]-Ax_n);
        }}}
      } // ss-loop
    } // box-loop
    level->cycles.smooth += (uint64_t)(CycleTime()-_timeStart);
  } // s-loop
}
