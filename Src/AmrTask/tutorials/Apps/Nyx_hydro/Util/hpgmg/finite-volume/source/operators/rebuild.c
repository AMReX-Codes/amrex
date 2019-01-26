//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
/*
// power method for calculating the dominant eigenvalue of D^{-1}A
double power_method(level_type * level, double a, double b, int max_iterations){
  int i;
  int  x_id = VECTOR_U;
  int Ax_id = VECTOR_TEMP;
  double lambda_max = 0;

  #ifdef USE_MPI
  double lmax_start = MPI_Wtime();
  #endif
  if(level->my_rank==0){fprintf(stdout,"  calculating lambda_max...");fflush(stdout);}

  random_vector(level,x_id);
  for(i=0;i<max_iterations;i++){
   apply_op(level,Ax_id, x_id,a,b);
   mul_vectors(level,Ax_id,1.0,VECTOR_DINV,Ax_id); // D^{-1}Ax
   double   x_dot_x = dot(level, x_id,x_id);
   double DAx_dot_x = dot(level,Ax_id,x_id);
   lambda_max = DAx_dot_x / x_dot_x;
   double Ax_max = norm(level,Ax_id); // renormalize Ax (== new x)
   scale_vector(level,x_id,1.0/Ax_max,Ax_id); 
  }
  #ifdef USE_MPI
  if(level->my_rank==0){fprintf(stdout,"  %1.15e (%0.6f seconds)\n",lambda_max,MPI_Wtime()-lmax_start);}
  #else
  if(level->my_rank==0){fprintf(stdout,"  %1.15e\n",lambda_max);}
  #endif
  return(lambda_max);
}
*/


//------------------------------------------------------------------------------------------------------------------------------
// Accurate estimates of D^{-1} are essential in realizing high-performance and stable smoothers.
// Unfortunately, complex boundary conditions can make it difficult to express D^{-1} analytically
// As such, this black-box routine will calculate D^{-1}, l1 norm, the dominant eigenvalue using only the apply_op_ijk macro
// colors_in_each_dim should be sufficiently large as to decouple the boundary condition from the operator
// e.g. with quartic BC's, colors_in_each_dim==4 (total of 64 colors in 3D)
// If using periodic BCs, one should be able to set colors_in_each_dim to stencil_get_radius();
// NOTE, as this function is not timed, it has not been optimized for performance.
void rebuild_operator_blackbox(level_type * level, double a, double b, int colors_in_each_dim){

  // trying to color a 1^3 grid with 8 colors won't work... reduce the number of colors...
  if(level->dim.i<colors_in_each_dim)colors_in_each_dim=level->dim.i;
  if(level->dim.j<colors_in_each_dim)colors_in_each_dim=level->dim.j;
  if(level->dim.k<colors_in_each_dim)colors_in_each_dim=level->dim.k;

  if(level->my_rank==0){fprintf(stdout,"  calculating D^{-1} exactly for level h=%e using %d colors...  ",level->h,colors_in_each_dim*colors_in_each_dim*colors_in_each_dim);fflush(stdout);}
  #ifdef USE_MPI
  double dinv_start = MPI_Wtime();
  #endif

  #if 0 // naive version using existing routines.  Doesn't calculate l1inv or estimate the dominant eigenvalue
  int         x_id = VECTOR_U;
  int        Ax_id = VECTOR_TEMP;
  int icolor,jcolor,kcolor;
  zero_vector(level,VECTOR_DINV);
  zero_vector(level,VECTOR_L1INV);
  for(kcolor=0;kcolor<colors_in_each_dim;kcolor++){
  for(jcolor=0;jcolor<colors_in_each_dim;jcolor++){
  for(icolor=0;icolor<colors_in_each_dim;icolor++){
    color_vector(level,x_id,colors_in_each_dim,icolor,jcolor,kcolor);  // color the grid as 1's and 0's
        apply_op(level,Ax_id,x_id,a,b);                                // includes effects of boundary conditions on Aii
     mul_vectors(level,Ax_id,1.0,x_id,Ax_id);                          // zero out the off-diagonal contributions 
     add_vectors(level,VECTOR_DINV,1.0,Ax_id,1.0,VECTOR_DINV);         // add to running sum of Aii
  }}}
  invert_vector(level,VECTOR_DINV,1.0,VECTOR_DINV);
  #else

  int         x_id = VECTOR_TEMP;
  int       Aii_id = VECTOR_DINV;
  int sumAbsAij_id = VECTOR_L1INV;
  int icolor,jcolor,kcolor;
  double dominant_eigenvalue = -1e9;
  int block;

  // initialize Aii[] = subAbsAij[] = 0's
  zero_vector(level,      Aii_id);
  zero_vector(level,sumAbsAij_id);

  // loop over all colors...
  for(kcolor=0;kcolor<colors_in_each_dim;kcolor++){
  for(jcolor=0;jcolor<colors_in_each_dim;jcolor++){
  for(icolor=0;icolor<colors_in_each_dim;icolor++){
    // color the grid as 1's and 0's
    color_vector(level,x_id,colors_in_each_dim,icolor,jcolor,kcolor);

    // exchange the boundary of x in preparation for Ax
    exchange_boundary(level,x_id,stencil_get_shape());
            apply_BCs(level,x_id,stencil_get_shape());
 
    // apply the operator and add to Aii and AbsAij 
    PRAGMA_THREAD_ACROSS_BLOCKS(level,block,level->num_my_blocks)
    for(block=0;block<level->num_my_blocks;block++){
      const int box = level->my_blocks[block].read.box;
      const int ilo = level->my_blocks[block].read.i;
      const int jlo = level->my_blocks[block].read.j;
      const int klo = level->my_blocks[block].read.k;
      const int ihi = level->my_blocks[block].dim.i + ilo;
      const int jhi = level->my_blocks[block].dim.j + jlo;
      const int khi = level->my_blocks[block].dim.k + klo;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int  ghosts = level->my_boxes[box].ghosts;
      const double h2inv = 1.0/(level->h*level->h);
      const double * __restrict__         x = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
      const double * __restrict__     alpha = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride);
      const double * __restrict__    beta_i = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride);
      const double * __restrict__    beta_j = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride);
      const double * __restrict__    beta_k = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride);
            double * __restrict__       Aii = level->my_boxes[box].vectors[       Aii_id] + ghosts*(1+jStride+kStride);
            double * __restrict__ sumAbsAij = level->my_boxes[box].vectors[ sumAbsAij_id] + ghosts*(1+jStride+kStride);
  
      int i,j,k;
      for(k=klo;k<khi;k++){
      for(j=jlo;j<jhi;j++){
      for(i=ilo;i<ihi;i++){
        int ijk = i + j*jStride + k*kStride;
        double Ax = apply_op_ijk(x);
              Aii[ijk] +=      (    x[ijk])*Ax; // add the effect of setting one grid point (i) to 1.0 to Aii
        sumAbsAij[ijk] += fabs((1.0-x[ijk])*Ax);
      }}}
    }
  }}}


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // take Aii and the row sum sumAbsAij and calculate D^{-1} and L1^{-1}...
  PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,block,level->num_my_blocks,dominant_eigenvalue)
  for(block=0;block<level->num_my_blocks;block++){
    const int box = level->my_blocks[block].read.box;
    const int ilo = level->my_blocks[block].read.i;
    const int jlo = level->my_blocks[block].read.j;
    const int klo = level->my_blocks[block].read.k;
    const int ihi = level->my_blocks[block].dim.i + ilo;
    const int jhi = level->my_blocks[block].dim.j + jlo;
    const int khi = level->my_blocks[block].dim.k + klo;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const double h2inv = 1.0/(level->h*level->h);
    double * __restrict__       Aii = level->my_boxes[box].vectors[      Aii_id] + ghosts*(1+jStride+kStride);
    double * __restrict__ sumAbsAij = level->my_boxes[box].vectors[sumAbsAij_id] + ghosts*(1+jStride+kStride);

    double block_eigenvalue = -1e9;
    int i,j,k;
    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
      int ijk = i + j*jStride + k*kStride;

      // catch failure...
      if(Aii[ijk]==0.0){
        printf("Aii[%d,%d,%d]==0.0 !!!\n",i+level->my_boxes[box].low.i,j+level->my_boxes[box].low.j,k+level->my_boxes[box].low.k);
        Aii[ijk] = a+b*h2inv; // FIX !!!
      }

      // upper limit to Gershgorin disc == bound on dominant eigenvalue
      double Di = (Aii[ijk] + sumAbsAij[ijk])/Aii[ijk];if(Di>block_eigenvalue)block_eigenvalue=Di;

      // inverse of the L1 row norm... L1inv = ( D+D^{L1} )^{-1}
      // sumAbsAij[ijk] = 1.0/(Aii[ijk]+sumAbsAij[ijk]);
      // alternately, as suggested by eq 6.5 in Baker et al, "Multigrid smoothers for ultra-parallel computing: additional theory and discussion"...
      if(Aii[ijk]>=1.5*sumAbsAij[ijk])sumAbsAij[ijk] = 1.0/(Aii[ijk]                   ); // VECTOR_L1INV = ...
                                 else sumAbsAij[ijk] = 1.0/(Aii[ijk]+0.5*sumAbsAij[ijk]); // VECTOR_L1INV = ...

      // inverse of the diagonal...
      Aii[ijk] = 1.0/Aii[ijk]; // VECTOR_DINV = ...

    }}}
    if(block_eigenvalue>dominant_eigenvalue){dominant_eigenvalue = block_eigenvalue;}
  }
  #ifdef USE_MPI
  if(level->my_rank==0){fprintf(stdout,"done (%0.6f seconds)\n",MPI_Wtime()-dinv_start);}
  #else
  if(level->my_rank==0){fprintf(stdout,"done\n");}
  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Reduce the local estimate of the dominant eigenvalue to a global estimate
  #ifdef USE_MPI
  double _timeStartAllReduce = getTime();
  double send = dominant_eigenvalue;
  MPI_Allreduce(&send,&dominant_eigenvalue,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  double _timeEndAllReduce = getTime();
  level->timers.collectives   += (double)(_timeEndAllReduce-_timeStartAllReduce);
  #endif
  if(level->my_rank==0){fprintf(stdout,"  estimating  lambda_max... <%1.15e\n",dominant_eigenvalue);fflush(stdout);}
  level->dominant_eigenvalue_of_DinvA = dominant_eigenvalue;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //level->dominant_eigenvalue_of_DinvA = power_method(level,a,b,10);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #endif
}
//------------------------------------------------------------------------------------------------------------------------------
