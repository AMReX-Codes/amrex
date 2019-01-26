//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// This routines calculates the residual (res=rhs-Ax) using the linear operator specified in the apply_op_ijk macro
// This requires exchanging a ghost zone and/or enforcing a boundary condition.
// NOTE, x_id must be distinct from rhs_id and res_id
void residual(level_type * level, int res_id, int x_id, int rhs_id, double a, double b){
  // exchange the boundary for x in prep for Ax...
  exchange_boundary(level,x_id,stencil_get_shape());
          apply_BCs(level,x_id,stencil_get_shape());

  // now do residual/restriction proper...
  double _timeStart = getTime();
  int block;

  PRAGMA_THREAD_ACROSS_BLOCKS(level,block,level->num_my_blocks)
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
    const double h2inv = 1.0/(level->h*level->h);
    const double * __restrict__ x      = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    const double * __restrict__ rhs    = level->my_boxes[box].vectors[       rhs_id] + ghosts*(1+jStride+kStride);
    const double * __restrict__ alpha  = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_i = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_j = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_k = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride);
          double * __restrict__ res    = level->my_boxes[box].vectors[       res_id] + ghosts*(1+jStride+kStride);

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
      int ijk = i + j*jStride + k*kStride;
      double Ax = apply_op_ijk(x);
      res[ijk] = rhs[ijk]-Ax;
    }}}
  }
  level->timers.residual += (double)(getTime()-_timeStart);
}

