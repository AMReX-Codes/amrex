//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void apply_op(level_type * level, int Ax_id, int x_id, double a, double b){  // y=Ax
  // exchange the boundary of x in preparation for Ax
  exchange_boundary(level,x_id,stencil_get_shape());
          apply_BCs(level,x_id,stencil_get_shape());

  // now do Ax proper...
  double _timeStart = getTime();
  const int  ghosts = level->box_ghosts;
  const int jStride = level->box_jStride;
  const int kStride = level->box_kStride;
  const int     dim = level->box_dim;
  const double h2inv = 1.0/(level->h*level->h);
  int box;

  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const double * __restrict__ x      = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
          double * __restrict__ Ax     = level->my_boxes[box].vectors[        Ax_id] + ghosts*(1+jStride+kStride); 
    const double * __restrict__ alpha  = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_i = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_j = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_k = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride);

    PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
      int ijk = i + j*jStride + k*kStride;
      Ax[ijk] = apply_op_ijk(x);
    }}}
  }
  level->timers.apply_op += (double)(getTime()-_timeStart);
}
//------------------------------------------------------------------------------------------------------------------------------
