//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void matmul(level_type * level, double *C, int * id_A, int * id_B, int rows, int cols, int A_equals_B_transpose){
  // *id_A = m vector_id's (conceptually pointers to the rows    of a m x level->num_my_boxes*volume matrix)
  // *id_B = n vector_id's (conceptually pointers to the columns of a level->num_my_boxes*volume matrix x n)
  // *C is a mxn matrix where C[rows][cols] = dot(id_A[rows],id_B[cols])

  // FIX, id_A and id_B are likely the same and thus C[][] will be symmetric (modulo missing row?)
  // if(A_equals_B_transpose && (cols>=rows)) then use id_B and only run for nn>=mm // common case for s-step Krylov methods
  // C_is_symmetric && cols< rows (use id_A)
  int mm,nn;


  double _timeStart = getTime();
  // FIX... rather than performing an all_reduce on the essentially symmetric [G,g], do the all_reduce on the upper triangle and then duplicate (saves BW)
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static,1) collapse(2)
  #endif
  for(mm=0;mm<rows;mm++){
  for(nn=0;nn<cols;nn++){
  if(nn>=mm){ // upper triangular
    int box;
    double a_dot_b_level =  0.0;
    for(box=0;box<level->num_my_boxes;box++){
      int i,j,k;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int  ghosts = level->my_boxes[box].ghosts;
      const int     dim = level->my_boxes[box].dim;
      double * __restrict__ grid_a = level->my_boxes[box].vectors[id_A[mm]] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
      double * __restrict__ grid_b = level->my_boxes[box].vectors[id_B[nn]] + ghosts*(1+jStride+kStride); 
      double a_dot_b_box = 0.0;
      for(k=0;k<dim;k++){
      for(j=0;j<dim;j++){
      for(i=0;i<dim;i++){
        int ijk = i + j*jStride + k*kStride;
        a_dot_b_box += grid_a[ijk]*grid_b[ijk];
      }}}
      a_dot_b_level+=a_dot_b_box;
    }
                             C[mm*cols + nn] = a_dot_b_level; // C[mm][nn]
    if((mm<cols)&&(nn<rows)){C[nn*cols + mm] = a_dot_b_level;}// C[nn][mm] 
  }
  }}
  level->timers.blas3 += (double)(getTime()-_timeStart);

  #ifdef USE_MPI
  double *send_buffer = (double*)malloc(rows*cols*sizeof(double));
  for(mm=0;mm<rows;mm++){
  for(nn=0;nn<cols;nn++){
    send_buffer[mm*cols + nn] = C[mm*cols + nn];
  }}
  double _timeStartAllReduce = getTime();
  MPI_Allreduce(send_buffer,C,rows*cols,MPI_DOUBLE,MPI_SUM,level->MPI_COMM_ALLREDUCE);
  double _timeEndAllReduce = getTime();
  level->timers.collectives   += (double)(_timeEndAllReduce-_timeStartAllReduce);
  free(send_buffer);
  #endif

}

