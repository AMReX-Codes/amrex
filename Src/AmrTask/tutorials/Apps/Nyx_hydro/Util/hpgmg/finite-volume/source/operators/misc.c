//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void zero_vector(level_type * level, int id_a){
  // zero's the entire grid INCLUDING ghost zones...
  double _timeStart = getTime();
  int block;

  PRAGMA_THREAD_ACROSS_BLOCKS(level,block,level->num_my_blocks)
  for(block=0;block<level->num_my_blocks;block++){
    const int box = level->my_blocks[block].read.box;
          int ilo = level->my_blocks[block].read.i;
          int jlo = level->my_blocks[block].read.j;
          int klo = level->my_blocks[block].read.k;
          int ihi = level->my_blocks[block].dim.i + ilo;
          int jhi = level->my_blocks[block].dim.j + jlo;
          int khi = level->my_blocks[block].dim.k + klo;
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;

    // expand the size of the block to include the ghost zones...
    if(ilo<=  0)ilo-=ghosts; 
    if(jlo<=  0)jlo-=ghosts; 
    if(klo<=  0)klo-=ghosts; 
    if(ihi>=dim)ihi+=ghosts; 
    if(jhi>=dim)jhi+=ghosts; 
    if(khi>=dim)khi+=ghosts; 

    double * __restrict__ grid = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
      int ijk = i + j*jStride + k*kStride;
      grid[ijk] = 0.0;
    }}}
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
void init_vector(level_type * level, int id_a, double scalar){
  // initializes the grid to a scalar while zero'ing the ghost zones...
  double _timeStart = getTime();
  int block;

  PRAGMA_THREAD_ACROSS_BLOCKS(level,block,level->num_my_blocks)
  for(block=0;block<level->num_my_blocks;block++){
    const int box = level->my_blocks[block].read.box;
          int ilo = level->my_blocks[block].read.i;
          int jlo = level->my_blocks[block].read.j;
          int klo = level->my_blocks[block].read.k;
          int ihi = level->my_blocks[block].dim.i + ilo;
          int jhi = level->my_blocks[block].dim.j + jlo;
          int khi = level->my_blocks[block].dim.k + klo;
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;

    // expand the size of the block to include the ghost zones...
    if(ilo<=  0)ilo-=ghosts; 
    if(jlo<=  0)jlo-=ghosts; 
    if(klo<=  0)klo-=ghosts; 
    if(ihi>=dim)ihi+=ghosts; 
    if(jhi>=dim)jhi+=ghosts; 
    if(khi>=dim)khi+=ghosts; 

    double * __restrict__ grid = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
        int ijk = i + j*jStride + k*kStride;
        int ghostZone = (i<0) || (j<0) || (k<0) || (i>=dim) || (j>=dim) || (k>=dim);
        grid[ijk] = ghostZone ? 0.0 : scalar;
    }}}
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
// add vectors id_a (scaled by scale_a) and id_b (scaled by scale_b) and store the result in vector id_c
// i.e. c[] = scale_a*a[] + scale_b*b[]
// note, only non ghost zone values are included in this calculation
void add_vectors(level_type * level, int id_c, double scale_a, int id_a, double scale_b, int id_b){
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
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_b = level->my_boxes[box].vectors[id_b] + ghosts*(1+jStride+kStride);

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
        int ijk = i + j*jStride + k*kStride;
        grid_c[ijk] = scale_a*grid_a[ijk] + scale_b*grid_b[ijk];
    }}}
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
// multiply each element of vector id_a by vector id_b and scale, and place the result in vector id_c
// i.e. c[]=scale*a[]*b[]
// note, only non ghost zone values are included in this calculation
void mul_vectors(level_type * level, int id_c, double scale, int id_a, int id_b){
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
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_b = level->my_boxes[box].vectors[id_b] + ghosts*(1+jStride+kStride);

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
        int ijk = i + j*jStride + k*kStride;
        grid_c[ijk] = scale*grid_a[ijk]*grid_b[ijk];
    }}}
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
// invert each element of vector id_a, scale by scale_a, and place the result in vector id_c
// i.e. c[]=scale_a/a[]
// note, only non ghost zone values are included in this calculation
void invert_vector(level_type * level, int id_c, double scale_a, int id_a){
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
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
        int ijk = i + j*jStride + k*kStride;
        grid_c[ijk] = scale_a/grid_a[ijk];
    }}}
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
// scale vector id_a by scale_a and place the result in vector id_c
// i.e. c[]=scale_a*a[]
// note, only non ghost zone values are included in this calculation
void scale_vector(level_type * level, int id_c, double scale_a, int id_a){
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
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
        int ijk = i + j*jStride + k*kStride;
        grid_c[ijk] = scale_a*grid_a[ijk];
    }}}
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
// return the dot product of vectors id_a and id_b
// note, only non ghost zone values are included in this calculation
double dot(level_type * level, int id_a, int id_b){
  double _timeStart = getTime();


  int block;
  double a_dot_b_level =  0.0;

  PRAGMA_THREAD_ACROSS_BLOCKS_SUM(level,block,level->num_my_blocks,a_dot_b_level)
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
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double * __restrict__ grid_b = level->my_boxes[box].vectors[id_b] + ghosts*(1+jStride+kStride);
    double a_dot_b_block = 0.0;

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
      int ijk = i + j*jStride + k*kStride;
      a_dot_b_block += grid_a[ijk]*grid_b[ijk];
    }}}
    a_dot_b_level+=a_dot_b_block;
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);

  #ifdef USE_MPI
  double _timeStartAllReduce = getTime();
  double send = a_dot_b_level;
  MPI_Allreduce(&send,&a_dot_b_level,1,MPI_DOUBLE,MPI_SUM,level->MPI_COMM_ALLREDUCE);
  double _timeEndAllReduce = getTime();
  level->timers.collectives   += (double)(_timeEndAllReduce-_timeStartAllReduce);
  #endif

  return(a_dot_b_level);
}

//------------------------------------------------------------------------------------------------------------------------------
// return the max (infinity) norm of the vector id_a.
// note, only non ghost zone values are included in this calculation
double norm(level_type * level, int id_a){ // implements the max norm
  double _timeStart = getTime();

  int block;
  double max_norm =  0.0;

  PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,block,level->num_my_blocks,max_norm)
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
    double * __restrict__ grid   = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double block_norm = 0.0;

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){ 
      int ijk = i + j*jStride + k*kStride;
      double fabs_grid_ijk = fabs(grid[ijk]);
      if(fabs_grid_ijk>block_norm){block_norm=fabs_grid_ijk;} // max norm
    }}}

    if(block_norm>max_norm){max_norm = block_norm;}
  } // block list
  level->timers.blas1 += (double)(getTime()-_timeStart);

  #ifdef USE_MPI
  double _timeStartAllReduce = getTime();
  double send = max_norm;
  MPI_Allreduce(&send,&max_norm,1,MPI_DOUBLE,MPI_MAX,level->MPI_COMM_ALLREDUCE);
  double _timeEndAllReduce = getTime();
  level->timers.collectives   += (double)(_timeEndAllReduce-_timeStartAllReduce);
  #endif
  return(max_norm);
}


//------------------------------------------------------------------------------------------------------------------------------
// return the mean (arithmetic average value) of vector id_a
// essentially, this is a l1 norm by a scaling by the inverse of the total (global) number of cells
// note, only non ghost zone values are included in this calculation
double mean(level_type * level, int id_a){
  double _timeStart = getTime();


  int block;
  double sum_level =  0.0;

  PRAGMA_THREAD_ACROSS_BLOCKS_SUM(level,block,level->num_my_blocks,sum_level)
  for(block=0;block<level->num_my_blocks;block++){
    const int box = level->my_blocks[block].read.box;
    const int ilo = level->my_blocks[block].read.i;
    const int jlo = level->my_blocks[block].read.j;
    const int klo = level->my_blocks[block].read.k;
    const int ihi = level->my_blocks[block].dim.i + ilo;
    const int jhi = level->my_blocks[block].dim.j + jlo;
    const int khi = level->my_blocks[block].dim.k + klo;
    int i,j,k;
    int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double sum_block = 0.0;

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
      int ijk = i + j*jStride + k*kStride;
      sum_block += grid_a[ijk];
    }}}
    sum_level+=sum_block;
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);
  double ncells_level = (double)level->dim.i*(double)level->dim.j*(double)level->dim.k;

  #ifdef USE_MPI
  double _timeStartAllReduce = getTime();
  double send = sum_level;
  MPI_Allreduce(&send,&sum_level,1,MPI_DOUBLE,MPI_SUM,level->MPI_COMM_ALLREDUCE);
  double _timeEndAllReduce = getTime();
  level->timers.collectives   += (double)(_timeEndAllReduce-_timeStartAllReduce);
  #endif

  double mean_level = sum_level / ncells_level;
  return(mean_level);
}


//------------------------------------------------------------------------------------------------------------------------------
// add the scalar value shift_a to each element of vector id_a and store the result in vector id_c
// note, only non ghost zone values are included in this calculation
void shift_vector(level_type * level, int id_c, int id_a, double shift_a){
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
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point


    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
      int ijk = i + j*jStride + k*kStride;
      grid_c[ijk] = grid_a[ijk] + shift_a;
    }}}
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);
}

//------------------------------------------------------------------------------------------------------------------------------
// calculate the error between two vectors (id_a and id_b) using either the max (infinity) norm or the L2 norm
// note, only non ghost zone values are included in this calculation
double error(level_type * level, int id_a, int id_b){
  double h3 = level->h * level->h * level->h;
               add_vectors(level,VECTOR_TEMP,1.0,id_a,-1.0,id_b);            // VECTOR_TEMP = id_a - id_b
  double   max =      norm(level,VECTOR_TEMP);                return(max);   // max norm of error function
  double    L2 = sqrt( dot(level,VECTOR_TEMP,VECTOR_TEMP)*h3);return( L2);   // normalized L2 error ?
}


//------------------------------------------------------------------------------------------------------------------------------
// Color the vector id_a with 1's and 0's
// The pattern is dictated by the number of colors in each dimension and the 'active' color (i,j,kcolor)
// note, only non ghost zone values are included in this calculation
//   e.g. colors_in_each_dim=3, icolor=1, jcolor=2...
//   -+---+---+---+-
//    | 0 | 1 | 0 |
//   -+---+---+---+-
//    | 0 | 0 | 0 |
//   -+---+---+---+-
//    | 0 | 0 | 0 |
//   -+---+---+---+-
//
void color_vector(level_type * level, int id_a, int colors_in_each_dim, int icolor, int jcolor, int kcolor){
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
    const int boxlowi = level->my_boxes[box].low.i;
    const int boxlowj = level->my_boxes[box].low.j;
    const int boxlowk = level->my_boxes[box].low.k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    double * __restrict__ grid = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    int i,j,k;

    for(k=klo;k<khi;k++){double sk=0.0;if( ((k+boxlowk+kcolor)%colors_in_each_dim) == 0 )sk=1.0; // if colors_in_each_dim==1 (don't color), all cells are set to 1.0
    for(j=jlo;j<jhi;j++){double sj=0.0;if( ((j+boxlowj+jcolor)%colors_in_each_dim) == 0 )sj=1.0;
    for(i=ilo;i<ihi;i++){double si=0.0;if( ((i+boxlowi+icolor)%colors_in_each_dim) == 0 )si=1.0;
      int ijk = i + j*jStride + k*kStride;
      grid[ijk] = si*sj*sk;
    }}}
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
// Initialize each element of vector id_a with a "random" value.  
// For simplicity, random is defined as -1.0 or +1.0 and is based on whether the coordinates of the element are even or odd
// note, only non ghost zone values are included in this calculation
void random_vector(level_type * level, int id_a){
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
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    double * __restrict__ grid = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    int i,j,k;

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){
      int ijk = i + j*jStride + k*kStride;
      grid[ijk] = -1.000 + 2.0*(i^j^k^0x1);
    }}}
  }
  level->timers.blas1 += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
