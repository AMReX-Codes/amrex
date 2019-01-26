//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
static inline void CopyBlock(level_type *level, int id, blockCopy_type *block){
  // copy 3D array from read_i,j,k of read[] to write_i,j,k in write[]
  int   dim_i       = block->dim.i;
  int   dim_j       = block->dim.j;
  int   dim_k       = block->dim.k;

  int  read_i       = block->read.i;
  int  read_j       = block->read.j;
  int  read_k       = block->read.k;
  int  read_jStride = block->read.jStride;
  int  read_kStride = block->read.kStride;

  int write_i       = block->write.i;
  int write_j       = block->write.j;
  int write_k       = block->write.k;
  int write_jStride = block->write.jStride;
  int write_kStride = block->write.kStride;

  const double * __restrict__  read = block->read.ptr;
        double * __restrict__ write = block->write.ptr;

  if(block->read.box >=0){
     read_jStride = level->my_boxes[block->read.box ].jStride;
     read_kStride = level->my_boxes[block->read.box ].kStride;
     read = level->my_boxes[ block->read.box].vectors[id] + level->box_ghosts*(1+ read_jStride+ read_kStride);
  }
  if(block->write.box>=0){
    write_jStride = level->my_boxes[block->write.box].jStride;
    write_kStride = level->my_boxes[block->write.box].kStride;
    write = level->my_boxes[block->write.box].vectors[id] + level->box_ghosts*(1+write_jStride+write_kStride);
  }


  int i,j,k;
  if(dim_i==1){ // be smart and don't have an inner loop from 0 to 0
    for(k=0;k<dim_k;k++){
    for(j=0;j<dim_j;j++){
      int  read_ijk = ( read_i) + (j+ read_j)* read_jStride + (k+ read_k)* read_kStride;
      int write_ijk = (write_i) + (j+write_j)*write_jStride + (k+write_k)*write_kStride;
      write[write_ijk] = read[read_ijk];
    }}
  }else if(dim_i==2){ // be smart and don't have an inner loop from 0 to 1
    for(k=0;k<dim_k;k++){
    for(j=0;j<dim_j;j++){
      int  read_ijk = ( read_i) + (j+ read_j)* read_jStride + (k+ read_k)* read_kStride;
      int write_ijk = (write_i) + (j+write_j)*write_jStride + (k+write_k)*write_kStride;
      write[write_ijk+0] = read[read_ijk+0];
      write[write_ijk+1] = read[read_ijk+1];
    }}
  }else if(dim_i==4){ // be smart and don't have an inner loop from 0 to 3
    for(k=0;k<dim_k;k++){
    for(j=0;j<dim_j;j++){
      int  read_ijk = ( read_i) + (j+ read_j)* read_jStride + (k+ read_k)* read_kStride;
      int write_ijk = (write_i) + (j+write_j)*write_jStride + (k+write_k)*write_kStride;
      write[write_ijk+0] = read[read_ijk+0];
      write[write_ijk+1] = read[read_ijk+1];
      write[write_ijk+2] = read[read_ijk+2];
      write[write_ijk+3] = read[read_ijk+3];
    }}
  }else if(dim_j==1){ // don't have a 0..0 loop
    for(k=0;k<dim_k;k++){
    for(i=0;i<dim_i;i++){
      int  read_ijk = (i+ read_i) + ( read_j)* read_jStride + (k+ read_k)* read_kStride;
      int write_ijk = (i+write_i) + (write_j)*write_jStride + (k+write_k)*write_kStride;
      write[write_ijk] = read[read_ijk];
    }}
  }else if(dim_k==1){ // don't have a 0..0 loop
    for(j=0;j<dim_j;j++){
    for(i=0;i<dim_i;i++){
      int  read_ijk = (i+ read_i) + (j+ read_j)* read_jStride + ( read_k)* read_kStride;
      int write_ijk = (i+write_i) + (j+write_j)*write_jStride + (write_k)*write_kStride;
      write[write_ijk] = read[read_ijk];
    }}
  }else{ // general case...
    for(k=0;k<dim_k;k++){
    for(j=0;j<dim_j;j++){
    for(i=0;i<dim_i;i++){
      int  read_ijk = (i+ read_i) + (j+ read_j)* read_jStride + (k+ read_k)* read_kStride;
      int write_ijk = (i+write_i) + (j+write_j)*write_jStride + (k+write_k)*write_kStride;
      write[write_ijk] = read[read_ijk];
    }}}
  }

}


//------------------------------------------------------------------------------------------------------------------------------
static inline void IncrementBlock(level_type *level, int id, double prescale, blockCopy_type *block){
  // copy 3D array from read_i,j,k of read[] to write_i,j,k in write[]
  int   dim_i       = block->dim.i;
  int   dim_j       = block->dim.j;
  int   dim_k       = block->dim.k;

  int  read_i       = block->read.i;
  int  read_j       = block->read.j;
  int  read_k       = block->read.k;
  int  read_jStride = block->read.jStride;
  int  read_kStride = block->read.kStride;

  int write_i       = block->write.i;
  int write_j       = block->write.j;
  int write_k       = block->write.k;
  int write_jStride = block->write.jStride;
  int write_kStride = block->write.kStride;

  const double * __restrict__  read = block->read.ptr;
        double * __restrict__ write = block->write.ptr;

  if(block->read.box >=0){
     read_jStride = level->my_boxes[block->read.box ].jStride;
     read_kStride = level->my_boxes[block->read.box ].kStride;
     read = level->my_boxes[ block->read.box].vectors[id] + level->box_ghosts*(1+ read_jStride+ read_kStride);
  }
  if(block->write.box>=0){
    write_jStride = level->my_boxes[block->write.box].jStride;
    write_kStride = level->my_boxes[block->write.box].kStride;
    write = level->my_boxes[block->write.box].vectors[id] + level->box_ghosts*(1+write_jStride+write_kStride);
  }

  int i,j,k;
  for(k=0;k<dim_k;k++){
  for(j=0;j<dim_j;j++){
  for(i=0;i<dim_i;i++){
    int  read_ijk = (i+ read_i) + (j+ read_j)* read_jStride + (k+ read_k)* read_kStride;
    int write_ijk = (i+write_i) + (j+write_j)*write_jStride + (k+write_k)*write_kStride;
    write[write_ijk] = prescale*write[write_ijk] + read[read_ijk]; // CAREFUL !!!  you must guarantee you zero'd the MPI buffers(write[]) and destination boxes at some point to avoid 0.0*NaN or 0.0*inf
  }}}

}

//------------------------------------------------------------------------------------------------------------------------------
