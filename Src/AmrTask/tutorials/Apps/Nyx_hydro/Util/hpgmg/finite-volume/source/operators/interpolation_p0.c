//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
static inline void interpolation_p0_block(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, blockCopy_type *block){
  // interpolate 3D array from read_i,j,k of read[] to write_i,j,k in write[]
  int   dim_i       = block->dim.i<<1; // calculate the dimensions of the resultant fine block
  int   dim_j       = block->dim.j<<1;
  int   dim_k       = block->dim.k<<1;

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

  double * __restrict__  read = block->read.ptr;
  double * __restrict__ write = block->write.ptr;
  if(block->read.box >=0){
     read = level_c->my_boxes[ block->read.box].vectors[id_c] + level_c->my_boxes[ block->read.box].ghosts*(1+level_c->my_boxes[ block->read.box].jStride+level_c->my_boxes[ block->read.box].kStride);
     read_jStride = level_c->my_boxes[block->read.box ].jStride;
     read_kStride = level_c->my_boxes[block->read.box ].kStride;
  }
  if(block->write.box>=0){
    write = level_f->my_boxes[block->write.box].vectors[id_f] + level_f->my_boxes[block->write.box].ghosts*(1+level_f->my_boxes[block->write.box].jStride+level_f->my_boxes[block->write.box].kStride);
    write_jStride = level_f->my_boxes[block->write.box].jStride;
    write_kStride = level_f->my_boxes[block->write.box].kStride;
  }
 
 
  int i,j,k;
  for(k=0;k<dim_k;k++){
  for(j=0;j<dim_j;j++){
  for(i=0;i<dim_i;i++){
    int write_ijk = ((i   )+write_i) + (((j   )+write_j)*write_jStride) + (((k   )+write_k)*write_kStride);
    int  read_ijk = ((i>>1)+ read_i) + (((j>>1)+ read_j)* read_jStride) + (((k>>1)+ read_k)* read_kStride);
    write[write_ijk] = prescale_f*write[write_ijk] + read[read_ijk]; // CAREFUL !!!  you must guarantee you zero'd the MPI buffers(write[]) and destination boxes at some point to avoid 0.0*NaN or 0.0*inf
  }}}

}


//------------------------------------------------------------------------------------------------------------------------------
// perform a (inter-level) piecewise constant interpolation
void interpolation_p0(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){
  double _timeCommunicationStart = getTime();
  double _timeStart,_timeEnd;
  int my_tag = (level_f->tag<<4) | 0x6;
  int buffer=0;
  int n;


  #ifdef USE_MPI
  // by convention, level_f allocates a combined array of requests for both level_f recvs and level_c sends...
  int nMessages = level_c->interpolation.num_sends + level_f->interpolation.num_recvs;
  MPI_Request *recv_requests = level_f->interpolation.requests;
  MPI_Request *send_requests = level_f->interpolation.requests + level_f->interpolation.num_recvs;


  // loop through packed list of MPI receives and prepost Irecv's...
  if(level_f->interpolation.num_recvs>0){
    _timeStart = getTime();
    #ifdef USE_MPI_THREAD_MULTIPLE
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(n=0;n<level_f->interpolation.num_recvs;n++){
      MPI_Irecv(level_f->interpolation.recv_buffers[n],
                level_f->interpolation.recv_sizes[n],
                MPI_DOUBLE,
                level_f->interpolation.recv_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &recv_requests[n]
      );
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_recv += (_timeEnd-_timeStart);
  }


  // pack MPI send buffers...
  if(level_c->interpolation.num_blocks[0]>0){
    _timeStart = getTime();
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_c->interpolation.num_blocks[0])
    for(buffer=0;buffer<level_c->interpolation.num_blocks[0];buffer++){
      // !!! prescale==0 because you don't want to increment the MPI buffer
      interpolation_p0_block(level_f,id_f,0.0,level_c,id_c,&level_c->interpolation.blocks[0][buffer]);
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_pack += (_timeEnd-_timeStart);
  }


  // loop through MPI send buffers and post Isend's...
  if(level_c->interpolation.num_sends>0){
    _timeStart = getTime();
    #ifdef USE_MPI_THREAD_MULTIPLE
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(n=0;n<level_c->interpolation.num_sends;n++){
      MPI_Isend(level_c->interpolation.send_buffers[n],
                level_c->interpolation.send_sizes[n],
                MPI_DOUBLE,
                level_c->interpolation.send_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &send_requests[n]
      );
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_send += (_timeEnd-_timeStart);
  }
  #endif


  // perform local interpolation... try and hide within Isend latency... 
  if(level_c->interpolation.num_blocks[1]>0){
    _timeStart = getTime();
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_c->interpolation.num_blocks[1])
    for(buffer=0;buffer<level_c->interpolation.num_blocks[1];buffer++){
      interpolation_p0_block(level_f,id_f,prescale_f,level_c,id_c,&level_c->interpolation.blocks[1][buffer]);
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_local += (_timeEnd-_timeStart);
  }


  // wait for MPI to finish...
  #ifdef USE_MPI 
  if(nMessages>0){
    _timeStart = getTime();
    MPI_Waitall(nMessages,level_f->interpolation.requests,level_f->interpolation.status);
    _timeEnd = getTime();
    level_f->timers.interpolation_wait += (_timeEnd-_timeStart);
  }


  // unpack MPI receive buffers 
  if(level_f->interpolation.num_blocks[2]>0){
    _timeStart = getTime();
    PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_f->interpolation.num_blocks[2])
    for(buffer=0;buffer<level_f->interpolation.num_blocks[2];buffer++){
      IncrementBlock(level_f,id_f,prescale_f,&level_f->interpolation.blocks[2][buffer]);
    }
    _timeEnd = getTime();
    level_f->timers.interpolation_unpack += (_timeEnd-_timeStart);
  }
  #endif 
 
 
  level_f->timers.interpolation_total += (double)(getTime()-_timeCommunicationStart);
}
