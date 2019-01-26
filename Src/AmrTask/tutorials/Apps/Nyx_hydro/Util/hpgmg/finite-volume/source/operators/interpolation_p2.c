//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
static inline void interpolation_p2_block(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, blockCopy_type *block){
  // interpolate 3D array from read_i,j,k of read[] to write_i,j,k in write[]
  int write_dim_i   = block->dim.i<<1; // calculate the dimensions of the resultant fine block
  int write_dim_j   = block->dim.j<<1;
  int write_dim_k   = block->dim.k<<1;

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
     read_jStride = level_c->my_boxes[block->read.box ].jStride;
     read_kStride = level_c->my_boxes[block->read.box ].kStride;
     read = level_c->my_boxes[ block->read.box].vectors[id_c] + level_c->box_ghosts*(1+ read_jStride+ read_kStride);
  }
  if(block->write.box>=0){
    write_jStride = level_f->my_boxes[block->write.box].jStride;
    write_kStride = level_f->my_boxes[block->write.box].kStride;
    write = level_f->my_boxes[block->write.box].vectors[id_f] + level_f->box_ghosts*(1+write_jStride+write_kStride);
  }
 

  #ifdef USE_NAIVE_INTERP
  int i,j,k;
  double OneOver32Cubed = 1.0/32768.0;
  for(k=0;k<write_dim_k;k++){int delta_k=-read_kStride;if(k&0x1)delta_k=read_kStride;
  for(j=0;j<write_dim_j;j++){int delta_j=-read_jStride;if(j&0x1)delta_j=read_jStride;
  for(i=0;i<write_dim_i;i++){int delta_i=           -1;if(i&0x1)delta_i=           1; // i.e. even points look backwards while odd points look forward
    int write_ijk = ((i   )+write_i) + (((j   )+write_j)*write_jStride) + (((k   )+write_k)*write_kStride);
    int  read_ijk = ((i>>1)+ read_i) + (((j>>1)+ read_j)* read_jStride) + (((k>>1)+ read_k)* read_kStride);
    //
    // | -3/32 | 30/32 |  5/32 |
    // |---+---|---+---|---+---|
    // |   |   |   | x |   |   |
    //
    write[write_ijk] = prescale_f*write[write_ijk] +
                       OneOver32Cubed*(
                         -27.0*read[read_ijk-delta_i-delta_j-delta_k] +
                         270.0*read[read_ijk        -delta_j-delta_k] +
                          45.0*read[read_ijk+delta_i-delta_j-delta_k] +
                         270.0*read[read_ijk-delta_i        -delta_k] +
                       -2700.0*read[read_ijk                -delta_k] +
                        -450.0*read[read_ijk+delta_i        -delta_k] +
                          45.0*read[read_ijk-delta_i+delta_j-delta_k] +
                        -450.0*read[read_ijk        +delta_j-delta_k] +
                         -75.0*read[read_ijk+delta_i+delta_j-delta_k] +

                         270.0*read[read_ijk-delta_i-delta_j        ] +
                       -2700.0*read[read_ijk        -delta_j        ] +
                        -450.0*read[read_ijk+delta_i-delta_j        ] +
                       -2700.0*read[read_ijk-delta_i                ] +
                       27000.0*read[read_ijk                        ] +
                        4500.0*read[read_ijk+delta_i                ] +
                        -450.0*read[read_ijk-delta_i+delta_j        ] +
                        4500.0*read[read_ijk        +delta_j        ] +
                         750.0*read[read_ijk+delta_i+delta_j        ] +
                       
                          45.0*read[read_ijk-delta_i-delta_j+delta_k] +
                        -450.0*read[read_ijk        -delta_j+delta_k] +
                         -75.0*read[read_ijk+delta_i-delta_j+delta_k] +
                        -450.0*read[read_ijk-delta_i        +delta_k] +
                        4500.0*read[read_ijk                +delta_k] +
                         750.0*read[read_ijk+delta_i        +delta_k] +
                         -75.0*read[read_ijk-delta_i+delta_j+delta_k] +
                         750.0*read[read_ijk        +delta_j+delta_k] +
                         125.0*read[read_ijk+delta_i+delta_j+delta_k] 
                       );

  }}}
  #else
  int i,j,k;
  int ii,jj,kk;
  double w0 =  5.0/32.0;
  double w1 = 30.0/32.0;
  double w2 = -3.0/32.0;
  for(k=0,kk=0;k<write_dim_k;k+=2,kk++){
  for(j=0,jj=0;j<write_dim_j;j+=2,jj++){
  // compiler cannot infer/speculate write[ijk+write_jStride] is disjoint from write[ijk], so create a unique restrict pointers for each nonliteral offset...
  double * __restrict__ write00 = write + write_i + (write_j+j+0)*write_jStride + (write_k+k+0)*write_kStride;
  double * __restrict__ write10 = write + write_i + (write_j+j+1)*write_jStride + (write_k+k+0)*write_kStride;
  double * __restrict__ write01 = write + write_i + (write_j+j+0)*write_jStride + (write_k+k+1)*write_kStride;
  double * __restrict__ write11 = write + write_i + (write_j+j+1)*write_jStride + (write_k+k+1)*write_kStride;
  for(i=0,ii=0;i<write_dim_i;i+=2,ii++){
    int write_ijk = ( i+write_i) + ( j+write_j)*write_jStride + ( k+write_k)*write_kStride;
    int  read_ijk = (ii+ read_i) + (jj+ read_j)* read_jStride + (kk+ read_k)* read_kStride;
    //
    // |  5/32 | 30/32 | -3/32 | coarse grid
    // |---+---|---+---|---+---|
    // |   |   | ? |   |   |   | fine grid
    //

    // grab all coarse grid points...
    const double c000=read[read_ijk-1-read_jStride-read_kStride], c100=read[read_ijk  -read_jStride-read_kStride], c200=read[read_ijk+1-read_jStride-read_kStride];
    const double c010=read[read_ijk-1             -read_kStride], c110=read[read_ijk               -read_kStride], c210=read[read_ijk+1             -read_kStride];
    const double c020=read[read_ijk-1+read_jStride-read_kStride], c120=read[read_ijk  +read_jStride-read_kStride], c220=read[read_ijk+1+read_jStride-read_kStride];
    const double c001=read[read_ijk-1-read_jStride             ], c101=read[read_ijk  -read_jStride             ], c201=read[read_ijk+1-read_jStride             ];
    const double c011=read[read_ijk-1                          ], c111=read[read_ijk                            ], c211=read[read_ijk+1                          ];
    const double c021=read[read_ijk-1+read_jStride             ], c121=read[read_ijk  +read_jStride             ], c221=read[read_ijk+1+read_jStride             ];
    const double c002=read[read_ijk-1-read_jStride+read_kStride], c102=read[read_ijk  -read_jStride+read_kStride], c202=read[read_ijk+1-read_jStride+read_kStride];
    const double c012=read[read_ijk-1             +read_kStride], c112=read[read_ijk               +read_kStride], c212=read[read_ijk+1             +read_kStride];
    const double c022=read[read_ijk-1+read_jStride+read_kStride], c122=read[read_ijk  +read_jStride+read_kStride], c222=read[read_ijk+1+read_jStride+read_kStride];

    // interpolate in i to create fine i / coarse jk points...
    //
    // +-------+-------+-------+      :.......+---+---+.......:
    // |       |       |       |      :       |   |   |       :
    // |   c   |   c   |   c   |      :       | f | f |       :
    // |       |       |       |      :       |   |   |       :
    // +-------+-------+-------+      :.......+---+---+.......:
    // |       |       |       |      :       |   |   |       :
    // |   c   |   c   |   c   |  ->  :       | f | f |       :
    // |       |       |       |      :       |   |   |       :
    // +-------+-------+-------+      :.......+---+---+.......:
    // |       |       |       |      :       |   |   |       :
    // |   c   |   c   |   c   |      :       | f | f |       :
    // |       |       |       |      :       |   |   |       :
    // +-------+-------+-------+      :.......+---+---+.......:
    //
    const double f0c00 = ( w1*c100 + w0*c000 + w2*c200 );
    const double f1c00 = ( w1*c100 + w2*c000 + w0*c200 );
    const double f0c10 = ( w1*c110 + w0*c010 + w2*c210 );
    const double f1c10 = ( w1*c110 + w2*c010 + w0*c210 );
    const double f0c20 = ( w1*c120 + w0*c020 + w2*c220 );
    const double f1c20 = ( w1*c120 + w2*c020 + w0*c220 );

    const double f0c01 = ( w1*c101 + w0*c001 + w2*c201 );
    const double f1c01 = ( w1*c101 + w2*c001 + w0*c201 );
    const double f0c11 = ( w1*c111 + w0*c011 + w2*c211 );
    const double f1c11 = ( w1*c111 + w2*c011 + w0*c211 );
    const double f0c21 = ( w1*c121 + w0*c021 + w2*c221 );
    const double f1c21 = ( w1*c121 + w2*c021 + w0*c221 );

    const double f0c02 = ( w1*c102 + w0*c002 + w2*c202 );
    const double f1c02 = ( w1*c102 + w2*c002 + w0*c202 );
    const double f0c12 = ( w1*c112 + w0*c012 + w2*c212 );
    const double f1c12 = ( w1*c112 + w2*c012 + w0*c212 );
    const double f0c22 = ( w1*c122 + w0*c022 + w2*c222 );
    const double f1c22 = ( w1*c122 + w2*c022 + w0*c222 );

    // interpolate in j to create fine ij / coarse k points...
    //
    // :.......+---+---+.......:      :.......:.......:.......:
    // :       |   |   |       :      :       :       :       :
    // :       |   |   |       :      :       :       :       :
    // :       |   |   |       :      :       :       :       :
    // :.......+---+---+.......:      :.......+---+---+.......:
    // :       |   |   |       :      :       |   |   |       :
    // :       |   |   |       :  ->  :       +---+---+       :
    // :       |   |   |       :      :       |   |   |       :
    // :.......+---+---+.......:      :.......+---+---+.......:
    // :       |   |   |       :      :       :       :       :
    // :       |   |   |       :      :       :       :       :
    // :       |   |   |       :      :       :       :       :
    // :.......+---+---+.......:      :.......:.......:.......:
    //
    const double f00c0 = ( w1*f0c10 + w0*f0c00 + w2*f0c20 );
    const double f10c0 = ( w1*f1c10 + w0*f1c00 + w2*f1c20 );
    const double f01c0 = ( w1*f0c10 + w2*f0c00 + w0*f0c20 );
    const double f11c0 = ( w1*f1c10 + w2*f1c00 + w0*f1c20 );

    const double f00c1 = ( w1*f0c11 + w0*f0c01 + w2*f0c21 );
    const double f10c1 = ( w1*f1c11 + w0*f1c01 + w2*f1c21 );
    const double f01c1 = ( w1*f0c11 + w2*f0c01 + w0*f0c21 );
    const double f11c1 = ( w1*f1c11 + w2*f1c01 + w0*f1c21 );

    const double f00c2 = ( w1*f0c12 + w0*f0c02 + w2*f0c22 );
    const double f10c2 = ( w1*f1c12 + w0*f1c02 + w2*f1c22 );
    const double f01c2 = ( w1*f0c12 + w2*f0c02 + w0*f0c22 );
    const double f11c2 = ( w1*f1c12 + w2*f1c02 + w0*f1c22 );

    // interpolate in k to create fine ijk points...
    const double f000 = ( w1*f00c1 + w0*f00c0 + w2*f00c2 );
    const double f100 = ( w1*f10c1 + w0*f10c0 + w2*f10c2 );
    const double f010 = ( w1*f01c1 + w0*f01c0 + w2*f01c2 );
    const double f110 = ( w1*f11c1 + w0*f11c0 + w2*f11c2 );
    const double f001 = ( w1*f00c1 + w2*f00c0 + w0*f00c2 );
    const double f101 = ( w1*f10c1 + w2*f10c0 + w0*f10c2 );
    const double f011 = ( w1*f01c1 + w2*f01c0 + w0*f01c2 );
    const double f111 = ( w1*f11c1 + w2*f11c0 + w0*f11c2 );

    // commit to memory...
    #if 0 // compiler cannot infer/speculate write[ijk+write_jStride] is disjoint from write[ijk], and thus cannot vectorize...
    write[write_ijk                              ] = prescale_f*write[write_ijk                              ] + f000;
    write[write_ijk+1                            ] = prescale_f*write[write_ijk+1                            ] + f100;
    write[write_ijk  +write_jStride              ] = prescale_f*write[write_ijk  +write_jStride              ] + f010;
    write[write_ijk+1+write_jStride              ] = prescale_f*write[write_ijk+1+write_jStride              ] + f110;
    write[write_ijk                +write_kStride] = prescale_f*write[write_ijk                +write_kStride] + f001;
    write[write_ijk+1              +write_kStride] = prescale_f*write[write_ijk+1              +write_kStride] + f101;
    write[write_ijk  +write_jStride+write_kStride] = prescale_f*write[write_ijk  +write_jStride+write_kStride] + f011;
    write[write_ijk+1+write_jStride+write_kStride] = prescale_f*write[write_ijk+1+write_jStride+write_kStride] + f111;
    #else // use a unique restrict pointer for each pencil...
    write00[i  ] = prescale_f*write00[i  ] + f000;
    write00[i+1] = prescale_f*write00[i+1] + f100;
    write10[i  ] = prescale_f*write10[i  ] + f010;
    write10[i+1] = prescale_f*write10[i+1] + f110;
    write01[i  ] = prescale_f*write01[i  ] + f001;
    write01[i+1] = prescale_f*write01[i+1] + f101;
    write11[i  ] = prescale_f*write11[i  ] + f011;
    write11[i+1] = prescale_f*write11[i+1] + f111;
    #endif

  }}}
  #endif

}


//------------------------------------------------------------------------------------------------------------------------------
// perform a (inter-level) piecewise quadratic interpolation
void interpolation_p2(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){
    exchange_boundary(level_c,id_c,STENCIL_SHAPE_BOX);
         apply_BCs_p2(level_c,id_c,STENCIL_SHAPE_BOX);

  double _timeCommunicationStart = getTime();
  double _timeStart,_timeEnd;
  int buffer=0;
  int n;
  int my_tag = (level_f->tag<<4) | 0x7;


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
      interpolation_p2_block(level_f,id_f,0.0,level_c,id_c,&level_c->interpolation.blocks[0][buffer]);
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
      interpolation_p2_block(level_f,id_f,prescale_f,level_c,id_c,&level_c->interpolation.blocks[1][buffer]);
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
