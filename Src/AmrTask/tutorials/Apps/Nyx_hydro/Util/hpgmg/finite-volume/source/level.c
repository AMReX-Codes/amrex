//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
//------------------------------------------------------------------------------------------------------------------------------
#include "timers.h"
#include "defines.h"
#include "level.h"
#include "operators.h"
//------------------------------------------------------------------------------------------------------------------------------
void print_communicator(int printSendRecv, int rank, int level, communicator_type *comm){
  int i;
  printf("rank=%2d level=%d ",rank,level);
  if(printSendRecv & 0x1){
  printf("num_sends=%2d ",comm->num_sends);
  printf("send_ranks=[ ");for(i=0;i<comm->num_sends;i++)printf("%2d ",comm->send_ranks[i]);printf("] ");
  printf("send_sizes=[ ");for(i=0;i<comm->num_sends;i++)printf("%2d ",comm->send_sizes[i]);printf("] ");
  printf("send_buffers=[ ");for(i=0;i<comm->num_sends;i++)printf("%08lx ",(uint64_t)comm->send_buffers[i]);printf("] ");
  for(i=0;i<comm->num_blocks[0];i++)printf("[ %dx%dx%d from %d %d %d %d %d to %d %d %d %d %d ] ",comm->blocks[0][i].dim.i,comm->blocks[0][i].dim.j,comm->blocks[0][i].dim.k,comm->blocks[0][i].read.i,comm->blocks[0][i].read.j,comm->blocks[0][i].read.k,comm->blocks[0][i].read.jStride,comm->blocks[0][i].read.kStride,comm->blocks[0][i].write.i,comm->blocks[0][i].write.j,comm->blocks[0][i].write.k,comm->blocks[0][i].write.jStride,comm->blocks[0][i].write.kStride);
  printf("\n");
  }
  if(printSendRecv & 0x2){
  for(i=0;i<comm->num_blocks[1];i++)printf("[ %dx%dx%d from %d %d %d %d %d to %d %d %d %d %d ] ",comm->blocks[1][i].dim.i,comm->blocks[1][i].dim.j,comm->blocks[1][i].dim.k,comm->blocks[1][i].read.i,comm->blocks[1][i].read.j,comm->blocks[1][i].read.k,comm->blocks[1][i].read.jStride,comm->blocks[1][i].read.kStride,comm->blocks[1][i].write.i,comm->blocks[1][i].write.j,comm->blocks[1][i].write.k,comm->blocks[1][i].write.jStride,comm->blocks[1][i].write.kStride);
  printf("\n");
  }
  if(printSendRecv & 0x4){
  printf("num_recvs=%2d ",comm->num_recvs);
  printf("recv_ranks=[ ");for(i=0;i<comm->num_recvs;i++)printf("%2d ",comm->recv_ranks[i]);printf("] ");
  printf("recv_sizes=[ ");for(i=0;i<comm->num_recvs;i++)printf("%2d ",comm->recv_sizes[i]);printf("] ");
  printf("recv_buffers=[ ");for(i=0;i<comm->num_recvs;i++)printf("%08lx ",(uint64_t)comm->recv_buffers[i]);printf("] ");
  for(i=0;i<comm->num_blocks[2];i++)printf("[ %dx%dx%d from %d %d %d %d %d to %d %d %d %d %d ] ",comm->blocks[2][i].dim.i,comm->blocks[2][i].dim.j,comm->blocks[2][i].dim.k,comm->blocks[2][i].read.i,comm->blocks[2][i].read.j,comm->blocks[2][i].read.k,comm->blocks[2][i].read.jStride,comm->blocks[2][i].read.kStride,comm->blocks[2][i].write.i,comm->blocks[2][i].write.j,comm->blocks[2][i].write.k,comm->blocks[2][i].write.jStride,comm->blocks[2][i].write.kStride);
  printf("\n");
  }
  fflush(stdout);
}
//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  int sendRank;
  int sendBoxID;
  int sendBox;
  int sendDir;
  int recvRank;
  int recvBoxID;
  int recvBox;
} GZ_type;


int qsortGZ(const void *a, const void*b){
  GZ_type *gza = (GZ_type*)a;
  GZ_type *gzb = (GZ_type*)b;
  // by convention, MPI buffers are first sorted by sendRank
  if(gza->sendRank < gzb->sendRank)return(-1);
  if(gza->sendRank > gzb->sendRank)return( 1);
  // then by sendBoxID
  if(gza->sendBoxID < gzb->sendBoxID)return(-1);
  if(gza->sendBoxID > gzb->sendBoxID)return( 1);
  // and finally by the direction sent
  if(gza->sendDir < gzb->sendDir)return(-1);
  if(gza->sendDir > gzb->sendDir)return( 1);
  return(0);
}


int qsortInt(const void *a, const void *b){
  int *ia = (int*)a;
  int *ib = (int*)b;
  if(*ia < *ib)return(-1);
  if(*ia > *ib)return( 1);
               return( 0);
}

int qsortBlock(const void *a, const void *b){
  blockCopy_type *ba = (blockCopy_type*)a;
  blockCopy_type *bb = (blockCopy_type*)b;

  if(ba->write.box >= 0){
    // sort by box...
    if(ba->write.box < bb->write.box)return(-1);
    if(ba->write.box > bb->write.box)return( 1);
    // now sort by k
    if(ba->write.k   < bb->write.k  )return(-1);
    if(ba->write.k   > bb->write.k  )return( 1);
    // now sort by j
    if(ba->write.j   < bb->write.j  )return(-1);
    if(ba->write.j   > bb->write.j  )return( 1);
    // now sort by i
    if(ba->write.i   < bb->write.i  )return(-1);
    if(ba->write.i   > bb->write.i  )return( 1);
  }else if(ba->read.box >= 0){
    // sort by box...
    if(ba->read.box  < bb->read.box )return(-1);
    if(ba->read.box  > bb->read.box )return( 1);
    // now sort by k
    if(ba->read.k    < bb->read.k   )return(-1);
    if(ba->read.k    > bb->read.k   )return( 1);
    // now sort by j
    if(ba->read.j    < bb->read.j   )return(-1);
    if(ba->read.j    > bb->read.j   )return( 1);
    // now sort by i
    if(ba->read.i    < bb->read.i   )return(-1);
    if(ba->read.i    > bb->read.i   )return( 1);
  }
                                     return( 0);
}


//------------------------------------------------------------------------------------------------------------------------------
void decompose_level_lex(int *rank_of_box, int idim, int jdim, int kdim, int ranks){
  // simple lexicographical decomposition of the domain (i-j-k ordering)
  // load balancing is easily realized
  // unfortunately, each process will likely receive one or two long pencils of boxes. 
  // as such, the resultant surface:volum ratio will likely be poor
  int boxes = idim*jdim*kdim;
  int i,j,k;
  for(k=0;k<kdim;k++){
  for(j=0;j<jdim;j++){
  for(i=0;i<idim;i++){
    int b = k*jdim*idim + j*idim + i;
    rank_of_box[b] = ((uint64_t)ranks*(uint64_t)b)/(uint64_t)boxes; // ranks*b can be as larger than ranks^2... can over flow int
  }}} 
}


//---------------------------------------------------------------------------------------------------------------------------------------------------
void decompose_level_bisection_special(int *rank_of_box, int jStride, int kStride, int ilo, int jlo, int klo, int idim, int jdim, int kdim, int rank_lo, int ranks){
  // if possible, recursively partition the domain by a prime number (e.g. try an parition a 9^3 array into 3 equal pieces instead of 5x9^2 and 4x9^2)
  // if not, default to simple bisection
  // this function should ensure that each process receives a compact rectahedral collection of boxes
  // however, load imbalance can occur
  // the choice of whether to try and partition with the largest prime or smallest prime first is up to the user

  #define numPrimes 13
  //int primes[numPrimes] = {41,37,31,29,23,19,17,13,11,7,5,3,2};
  int primes[numPrimes] = {2,3,5,7,11,13,17,19,23,29,31,37,41};
  int i,j,k,p,f,ff;


  // base case, no further recursion...
  if( (ranks==1)|| ((idim==1)&&(jdim==1)&&(kdim==1)) ){
    for(i=ilo;i<ilo+idim;i++){
    for(j=jlo;j<jlo+jdim;j++){
    for(k=klo;k<klo+kdim;k++){
      int b = i + j*jStride + k*kStride;
      rank_of_box[b] = rank_lo;
    }}}
    return;
  }


  // special cases for perfectly matched problem sizes with numbers of processes (but not powers of 2)...
  for(p=0;p<numPrimes;p++){
    f=primes[p];
    if( (kdim>=idim)&&(kdim>=jdim) ){if( (kdim%f==0) && (ranks%f==0) ){for(ff=0;ff<f;ff++)decompose_level_bisection_special(rank_of_box,jStride,kStride,ilo,jlo,klo+ff*kdim/f,idim,jdim,kdim/f,rank_lo+ff*ranks/f,ranks/f);return;}}
    if( (jdim>=idim)&&(jdim>=kdim) ){if( (jdim%f==0) && (ranks%f==0) ){for(ff=0;ff<f;ff++)decompose_level_bisection_special(rank_of_box,jStride,kStride,ilo,jlo+ff*jdim/f,klo,idim,jdim/f,kdim,rank_lo+ff*ranks/f,ranks/f);return;}}
    if( (idim>=jdim)&&(idim>=kdim) ){if( (idim%f==0) && (ranks%f==0) ){for(ff=0;ff<f;ff++)decompose_level_bisection_special(rank_of_box,jStride,kStride,ilo+ff*idim/f,jlo,klo,idim/f,jdim,kdim,rank_lo+ff*ranks/f,ranks/f);return;}}
  }


  // try and bisect the domain in the i-dimension
  if( (idim>=jdim)&&(idim>=kdim) ){
    int dim0 = (int)(0.5*(double)idim + 0.50);
    int dim1 = idim-dim0;
    int r0 = (int)( 0.5 + (double)ranks*(double)dim0/(double)idim );
    int r1 = ranks-r0;
    decompose_level_bisection_special(rank_of_box,jStride,kStride,ilo     ,jlo,klo,dim0,jdim,kdim,rank_lo   ,r0); // lo
    decompose_level_bisection_special(rank_of_box,jStride,kStride,ilo+dim0,jlo,klo,dim1,jdim,kdim,rank_lo+r0,r1); // hi
    return;
  }
  // try and bisect the domain in the j-dimension
  if( (jdim>=idim)&&(jdim>=kdim) ){
    int dim0 = (int)(0.5*(double)jdim + 0.50);
    int dim1 = jdim-dim0;
    int r0 = (int)( 0.5 + (double)ranks*(double)dim0/(double)jdim );
    int r1 = ranks-r0;
    decompose_level_bisection_special(rank_of_box,jStride,kStride,ilo,jlo     ,klo,idim,dim0,kdim,rank_lo   ,r0); // lo
    decompose_level_bisection_special(rank_of_box,jStride,kStride,ilo,jlo+dim0,klo,idim,dim1,kdim,rank_lo+r0,r1); // hi
    return;
  }
  // try and bisect the domain in the k-dimension
  if( (kdim>=idim)&&(kdim>=jdim) ){
    int dim0 = (int)(0.5*(double)kdim + 0.50);
    int dim1 = kdim-dim0;
    int r0 = (int)( 0.5 + (double)ranks*(double)dim0/(double)kdim );
    int r1 = ranks-r0;
    decompose_level_bisection_special(rank_of_box,jStride,kStride,ilo,jlo,klo     ,idim,jdim,dim0,rank_lo   ,r0); // lo
    decompose_level_bisection_special(rank_of_box,jStride,kStride,ilo,jlo,klo+dim0,idim,jdim,dim1,rank_lo+r0,r1); // hi
    return;
  }
  fprintf(stderr,"decompose_level_bisection_special failed !!!\n");exit(0);
}


//---------------------------------------------------------------------------------------------------------------------------------------------------
void decompose_level_bisection(int *rank_of_box, int jStride, int kStride, int ilo, int jlo, int klo, int idim, int jdim, int kdim, int ranks, int sfc_offset, int sfc_max_length){

  // base case... 
  if( (idim==1) && (jdim==1) && (kdim==1) ){
    int b = ilo + jlo*jStride + klo*kStride;
    rank_of_box[b] = ((uint64_t)ranks*(uint64_t)sfc_offset)/(uint64_t)sfc_max_length; // sfc_max_length is the precomputed maximum length
    return;
  }

  // try and bisect the domain in the i-dimension
  if( (idim>=jdim)&&(idim>=kdim) ){
    int dim0 = (int)(0.5*(double)idim + 0.50);
    int dim1 = idim-dim0;
    int sfc_delta = dim0*jdim*kdim;
    decompose_level_bisection(rank_of_box,jStride,kStride,ilo     ,jlo,klo,dim0,jdim,kdim,ranks,sfc_offset          ,sfc_max_length); // lo
    decompose_level_bisection(rank_of_box,jStride,kStride,ilo+dim0,jlo,klo,dim1,jdim,kdim,ranks,sfc_offset+sfc_delta,sfc_max_length); // hi
    return;
  }

  // try and bisect the domain in the j-dimension
  if( (jdim>=idim)&&(jdim>=kdim) ){
    int dim0 = (int)(0.5*(double)jdim + 0.50);
    int dim1 = jdim-dim0;
    int sfc_delta = idim*dim0*kdim;
    decompose_level_bisection(rank_of_box,jStride,kStride,ilo,jlo     ,klo,idim,dim0,kdim,ranks,sfc_offset          ,sfc_max_length); // lo
    decompose_level_bisection(rank_of_box,jStride,kStride,ilo,jlo+dim0,klo,idim,dim1,kdim,ranks,sfc_offset+sfc_delta,sfc_max_length); // hi
    return;
  }

  // try and bisect the domain in the k-dimension
  if( (kdim>=idim)&&(kdim>=jdim) ){
    int dim0 = (int)(0.5*(double)kdim + 0.50);
    int dim1 = kdim-dim0;
    int sfc_delta = idim*jdim*dim0;
    decompose_level_bisection(rank_of_box,jStride,kStride,ilo,jlo,klo     ,idim,jdim,dim0,ranks,sfc_offset          ,sfc_max_length); // lo
    decompose_level_bisection(rank_of_box,jStride,kStride,ilo,jlo,klo+dim0,idim,jdim,dim1,ranks,sfc_offset+sfc_delta,sfc_max_length); // hi
    return;
  }

  // failure...
  fprintf(stderr,"decompose_level_bisection failed !!!\n");exit(0);
}


//---------------------------------------------------------------------------------------------------------------------------------------------------
// Given a bounding box (idim,jdim,kdim) use a Z-morton Space Filling Curve (SFC) to assign the boxes within the (boxes_in_i,boxes_in_j,boxes_in_k) valid region domain
//  sfc_offset is the current offset within the space filling curve (starts with 0)
//  this function returns the new offset based on how many actual boxes it found within (ilo,jlo,klo) + (idim,jdim,kdim)
//  sfc_max_length is the maximum length of the SFC.  Note, if this length exceeds boxes_in_i*boxes_in_j*boxes_in_k, then some processes with receive no work
int decompose_level_zmort(int *rank_of_box, int boxes_in_i, int boxes_in_j, int boxes_in_k, int ilo, int jlo, int klo, int idim, int jdim, int kdim, int ranks, int sfc_offset, int sfc_max_length){

  // invalid cases...
  if(idim<1)return(sfc_offset);
  if(jdim<1)return(sfc_offset);
  if(kdim<1)return(sfc_offset);
  if(ilo <0)return(sfc_offset);
  if(jlo <0)return(sfc_offset);
  if(klo <0)return(sfc_offset);

  // base case... 
  if( (idim==1) && (jdim==1) && (kdim==1) ){
    if( (ilo<boxes_in_i) && (jlo<boxes_in_j) && (klo<boxes_in_k) ){
      // deemed a valid box (could be augmented for irregular domains)
      int b = ilo + jlo*boxes_in_i + klo*boxes_in_i*boxes_in_j;
      rank_of_box[b] = ((uint64_t)ranks*(uint64_t)(sfc_offset))/(uint64_t)sfc_max_length; // sfc_max_length is the precomputed maximum length
      return(sfc_offset+1);
    }
    return(sfc_offset); // region outside valid domain;  sfc_offset is unchanged
  }

  // bisect in 3D...
  int imid = ilo + (idim/2);
  int jmid = jlo + (jdim/2);
  int kmid = klo + (kdim/2);

  sfc_offset=decompose_level_zmort(rank_of_box,boxes_in_i,boxes_in_j,boxes_in_k,ilo ,jlo ,klo ,     idim/2,     jdim/2,     kdim/2,ranks,sfc_offset,sfc_max_length);
  sfc_offset=decompose_level_zmort(rank_of_box,boxes_in_i,boxes_in_j,boxes_in_k,imid,jlo ,klo ,idim-idim/2,     jdim/2,     kdim/2,ranks,sfc_offset,sfc_max_length);
  sfc_offset=decompose_level_zmort(rank_of_box,boxes_in_i,boxes_in_j,boxes_in_k,ilo ,jmid,klo ,     idim/2,jdim-jdim/2,     kdim/2,ranks,sfc_offset,sfc_max_length);
  sfc_offset=decompose_level_zmort(rank_of_box,boxes_in_i,boxes_in_j,boxes_in_k,imid,jmid,klo ,idim-idim/2,jdim-jdim/2,     kdim/2,ranks,sfc_offset,sfc_max_length);
  sfc_offset=decompose_level_zmort(rank_of_box,boxes_in_i,boxes_in_j,boxes_in_k,ilo ,jlo ,kmid,     idim/2,     jdim/2,kdim-kdim/2,ranks,sfc_offset,sfc_max_length);
  sfc_offset=decompose_level_zmort(rank_of_box,boxes_in_i,boxes_in_j,boxes_in_k,imid,jlo ,kmid,idim-idim/2,     jdim/2,kdim-kdim/2,ranks,sfc_offset,sfc_max_length);
  sfc_offset=decompose_level_zmort(rank_of_box,boxes_in_i,boxes_in_j,boxes_in_k,ilo ,jmid,kmid,     idim/2,jdim-jdim/2,kdim-kdim/2,ranks,sfc_offset,sfc_max_length);
  sfc_offset=decompose_level_zmort(rank_of_box,boxes_in_i,boxes_in_j,boxes_in_k,imid,jmid,kmid,idim-idim/2,jdim-jdim/2,kdim-kdim/2,ranks,sfc_offset,sfc_max_length);
  return(sfc_offset);
}


//------------------------------------------------------------------------------------------------------------------------------
//int decompose_level_hilbert(int *rank_of_box, int boxes_in_i, int boxes_in_j, int boxes_in_k, int ilo, int jlo, int klo, int idim, int jdim, int kdim, int ranks, int sfc_offset, int sfc_max_length){
// implements a 3D hilbert curve on the non-power of two domain using a power of two bounding box
//}


//---------------------------------------------------------------------------------------------------------------------------------------------------
void print_decomposition(level_type *level){
  if(level->my_rank!=0)return;
  printf("\n");
  int i,j,k;
  int jStride = level->boxes_in.i;
  int kStride = level->boxes_in.i*level->boxes_in.j;
  for(k=level->boxes_in.k-1;k>=0;k--){ // (i,j,k)=(0,0,0) is bottom left corner
  for(j=level->boxes_in.j-1;j>=0;j--){ // (i,j)=(0,0) is bottom left corner
  for(i=0;i<j;i++)printf(" ");
  for(i=0;i<level->boxes_in.i;i++){
    int b = i + j*jStride + k*kStride;
    printf("%4d ",level->rank_of_box[b]);
  }printf("\n");
  }printf("\n\n");
  }
  fflush(stdout);
}


//------------------------------------------------------------------------------------------------------------------------------
// append the specified block (logical region) to the current list of blocks
// each block may be tiled to...
//  - create more parallelism across the list of blocks
//  - limit parallelism within a block
//  - limit the memory requirements for each block
#ifndef BLOCK_LIST_MIN_SIZE
#define BLOCK_LIST_MIN_SIZE 1000
#endif
void append_block_to_list(blockCopy_type ** blocks, int *allocated_blocks, int *num_blocks,
                          int dim_i, int dim_j, int dim_k,
                          int  read_box, double*  read_ptr, int  read_i, int  read_j, int  read_k, int  read_jStride, int  read_kStride, int  read_scale,
                          int write_box, double* write_ptr, int write_i, int write_j, int write_k, int write_jStride, int write_kStride, int write_scale,
                          int blockcopy_tile_i, int blockcopy_tile_j, int blockcopy_tile_k, 
                          int subtype
                         ){
  // Take a dim_j x dim_k iteration space and tile it into smaller faces of size blockcopy_tile_j x blockcopy_tile_k
  // This increases the number of blockCopies in the ghost zone exchange and thereby increases the thread-level parallelism

  #if 0
  // use recursive (z-mort) ordering of tiles in order to improve locality on deep memory hierarchies...
  int doRecursion=0;
  if(dim_i > blockcopy_tile_i)doRecursion=1;
  if(dim_j > blockcopy_tile_j)doRecursion=1;
  if(dim_k > blockcopy_tile_k)doRecursion=1;
  if( read_scale != 1)doRecursion=0; // disable recursion for restriction
  if(write_scale != 1)doRecursion=0; // disable recursion for interpolation
  if(doRecursion){
    int mid_i = (dim_i + 1)/2;
    int mid_j = (dim_j + 1)/2;
    int mid_k = (dim_k + 1)/2;
        mid_i = blockcopy_tile_i*( (mid_i+blockcopy_tile_i-1)/blockcopy_tile_i);
        mid_j = blockcopy_tile_j*( (mid_j+blockcopy_tile_j-1)/blockcopy_tile_j);
        mid_k = blockcopy_tile_k*( (mid_k+blockcopy_tile_k-1)/blockcopy_tile_k);
     if(mid_i>dim_i)mid_i=dim_i;
     if(mid_j>dim_j)mid_j=dim_j;
     if(mid_k>dim_k)mid_k=dim_k;
    append_block_to_list(blocks,allocated_blocks,num_blocks,      mid_i,      mid_j,      mid_k,
                          read_box, read_ptr, read_i      , read_j      , read_k      , read_jStride, read_kStride, read_scale,
                         write_box,write_ptr,write_i      ,write_j      ,write_k      ,write_jStride,write_kStride,write_scale,
                         blockcopy_tile_i,blockcopy_tile_j,blockcopy_tile_k,subtype);
    append_block_to_list(blocks,allocated_blocks,num_blocks,dim_i-mid_i,      mid_j,      mid_k,
                          read_box, read_ptr, read_i+mid_i, read_j      , read_k      , read_jStride, read_kStride, read_scale,
                         write_box,write_ptr,write_i+mid_i,write_j      ,write_k      ,write_jStride,write_kStride,write_scale,
                         blockcopy_tile_i,blockcopy_tile_j,blockcopy_tile_k,subtype);
    append_block_to_list(blocks,allocated_blocks,num_blocks,      mid_i,dim_j-mid_j,      mid_k,
                          read_box, read_ptr, read_i      , read_j+mid_j, read_k      , read_jStride, read_kStride, read_scale,
                         write_box,write_ptr,write_i      ,write_j+mid_j,write_k      ,write_jStride,write_kStride,write_scale,
                         blockcopy_tile_i,blockcopy_tile_j,blockcopy_tile_k,subtype);
    append_block_to_list(blocks,allocated_blocks,num_blocks,dim_i-mid_i,dim_j-mid_j,      mid_k,
                          read_box, read_ptr, read_i+mid_i, read_j+mid_j, read_k      , read_jStride, read_kStride, read_scale,
                         write_box,write_ptr,write_i+mid_i,write_j+mid_j,write_k      ,write_jStride,write_kStride,write_scale,
                         blockcopy_tile_i,blockcopy_tile_j,blockcopy_tile_k,subtype);
    append_block_to_list(blocks,allocated_blocks,num_blocks,      mid_i,      mid_j,dim_k-mid_k,
                          read_box, read_ptr, read_i      , read_j      , read_k+mid_k, read_jStride, read_kStride, read_scale,
                         write_box,write_ptr,write_i      ,write_j      ,write_k+mid_k,write_jStride,write_kStride,write_scale,
                         blockcopy_tile_i,blockcopy_tile_j,blockcopy_tile_k,subtype);
    append_block_to_list(blocks,allocated_blocks,num_blocks,dim_i-mid_i,      mid_j,dim_k-mid_k,
                          read_box, read_ptr, read_i+mid_i, read_j      , read_k+mid_k, read_jStride, read_kStride, read_scale,
                         write_box,write_ptr,write_i+mid_i,write_j      ,write_k+mid_k,write_jStride,write_kStride,write_scale,
                         blockcopy_tile_i,blockcopy_tile_j,blockcopy_tile_k,subtype);
    append_block_to_list(blocks,allocated_blocks,num_blocks,      mid_i,dim_j-mid_j,dim_k-mid_k,
                          read_box, read_ptr, read_i      , read_j+mid_j, read_k+mid_k, read_jStride, read_kStride, read_scale,
                         write_box,write_ptr,write_i      ,write_j+mid_j,write_k+mid_k,write_jStride,write_kStride,write_scale,
                         blockcopy_tile_i,blockcopy_tile_j,blockcopy_tile_k,subtype);
    append_block_to_list(blocks,allocated_blocks,num_blocks,dim_i-mid_i,dim_j-mid_j,dim_k-mid_k,
                          read_box, read_ptr, read_i+mid_i, read_j+mid_j, read_k+mid_k, read_jStride, read_kStride, read_scale,
                         write_box,write_ptr,write_i+mid_i,write_j+mid_j,write_k+mid_k,write_jStride,write_kStride,write_scale,
                         blockcopy_tile_i,blockcopy_tile_j,blockcopy_tile_k,subtype);
    return;
  }
  #endif
  // read_/write_scale are used to stride appropriately when read and write loop iterations spaces are different 
  // ghostZone:     read_scale=1, write_scale=1
  // interpolation: read_scale=1, write_scale=2
  // restriction:   read_scale=2, write_scale=1
  // FIX... dim_i,j,k -> read_dim_i,j,k, write_dim_i,j,k
  int ii,jj,kk;
  for(kk=0;kk<dim_k;kk+=blockcopy_tile_k){
  for(jj=0;jj<dim_j;jj+=blockcopy_tile_j){
  for(ii=0;ii<dim_i;ii+=blockcopy_tile_i){
    int dim_k_mod = dim_k-kk;if(dim_k_mod>blockcopy_tile_k)dim_k_mod=blockcopy_tile_k;
    int dim_j_mod = dim_j-jj;if(dim_j_mod>blockcopy_tile_j)dim_j_mod=blockcopy_tile_j;
    int dim_i_mod = dim_i-ii;if(dim_i_mod>blockcopy_tile_i)dim_i_mod=blockcopy_tile_i;
    if(*num_blocks >= *allocated_blocks){
      int oldSize = *allocated_blocks;
      if(*allocated_blocks == 0){*allocated_blocks=BLOCK_LIST_MIN_SIZE;*blocks=(blockCopy_type*) malloc(                 (*allocated_blocks)*sizeof(blockCopy_type));}
                            else{*allocated_blocks*=2;                 *blocks=(blockCopy_type*)realloc((void*)(*blocks),(*allocated_blocks)*sizeof(blockCopy_type));}
      if(*blocks == NULL){fprintf(stderr,"realloc failed - append_block_to_list (%d -> %d)\n",oldSize,*allocated_blocks);exit(0);}
    }
    (*blocks)[*num_blocks].subtype       = subtype;
    (*blocks)[*num_blocks].dim.i         = dim_i_mod;
    (*blocks)[*num_blocks].dim.j         = dim_j_mod;
    (*blocks)[*num_blocks].dim.k         = dim_k_mod;
    (*blocks)[*num_blocks].read.box      = read_box;
    (*blocks)[*num_blocks].read.ptr      = read_ptr;
    (*blocks)[*num_blocks].read.i        = read_i + read_scale*ii;
    (*blocks)[*num_blocks].read.j        = read_j + read_scale*jj;
    (*blocks)[*num_blocks].read.k        = read_k + read_scale*kk;
    (*blocks)[*num_blocks].read.jStride  = read_jStride;
    (*blocks)[*num_blocks].read.kStride  = read_kStride;
    (*blocks)[*num_blocks].write.box     = write_box;
    (*blocks)[*num_blocks].write.ptr     = write_ptr;
    (*blocks)[*num_blocks].write.i       = write_i + write_scale*ii;
    (*blocks)[*num_blocks].write.j       = write_j + write_scale*jj;
    (*blocks)[*num_blocks].write.k       = write_k + write_scale*kk;
    (*blocks)[*num_blocks].write.jStride = write_jStride;
    (*blocks)[*num_blocks].write.kStride = write_kStride;
             (*num_blocks)++;
  }}}
}


//----------------------------------------------------------------------------------------------------------------------------------------------------
// create a mini program that traverses the domain boundary intersecting with this process's boxes
// This includes faces, corners, and edges
void build_boundary_conditions(level_type *level, int shape){
  level->boundary_condition.blocks[shape]           = NULL;	// default for periodic (i.e. no BC's)
  level->boundary_condition.num_blocks[shape]       = 0;	// default for periodic (i.e. no BC's)
  level->boundary_condition.allocated_blocks[shape] = 0;	// default for periodic (i.e. no BC's)
  if(level->boundary_condition.type == BC_PERIODIC)return;

//int    faces[27] = {0,0,0,0,1,0,0,0,0,  0,1,0,1,0,1,0,1,0,  0,0,0,0,1,0,0,0,0};
  int    edges[27] = {0,1,0,1,0,1,0,1,0,  1,0,1,0,0,0,1,0,1,  0,1,0,1,0,1,0,1,0};
  int  corners[27] = {1,0,1,0,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,  1,0,1,0,0,0,1,0,1};

  int box, di,dj,dk;
  for(box=0;box<level->num_my_boxes;box++){	// traverse my list of boxes...
  for(dk=-1;dk<=1;dk++){			// for each box, examine its 26 neighbors...
  for(dj=-1;dj<=1;dj++){
  for(di=-1;di<=1;di++){
    int dir = 13+di+3*dj+9*dk; // face/edge/corner of *THIS* box (not the domain)

    // determine if this region (box's di,dj,dk ghost zone) is outside of the domain
    int regionIsOutside=0;
    int normal = 13; // normal effectively defines the normal vector to the *DOMAIN* for this region... 
                     // this addition is necessary for linearly interpolated BC's as a box's corner is not necessarily a domain's corner
    int myBox_i = level->my_boxes[box].low.i / level->box_dim;
    int myBox_j = level->my_boxes[box].low.j / level->box_dim;
    int myBox_k = level->my_boxes[box].low.k / level->box_dim;
    int neighborBox_i = (  myBox_i + di );
    int neighborBox_j = (  myBox_j + dj );
    int neighborBox_k = (  myBox_k + dk );
    if( neighborBox_i < 0                 ){regionIsOutside=1;normal-=1;}
    if( neighborBox_j < 0                 ){regionIsOutside=1;normal-=3;}
    if( neighborBox_k < 0                 ){regionIsOutside=1;normal-=9;}
    if( neighborBox_i >=level->boxes_in.i ){regionIsOutside=1;normal+=1;}
    if( neighborBox_j >=level->boxes_in.j ){regionIsOutside=1;normal+=3;}
    if( neighborBox_k >=level->boxes_in.k ){regionIsOutside=1;normal+=9;}

    // calculate ghost zone region size and coordinates relative to the first non-ghost zone element (0,0,0)
    int block_i=-1,block_j=-1,block_k=-1;
    int   dim_i=-1,  dim_j=-1,  dim_k=-1;
    switch(di){
      case -1:dim_i=level->box_ghosts;block_i=0-level->box_ghosts;break;
      case  0:dim_i=level->box_dim;   block_i=0;                  break;
      case  1:dim_i=level->box_ghosts;block_i=0+level->box_dim;   break;
    }
    switch(dj){
      case -1:dim_j=level->box_ghosts;block_j=0-level->box_ghosts;break;
      case  0:dim_j=level->box_dim;   block_j=0;                  break;
      case  1:dim_j=level->box_ghosts;block_j=0+level->box_dim;   break;
    }
    switch(dk){
      case -1:dim_k=level->box_ghosts;block_k=0-level->box_ghosts;break;
      case  0:dim_k=level->box_dim;   block_k=0;                  break;
      case  1:dim_k=level->box_ghosts;block_k=0+level->box_dim;   break;
    }

    // use regionIsOutside to short circuit logic and cull unnecessary regions...
    switch(shape){
      case STENCIL_SHAPE_STAR:      if(edges[dir]||corners[dir])regionIsOutside=0;break; // star-shaped stencils don't need BC's enforced on corners or edges
      case STENCIL_SHAPE_NO_CORNERS:if(            corners[dir])regionIsOutside=0;break; // these stencils don't need BC's enforced on edges
    }

    // default tile sizes...
    // NOTE, BC's may never tile smaller than the ghost zone depth
    int blockcopy_i = (BLOCKCOPY_TILE_I < level->box_ghosts) ? level->box_ghosts : BLOCKCOPY_TILE_I;
    int blockcopy_j = (BLOCKCOPY_TILE_J < level->box_ghosts) ? level->box_ghosts : BLOCKCOPY_TILE_J;
    int blockcopy_k = (BLOCKCOPY_TILE_K < level->box_ghosts) ? level->box_ghosts : BLOCKCOPY_TILE_K;

    #if 0
    // 2D tiling of faces
    // 1D tiling of edges
    // corners use defaults
    switch(dir){
      case  1:blockcopy_i=    8;blockcopy_j=10000;blockcopy_k=10000;break; //  i edge
      case  3:blockcopy_i=10000;blockcopy_j=    8;blockcopy_k=10000;break; //  j edge
      case  4:blockcopy_i=    8;blockcopy_j=    8;blockcopy_k=10000;break; // ij face
      case  5:blockcopy_i=10000;blockcopy_j=    8;blockcopy_k=10000;break; //  j edge
      case  7:blockcopy_i=    8;blockcopy_j=10000;blockcopy_k=10000;break; //  i edge

      case  9:blockcopy_i=10000;blockcopy_j=10000;blockcopy_k=    8;break; //  k edge
      case 10:blockcopy_i=    8;blockcopy_j=10000;blockcopy_k=    8;break; // ik face
      case 11:blockcopy_i=10000;blockcopy_j=10000;blockcopy_k=    8;break; //  k edge
      case 12:blockcopy_i=10000;blockcopy_j=    8;blockcopy_k=    8;break; // jk face

      case 14:blockcopy_i=10000;blockcopy_j=    8;blockcopy_k=    8;break; // jk face
      case 15:blockcopy_i=10000;blockcopy_j=10000;blockcopy_k=    8;break; //  k edge
      case 16:blockcopy_i=    8;blockcopy_j=10000;blockcopy_k=    8;break; // ik face
      case 17:blockcopy_i=10000;blockcopy_j=10000;blockcopy_k=    8;break; //  k edge

      case 19:blockcopy_i=    8;blockcopy_j=10000;blockcopy_k=10000;break; //  i edge
      case 21:blockcopy_i=10000;blockcopy_j=    8;blockcopy_k=10000;break; //  j edge
      case 22:blockcopy_i=    8;blockcopy_j=    8;blockcopy_k=10000;break; // ij face
      case 23:blockcopy_i=10000;blockcopy_j=    8;blockcopy_k=10000;break; //  j edge
      case 25:blockcopy_i=    8;blockcopy_j=10000;blockcopy_k=10000;break; //  i edge
    }
    #endif

    if(regionIsOutside){
    append_block_to_list(&(level->boundary_condition.blocks[shape]),&(level->boundary_condition.allocated_blocks[shape]),&(level->boundary_condition.num_blocks[shape]),
      /* dim.i         = */ dim_i,
      /* dim.j         = */ dim_j,
      /* dim.k         = */ dim_k,
      /* read.box      = */ box,
      /* read.ptr      = */ NULL,
      /* read.i        = */ block_i,
      /* read.j        = */ block_j,
      /* read.k        = */ block_k,
      /* read.jStride  = */ level->my_boxes[box].jStride,
      /* read.kStride  = */ level->my_boxes[box].kStride,
      /* read.scale    = */ 1,
      /* write.box     = */ box,
      /* write.ptr     = */ NULL,
      /* write.i       = */ block_i,
      /* write.j       = */ block_j,
      /* write.k       = */ block_k,
      /* write.jStride = */ level->my_boxes[box].jStride,
      /* write.kStride = */ level->my_boxes[box].kStride,
      /* write.scale   = */ 1,
      /* blockcopy_i   = */ blockcopy_i,
      /* blockcopy_j   = */ blockcopy_j,
      /* blockcopy_k   = */ blockcopy_k,
      /* subtype       = */ normal
    );
  }}}}}

  #ifdef BLOCK_SPATIAL_SORT
  // sort all the resultant blocks by box,k,j,i (good locality)
  qsort(level->boundary_condition.blocks[shape],level->boundary_condition.num_blocks[shape],sizeof(blockCopy_type),qsortBlock);
  #endif
}

//----------------------------------------------------------------------------------------------------------------------------------------------------
// create a mini program that packs data into MPI recv buffers, exchanges local data, and unpacks the MPI send buffers
//   broadly speaking... 
//   1. traverse my list of Boxes and create a list of ghosts that must be sent
//   2. create a list of neighbors to send to
//   3. allocate and populate the pack list and allocate the send buffers
//   4. allocate and populate the local exchange list
//   5. traverse my list of Boxes and create a list of ghosts that must be received
//   6. create a list of neighbors to receive from
//   7. allocate and populate the unpack list and allocate the recv buffers
//
//   thus a ghost zone exchange is
//   1. prepost a Irecv for each MPI recv buffer (1 per neighbor)
//   2. traverse the pack list
//   3. post the Isends for each MPI send buffer (1 per neighbor)
//   4. traverse the local copy list
//   5. waitall
//   6. traverse the unpack list
//
//     / 24 25 26 /
//    / 21 22 23 /	(k+1)
//   / 18 19 20 /
//
//     / 15 16 17 /
//    / 12 13 14 /	(k)
//   /  9 10 11 /
//
//     /  6  7  8 /
//    /  3  4  5 /	(k-1)
//   /  0  1  2 /
//
void build_exchange_ghosts(level_type *level, int shape){
  int    faces[27] = {0,0,0,0,1,0,0,0,0,  0,1,0,1,0,1,0,1,0,  0,0,0,0,1,0,0,0,0};
  int    edges[27] = {0,1,0,1,0,1,0,1,0,  1,0,1,0,0,0,1,0,1,  0,1,0,1,0,1,0,1,0};
  int  corners[27] = {1,0,1,0,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,  1,0,1,0,0,0,1,0,1};

  // initialize to defaults...
  level->exchange_ghosts[shape].num_recvs           = 0;
  level->exchange_ghosts[shape].num_sends           = 0;
  level->exchange_ghosts[shape].recv_ranks          = NULL;
  level->exchange_ghosts[shape].send_ranks          = NULL;
  level->exchange_ghosts[shape].recv_sizes          = NULL;
  level->exchange_ghosts[shape].send_sizes          = NULL;
  level->exchange_ghosts[shape].recv_buffers        = NULL;
  level->exchange_ghosts[shape].send_buffers        = NULL;
  level->exchange_ghosts[shape].blocks[0]           = NULL;
  level->exchange_ghosts[shape].blocks[1]           = NULL;
  level->exchange_ghosts[shape].blocks[2]           = NULL;
  level->exchange_ghosts[shape].num_blocks[0]       = 0;
  level->exchange_ghosts[shape].num_blocks[1]       = 0;
  level->exchange_ghosts[shape].num_blocks[2]       = 0;
  level->exchange_ghosts[shape].allocated_blocks[0] = 0;
  level->exchange_ghosts[shape].allocated_blocks[1] = 0;
  level->exchange_ghosts[shape].allocated_blocks[2] = 0;
  #ifdef USE_MPI
  level->exchange_ghosts[shape].requests            = NULL;
  level->exchange_ghosts[shape].status              = NULL;
  #endif

  int    n,CommunicateThisDir[27];for(n=0;n<27;n++)CommunicateThisDir[n] = faces[n] + edges[n] + corners[n];// to be safe, communicate everything
  switch(shape){
    case STENCIL_SHAPE_BOX:       for(n=0;n<27;n++)CommunicateThisDir[n] = faces[n] + edges[n] + corners[n];break;
    case STENCIL_SHAPE_STAR:      for(n=0;n<27;n++)CommunicateThisDir[n] = faces[n]                        ;break;
    case STENCIL_SHAPE_NO_CORNERS:for(n=0;n<27;n++)CommunicateThisDir[n] = faces[n] + edges[n]             ;break;
  }

  int sendBox,recvBox;
  int stage;
  int _rank;
  int ghost,numGhosts,numGhostsRemote;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // traverse my list of boxes and create a lists of neighboring boxes and neighboring ranks
  GZ_type *ghostsToSend = (GZ_type*)malloc(26*level->num_my_boxes*sizeof(GZ_type)); // There are at most 26 neighbors per box.
         int *sendRanks = (    int*)malloc(26*level->num_my_boxes*sizeof(    int)); // There are at most 26 neighbors per box.
  if(level->num_my_boxes>0){
  if(ghostsToSend == NULL){fprintf(stderr,"malloc failed - build_exchange_ghosts/ghostsToSend\n");exit(0);}
  if(sendRanks    == NULL){fprintf(stderr,"malloc failed - build_exchange_ghosts/sendRanks   \n");exit(0);}
  }
  numGhosts       = 0;
  numGhostsRemote = 0;
  for(sendBox=0;sendBox<level->num_my_boxes;sendBox++){
    int di,dj,dk;
    for(dk=-1;dk<=1;dk++){
    for(dj=-1;dj<=1;dj++){
    for(di=-1;di<=1;di++){
      int dir = 13+di+3*dj+9*dk;if(CommunicateThisDir[dir]){
      int       myBoxID = level->my_boxes[sendBox].global_box_id;
      int       myBox_i = level->my_boxes[sendBox].low.i / level->box_dim;
      int       myBox_j = level->my_boxes[sendBox].low.j / level->box_dim;
      int       myBox_k = level->my_boxes[sendBox].low.k / level->box_dim;
      int neighborBoxID = -1;
      if(level->boundary_condition.type == BC_PERIODIC){
        int neighborBox_i = (  myBox_i + di + level->boxes_in.i) % level->boxes_in.i;
        int neighborBox_j = (  myBox_j + dj + level->boxes_in.j) % level->boxes_in.j;
        int neighborBox_k = (  myBox_k + dk + level->boxes_in.k) % level->boxes_in.k;
            neighborBoxID =  neighborBox_i + neighborBox_j*level->boxes_in.i + neighborBox_k*level->boxes_in.i*level->boxes_in.j;
      }else{
        int neighborBox_i = (  myBox_i + di );
        int neighborBox_j = (  myBox_j + dj );
        int neighborBox_k = (  myBox_k + dk );
        if( (neighborBox_i>=0) && (neighborBox_i<level->boxes_in.i) && 
            (neighborBox_j>=0) && (neighborBox_j<level->boxes_in.j) && 
            (neighborBox_k>=0) && (neighborBox_k<level->boxes_in.k) ){  // i.e. the neighbor is a valid box
            neighborBoxID =  neighborBox_i + neighborBox_j*level->boxes_in.i + neighborBox_k*level->boxes_in.i*level->boxes_in.j;
        }
      }
      if(neighborBoxID>=0){
      if( level->rank_of_box[neighborBoxID] != -1 ){
        ghostsToSend[numGhosts].sendRank  = level->my_rank;
        ghostsToSend[numGhosts].sendBoxID = myBoxID;
        ghostsToSend[numGhosts].sendBox   = sendBox;
        ghostsToSend[numGhosts].sendDir   = dir;
        ghostsToSend[numGhosts].recvRank  = level->rank_of_box[neighborBoxID];
        ghostsToSend[numGhosts].recvBoxID = neighborBoxID;
        ghostsToSend[numGhosts].recvBox   = -1;
        if( level->rank_of_box[neighborBoxID] != level->my_rank ){
          sendRanks[numGhostsRemote++] = level->rank_of_box[neighborBoxID];
        }else{
          int recvBox=0;while(level->my_boxes[recvBox].global_box_id!=neighborBoxID)recvBox++; // search my list of boxes for the appropriate recvBox index
          ghostsToSend[numGhosts].recvBox   = recvBox;
        }
        numGhosts++;
      }}
    }}}}
  }
  // sort boxes by sendRank(==my rank) then by sendBoxID... ensures the sends and receive buffers are always sorted by sendBoxID...
  qsort(ghostsToSend,numGhosts      ,sizeof(GZ_type),qsortGZ );
  // sort the lists of neighboring ranks and remove duplicates...
  qsort(sendRanks   ,numGhostsRemote,sizeof(    int),qsortInt);
  int numSendRanks=0;_rank=-1;for(ghost=0;ghost<numGhostsRemote;ghost++)if(sendRanks[ghost] != _rank){_rank=sendRanks[ghost];sendRanks[numSendRanks++]=sendRanks[ghost];}


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // in a two-stage process, traverse the list of ghosts and allocate the pack/local lists as well as the MPI buffers, and then populate the pack/local lists
  level->exchange_ghosts[shape].num_sends     =                  numSendRanks;
  level->exchange_ghosts[shape].send_ranks    =     (int*)malloc(numSendRanks*sizeof(int));
  level->exchange_ghosts[shape].send_sizes    =     (int*)malloc(numSendRanks*sizeof(int));
  level->exchange_ghosts[shape].send_buffers  = (double**)malloc(numSendRanks*sizeof(double*));
  if(numSendRanks>0){
  if(level->exchange_ghosts[shape].send_ranks  ==NULL){fprintf(stderr,"malloc failed - exchange_ghosts[%d].send_ranks\n",shape);exit(0);}
  if(level->exchange_ghosts[shape].send_sizes  ==NULL){fprintf(stderr,"malloc failed - exchange_ghosts[%d].send_sizes\n",shape);exit(0);}
  if(level->exchange_ghosts[shape].send_buffers==NULL){fprintf(stderr,"malloc failed - exchange_ghosts[%d].send_buffers\n",shape);exit(0);}
  }
  level->exchange_ghosts[shape].blocks[0] = NULL;
  level->exchange_ghosts[shape].blocks[1] = NULL;
  level->exchange_ghosts[shape].num_blocks[0] = 0;
  level->exchange_ghosts[shape].num_blocks[1] = 0;
  level->exchange_ghosts[shape].allocated_blocks[0] = 0;
  level->exchange_ghosts[shape].allocated_blocks[1] = 0;
  for(stage=0;stage<=1;stage++){
    // stage=0... traverse the list and calculate the buffer sizes
    // stage=1... allocate MPI send buffers, traverse the list, and populate the unpack/local lists...
    int neighbor;
    for(neighbor=0;neighbor<numSendRanks;neighbor++){
      if(stage==1){
             level->exchange_ghosts[shape].send_buffers[neighbor] = (double*)malloc(level->exchange_ghosts[shape].send_sizes[neighbor]*sizeof(double));
          if(level->exchange_ghosts[shape].send_sizes[neighbor]>0)
          if(level->exchange_ghosts[shape].send_buffers[neighbor]==NULL){fprintf(stderr,"malloc failed - exchange_ghosts[%d].send_buffers[neighbor]\n",shape);exit(0);}
      memset(level->exchange_ghosts[shape].send_buffers[neighbor],                0,level->exchange_ghosts[shape].send_sizes[neighbor]*sizeof(double));
      }
      level->exchange_ghosts[shape].send_ranks[neighbor]=sendRanks[neighbor];
      level->exchange_ghosts[shape].send_sizes[neighbor]=0;
    }
    for(ghost=0;ghost<numGhosts;ghost++){
      int  dim_i=-1, dim_j=-1, dim_k=-1;
      int send_i=-1,send_j=-1,send_k=-1;
      int recv_i=-1,recv_j=-1,recv_k=-1;
  
      // decode ghostsToSend[ghost].sendDir (direction sent) into di/dj/dk 
      int di = ((ghostsToSend[ghost].sendDir % 3)  )-1;
      int dj = ((ghostsToSend[ghost].sendDir % 9)/3)-1;
      int dk = ((ghostsToSend[ghost].sendDir / 9)  )-1;
      switch(di){ // direction relative to sender
        case -1:send_i=0;                               dim_i=level->box_ghosts;recv_i=  level->box_dim;   break;
        case  0:send_i=0;                               dim_i=level->box_dim;   recv_i=0;                  break;
        case  1:send_i=level->box_dim-level->box_ghosts;dim_i=level->box_ghosts;recv_i=0-level->box_ghosts;break;
      }
      switch(dj){ // direction relative to sender
        case -1:send_j=0;                               dim_j=level->box_ghosts;recv_j=  level->box_dim;   break;
        case  0:send_j=0;                               dim_j=level->box_dim;   recv_j=0;                  break;
        case  1:send_j=level->box_dim-level->box_ghosts;dim_j=level->box_ghosts;recv_j=0-level->box_ghosts;break;
      }
      switch(dk){ // direction relative to sender
        case -1:send_k=0;                               dim_k=level->box_ghosts;recv_k=  level->box_dim;   break;
        case  0:send_k=0;                               dim_k=level->box_dim;   recv_k=0;                  break;
        case  1:send_k=level->box_dim-level->box_ghosts;dim_k=level->box_ghosts;recv_k=0-level->box_ghosts;break;
      }
 
      // determine if this ghost requires a pack or local exchange 
      int LocalExchange; // 0 = pack list, 1 = local exchange list
      if(ghostsToSend[ghost].recvRank != level->my_rank){
        LocalExchange=0; // pack
        neighbor=0;while(level->exchange_ghosts[shape].send_ranks[neighbor] != ghostsToSend[ghost].recvRank)neighbor++;
      }else{
        LocalExchange=1; // local
        neighbor=-1;
      }
   
      if(stage==1){ 
      if(LocalExchange) // append to the local exchange list...
      append_block_to_list(&(level->exchange_ghosts[shape].blocks[1]),&(level->exchange_ghosts[shape].allocated_blocks[1]),&(level->exchange_ghosts[shape].num_blocks[1]),
        /* dim.i         = */ dim_i,
        /* dim.j         = */ dim_j,
        /* dim.k         = */ dim_k,
        /* read.box      = */ ghostsToSend[ghost].sendBox,
        /* read.ptr      = */ NULL,
        /* read.i        = */ send_i,
        /* read.j        = */ send_j,
        /* read.k        = */ send_k,
        /* read.jStride  = */ level->my_boxes[ghostsToSend[ghost].sendBox].jStride,
        /* read.kStride  = */ level->my_boxes[ghostsToSend[ghost].sendBox].kStride,
        /* read.scale    = */ 1,
        /* write.box     = */ ghostsToSend[ghost].recvBox,
        /* write.ptr     = */ NULL,
        /* write.i       = */ recv_i,
        /* write.j       = */ recv_j,
        /* write.k       = */ recv_k,
        /* write.jStride = */ level->my_boxes[ghostsToSend[ghost].recvBox].jStride,
        /* write.kStride = */ level->my_boxes[ghostsToSend[ghost].recvBox].kStride,
        /* write.scale   = */ 1,
        /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
        /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
        /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
        /* subtype       = */ 0  
      );
      else // append to the MPI pack list...
      append_block_to_list(&(level->exchange_ghosts[shape].blocks[0]),&(level->exchange_ghosts[shape].allocated_blocks[0]),&(level->exchange_ghosts[shape].num_blocks[0]),
        /* dim.i         = */ dim_i,
        /* dim.j         = */ dim_j,
        /* dim.k         = */ dim_k,
        /* read.box      = */ ghostsToSend[ghost].sendBox,
        /* read.ptr      = */ NULL,
        /* read.i        = */ send_i,
        /* read.j        = */ send_j,
        /* read.k        = */ send_k,
        /* read.jStride  = */ level->my_boxes[ghostsToSend[ghost].sendBox].jStride,
        /* read.kStride  = */ level->my_boxes[ghostsToSend[ghost].sendBox].kStride,
        /* read.scale    = */ 1,
        /* write.box     = */ -1,
        /* write.ptr     = */ level->exchange_ghosts[shape].send_buffers[neighbor], // NOTE, 1. count _sizes, 2. allocate _buffers, 3. populate blocks
        /* write.i       = */ level->exchange_ghosts[shape].send_sizes[neighbor], // current offset in the MPI send buffer
        /* write.j       = */ 0,
        /* write.k       = */ 0,
        /* write.jStride = */ dim_i,       // contiguous block
        /* write.kStride = */ dim_i*dim_j, // contiguous block
        /* write.scale   = */ 1,
        /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
        /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
        /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
        /* subtype       = */ 0  
      );}
      if(neighbor>=0)level->exchange_ghosts[shape].send_sizes[neighbor]+=dim_i*dim_j*dim_k;
    } // ghost for-loop
  } // stage for-loop


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // free temporary storage...
  free(ghostsToSend);
  free(sendRanks);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // traverse my list of boxes and create a lists of neighboring boxes and neighboring ranks
  GZ_type *ghostsToRecv = (GZ_type*)malloc(26*level->num_my_boxes*sizeof(GZ_type)); // There are at most 26 neighbors per box.
         int *recvRanks = (    int*)malloc(26*level->num_my_boxes*sizeof(    int)); // There are at most 26 neighbors per box.
  if(level->num_my_boxes>0){
  if(ghostsToRecv == NULL){fprintf(stderr,"malloc failed - build_exchange_ghosts/ghostsToRecv\n");exit(0);}
  if(recvRanks    == NULL){fprintf(stderr,"malloc failed - build_exchange_ghosts/recvRanks   \n");exit(0);}
  }
  numGhosts       = 0;
  numGhostsRemote = 0;
  for(recvBox=0;recvBox<level->num_my_boxes;recvBox++){
    int di,dj,dk;
    for(dk=-1;dk<=1;dk++){
    for(dj=-1;dj<=1;dj++){
    for(di=-1;di<=1;di++){
      int dir = 13+di+3*dj+9*dk;if(CommunicateThisDir[dir]){
      int       myBoxID = level->my_boxes[recvBox].global_box_id;
      int       myBox_i = level->my_boxes[recvBox].low.i / level->box_dim;
      int       myBox_j = level->my_boxes[recvBox].low.j / level->box_dim;
      int       myBox_k = level->my_boxes[recvBox].low.k / level->box_dim;
      int neighborBoxID = -1;
      if(level->boundary_condition.type == BC_PERIODIC){
        int neighborBox_i = (  myBox_i + di + level->boxes_in.i) % level->boxes_in.i;
        int neighborBox_j = (  myBox_j + dj + level->boxes_in.j) % level->boxes_in.j;
        int neighborBox_k = (  myBox_k + dk + level->boxes_in.k) % level->boxes_in.k;
            neighborBoxID =  neighborBox_i + neighborBox_j*level->boxes_in.i + neighborBox_k*level->boxes_in.i*level->boxes_in.j;
      }else{
        int neighborBox_i = (  myBox_i + di );
        int neighborBox_j = (  myBox_j + dj );
        int neighborBox_k = (  myBox_k + dk );
        if( (neighborBox_i>=0) && (neighborBox_i<level->boxes_in.i) && 
            (neighborBox_j>=0) && (neighborBox_j<level->boxes_in.j) && 
            (neighborBox_k>=0) && (neighborBox_k<level->boxes_in.k) ){  // i.e. the neighbor is a valid box
            neighborBoxID =  neighborBox_i + neighborBox_j*level->boxes_in.i + neighborBox_k*level->boxes_in.i*level->boxes_in.j;
        }
      }
      if(neighborBoxID>=0){
      if( (level->rank_of_box[neighborBoxID] != -1) && (level->rank_of_box[neighborBoxID] != level->my_rank)  ){
        ghostsToRecv[numGhosts].sendRank  = level->rank_of_box[neighborBoxID];
        ghostsToRecv[numGhosts].sendBoxID = neighborBoxID;
        ghostsToRecv[numGhosts].sendBox   = -1;
        ghostsToRecv[numGhosts].sendDir   = 26-dir;
        ghostsToRecv[numGhosts].recvRank  = level->my_rank;
        ghostsToRecv[numGhosts].recvBoxID = myBoxID;
        ghostsToRecv[numGhosts].recvBox   = recvBox;
                     numGhosts++;
        recvRanks[numGhostsRemote++] = level->rank_of_box[neighborBoxID];
      }}
    }}}}
  }
  // sort boxes by sendRank then by sendBoxID... ensures the recvs and receive buffers are always sorted by sendBoxID...
  qsort(ghostsToRecv,numGhosts      ,sizeof(GZ_type),qsortGZ );
  // sort the lists of neighboring ranks and remove duplicates...
  qsort(recvRanks   ,numGhostsRemote,sizeof(    int),qsortInt);
  int numRecvRanks=0;_rank=-1;for(ghost=0;ghost<numGhostsRemote;ghost++)if(recvRanks[ghost] != _rank){_rank=recvRanks[ghost];recvRanks[numRecvRanks++]=recvRanks[ghost];}


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // in a two-stage process, traverse the list of ghosts and allocate the unpack lists as well as the MPI buffers, and then populate the unpack list
  level->exchange_ghosts[shape].num_recvs     =                  numRecvRanks;
  level->exchange_ghosts[shape].recv_ranks    =     (int*)malloc(numRecvRanks*sizeof(int));
  level->exchange_ghosts[shape].recv_sizes    =     (int*)malloc(numRecvRanks*sizeof(int));
  level->exchange_ghosts[shape].recv_buffers  = (double**)malloc(numRecvRanks*sizeof(double*));
  if(numRecvRanks>0){
  if(level->exchange_ghosts[shape].recv_ranks  ==NULL){fprintf(stderr,"malloc failed - exchange_ghosts[%d].recv_ranks\n",shape);exit(0);}
  if(level->exchange_ghosts[shape].recv_sizes  ==NULL){fprintf(stderr,"malloc failed - exchange_ghosts[%d].recv_sizes\n",shape);exit(0);}
  if(level->exchange_ghosts[shape].recv_buffers==NULL){fprintf(stderr,"malloc failed - exchange_ghosts[%d].recv_buffers\n",shape);exit(0);}
  }
  level->exchange_ghosts[shape].blocks[2] = NULL;
  level->exchange_ghosts[shape].num_blocks[2] = 0;
  level->exchange_ghosts[shape].allocated_blocks[2] = 0;
  for(stage=0;stage<=1;stage++){
    // stage=0... traverse the list and calculate the buffer sizes
    // stage=1... allocate MPI recv buffers, traverse the list, and populate the unpack/local lists...
    int neighbor;
    for(neighbor=0;neighbor<numRecvRanks;neighbor++){
      if(stage==1){
             level->exchange_ghosts[shape].recv_buffers[neighbor] = (double*)malloc(level->exchange_ghosts[shape].recv_sizes[neighbor]*sizeof(double));
          if(level->exchange_ghosts[shape].recv_sizes[neighbor]>0)
          if(level->exchange_ghosts[shape].recv_buffers[neighbor]==NULL){fprintf(stderr,"malloc failed - exchange_ghosts[%d].recv_buffers[neighbor]\n",shape);exit(0);}
      memset(level->exchange_ghosts[shape].recv_buffers[neighbor],                0,level->exchange_ghosts[shape].recv_sizes[neighbor]*sizeof(double));
      }
      level->exchange_ghosts[shape].recv_ranks[neighbor]=recvRanks[neighbor];
      level->exchange_ghosts[shape].recv_sizes[neighbor]=0;
    }
    for(ghost=0;ghost<numGhosts;ghost++){
      int  dim_i=-1, dim_j=-1, dim_k=-1;
    //int send_i=-1,send_j=-1,send_k=-1;
      int recv_i=-1,recv_j=-1,recv_k=-1;
  
      // decode ghostsToRecv[ghost].sendDir (direction sent) into di/dj/dk 
      int di = ((ghostsToRecv[ghost].sendDir % 3)  )-1;
      int dj = ((ghostsToRecv[ghost].sendDir % 9)/3)-1;
      int dk = ((ghostsToRecv[ghost].sendDir / 9)  )-1;
      switch(di){ // direction relative to sender
        case -1:dim_i=level->box_ghosts;recv_i=  level->box_dim;   break;
        case  0:dim_i=level->box_dim;   recv_i=0;                  break;
        case  1:dim_i=level->box_ghosts;recv_i=0-level->box_ghosts;break;
      }
      switch(dj){ // direction relative to sender
        case -1:dim_j=level->box_ghosts;recv_j=  level->box_dim;   break;
        case  0:dim_j=level->box_dim;   recv_j=0;                  break;
        case  1:dim_j=level->box_ghosts;recv_j=0-level->box_ghosts;break;
      }
      switch(dk){ // direction relative to sender
        case -1:dim_k=level->box_ghosts;recv_k=  level->box_dim;   break;
        case  0:dim_k=level->box_dim;   recv_k=0;                  break;
        case  1:dim_k=level->box_ghosts;recv_k=0-level->box_ghosts;break;
      }
 
      // determine if this ghost requires a pack or local exchange 
      neighbor=0;while(level->exchange_ghosts[shape].recv_ranks[neighbor] != ghostsToRecv[ghost].sendRank)neighbor++;
      if(stage==1)append_block_to_list(&(level->exchange_ghosts[shape].blocks[2]),&(level->exchange_ghosts[shape].allocated_blocks[2]),&(level->exchange_ghosts[shape].num_blocks[2]),
      /*dim.i         = */ dim_i,
      /*dim.j         = */ dim_j,
      /*dim.k         = */ dim_k,
      /*read.box      = */ -1,
      /*read.ptr      = */ level->exchange_ghosts[shape].recv_buffers[neighbor], // NOTE, 1. count _sizes, 2. allocate _buffers, 3. populate blocks
      /*read.i        = */ level->exchange_ghosts[shape].recv_sizes[neighbor], // current offset in the MPI recv buffer
      /*read.j        = */ 0,
      /*read.k        = */ 0,
      /*read.jStride  = */ dim_i,       // contiguous block
      /*read.kStride  = */ dim_i*dim_j, // contiguous block
      /*read.scale    = */ 1,
      /*write.box     = */ ghostsToRecv[ghost].recvBox,
      /*write.ptr     = */ NULL,
      /*write.i       = */ recv_i,
      /*write.j       = */ recv_j,
      /*write.k       = */ recv_k,
      /*write.jStride = */ level->my_boxes[ghostsToRecv[ghost].recvBox].jStride,
      /*write.kStride = */ level->my_boxes[ghostsToRecv[ghost].recvBox].kStride,
      /*write.scale   = */ 1,
      /* blockcopy_i  = */ BLOCKCOPY_TILE_I, // default
      /* blockcopy_j  = */ BLOCKCOPY_TILE_J, // default
      /* blockcopy_k  = */ BLOCKCOPY_TILE_K, // default
      /* subtype      = */ 0  
      );
      if(neighbor>=0)level->exchange_ghosts[shape].recv_sizes[neighbor]+=dim_i*dim_j*dim_k;
    } // ghost for-loop
  } // stage for-loop


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // free temporary storage...
  free(ghostsToRecv);
  free(recvRanks);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // malloc MPI requests/status arrays
  #ifdef USE_MPI
  level->exchange_ghosts[shape].requests = (MPI_Request*)malloc((level->exchange_ghosts[shape].num_sends+level->exchange_ghosts[shape].num_recvs)*sizeof(MPI_Request));
  level->exchange_ghosts[shape].status   = (MPI_Status *)malloc((level->exchange_ghosts[shape].num_sends+level->exchange_ghosts[shape].num_recvs)*sizeof(MPI_Status ));
  if((level->exchange_ghosts[shape].num_sends+level->exchange_ghosts[shape].num_recvs)>0){
  if(level->exchange_ghosts[shape].requests==NULL){fprintf(stderr,"malloc failed - exchange_ghosts[%d].requests\n",shape);exit(0);}
  if(level->exchange_ghosts[shape].status  ==NULL){fprintf(stderr,"malloc failed - exchange_ghosts[%d].status\n",shape);exit(0);}
  }
  #endif

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #ifdef BLOCK_SPATIAL_SORT
  // sort all the resultant blocks by box,k,j,i (good locality)
  qsort(level->exchange_ghosts[shape].blocks[0],level->exchange_ghosts[shape].num_blocks[0],sizeof(blockCopy_type),qsortBlock);
  qsort(level->exchange_ghosts[shape].blocks[1],level->exchange_ghosts[shape].num_blocks[1],sizeof(blockCopy_type),qsortBlock);
  qsort(level->exchange_ghosts[shape].blocks[2],level->exchange_ghosts[shape].num_blocks[2],sizeof(blockCopy_type),qsortBlock);
  #endif
}



//---------------------------------------------------------------------------------------------------------------------------------------------------
// create the pointers in level_type to the contiguous vector FP data (useful for bulk copies to/from accelerators)
// create the pointers in each box to their respective segment of the level's vector FP data (useful for box-relative operators)
// if( (level->numVectors > 0) && (numVectors > level->numVectors) ) then allocate additional space for (numVectors-level->numVectors) and copy old leve->numVectors data
void create_vectors(level_type *level, int numVectors){
  if(numVectors <= level->numVectors)return; // already have enough space
  double          * old_vectors_base = level->vectors_base; // save a pointer to the originally allocated data for subsequent free()
  double               * old_vector0 = NULL;
  if(level->numVectors>0)old_vector0 = level->vectors[0];   // save a pointer to old FP data to copy


  // calculate the size of each box...
  level->box_jStride =                    (level->box_dim+2*level->box_ghosts);while(level->box_jStride % BOX_ALIGN_JSTRIDE)level->box_jStride++; // pencil
  level->box_kStride = level->box_jStride*(level->box_dim+2*level->box_ghosts);while(level->box_kStride % BOX_ALIGN_KSTRIDE)level->box_kStride++; // plane
  level->box_volume  = level->box_kStride*(level->box_dim+2*level->box_ghosts);while(level->box_volume  % BOX_ALIGN_VOLUME )level->box_volume++;  // volume


  #define VECTOR_MALLOC_BULK
  #ifdef  VECTOR_MALLOC_BULK
    // allocate one aligned, double-precision array and divide it among vectors...
    uint64_t malloc_size = (uint64_t)numVectors*level->num_my_boxes*level->box_volume*sizeof(double) + 4096;
    level->vectors_base = (double*)malloc(malloc_size);
    if((numVectors>0)&&(level->vectors_base==NULL)){fprintf(stderr,"malloc failed - level->vectors_base\n");exit(0);}
    double * tmpbuf = level->vectors_base;
    while( (uint64_t)(tmpbuf+level->box_ghosts*(1+level->box_jStride+level->box_kStride)) & 0xff ){tmpbuf++;} // align first *non-ghost* zone element of first component to a 256-Byte boundary
    uint64_t ofs;
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(ofs=0;ofs<(uint64_t)numVectors*level->num_my_boxes*level->box_volume;ofs++){tmpbuf[ofs]=0.0;} // Faster in MPI+OpenMP environments, but not NUMA-aware
    // if there is existing FP data... copy it, then free old data and pointer array
    if(level->numVectors>0){
      memcpy(tmpbuf,old_vector0,(uint64_t)level->numVectors*level->num_my_boxes*level->box_volume*sizeof(double)); // FIX... omp thread ???
      if(old_vectors_base)free(old_vectors_base); // free old data...
    }
    // allocate an array of pointers which point to the union of boxes for each vector
    // NOTE, this requires just one copyin per vector to an accelerator rather than requiring one copyin per box per vector
    if(level->numVectors>0)free(level->vectors); // free any previously allocated vector array
    level->vectors = (double **)malloc(numVectors*sizeof(double*));
    if((numVectors>0)&&(level->vectors==NULL)){fprintf(stderr,"malloc failed - level->vectors\n");exit(0);}
    uint64_t c;for(c=0;c<numVectors;c++){level->vectors[c] = tmpbuf + (uint64_t)c*level->num_my_boxes*level->box_volume;}
  #else
    // allocate vectors individually (simple, but may cause conflict misses)
    double ** old_vectors = level->vectors;
    level->vectors = (double **)malloc(numVectors*sizeof(double*));
    uint64_t c;
    for(c=                0;c<level->numVectors;c++){level->vectors[c] = old_vectors[c];}
    for(c=level->numVectors;c<       numVectors;c++){
      level->vectors[c] = (double*)malloc((uint64_t)level->num_my_boxes*level->box_volume*sizeof(double));
      uint64_t ofs;
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for(ofs=0;ofs<(uint64_t)level->num_my_boxes*level->box_volume;ofs++){level->vectors[c][ofs]=0.0;} // Faster in MPI+OpenMP environments, but not NUMA-aware
    }
    free(old_vectors);
  #endif


  // build the list of boxes...
  int box=0;
  int i,j,k;
  for(k=0;k<level->boxes_in.k;k++){
  for(j=0;j<level->boxes_in.j;j++){
  for(i=0;i<level->boxes_in.i;i++){
    int jStride = level->boxes_in.i;
    int kStride = level->boxes_in.i*level->boxes_in.j;
    int b=i + j*jStride + k*kStride;
    if(level->rank_of_box[b]==level->my_rank){
      if(level->numVectors>0)free(level->my_boxes[box].vectors); // free previously allocated vector array
      level->my_boxes[box].vectors = (double **)malloc(numVectors*sizeof(double*));
      if((numVectors>0)&&(level->my_boxes[box].vectors==NULL)){fprintf(stderr,"malloc failed - level->my_boxes[box].vectors\n");exit(0);}
      uint64_t c;for(c=0;c<numVectors;c++){level->my_boxes[box].vectors[c] = level->vectors[c] + (uint64_t)box*level->box_volume;}
      level->my_boxes[box].numVectors = numVectors;
      level->my_boxes[box].dim        = level->box_dim;
      level->my_boxes[box].ghosts     = level->box_ghosts;
      level->my_boxes[box].jStride    = level->box_jStride;
      level->my_boxes[box].kStride    = level->box_kStride;
      level->my_boxes[box].volume     = level->box_volume;
      level->my_boxes[box].low.i      = i*level->box_dim;
      level->my_boxes[box].low.j      = j*level->box_dim;
      level->my_boxes[box].low.k      = k*level->box_dim;
      level->my_boxes[box].global_box_id = b;
      box++;
  }}}}

  // level now has created/initialized vector FP data
  level->numVectors = numVectors;
}


//---------------------------------------------------------------------------------------------------------------------------------------------------
// create a level by populating the basic data structure, distribute boxes within the level among processes, allocate memory, and create any auxilliaries
// box_ghosts must be >= stencil_get_radius()
// numVectors represents an estimate of the number of vectors needed in this level.  Additional vectors can be added via subsequent calls to create_vectors()
void create_level(level_type *level, int boxes_in_i, int box_dim, int box_ghosts, int numVectors, int domain_boundary_condition, int my_rank, int num_ranks, const MPI_Comm comm){
  int box;
  int TotalBoxes = boxes_in_i*boxes_in_i*boxes_in_i;

  if(my_rank==0){
  //if(domain_boundary_condition==BC_DIRICHLET)fprintf(stdout,"\nattempting to create a %d^3 level (with Dirichlet BC) using a %d^3 grid of %d^3 boxes and %d tasks...\n",box_dim*boxes_in_i,boxes_in_i,box_dim,num_ranks);
  //if(domain_boundary_condition==BC_PERIODIC )fprintf(stdout,"\nattempting to create a %d^3 level (with Periodic BC) using a %d^3 grid of %d^3 boxes and %d tasks...\n", box_dim*boxes_in_i,boxes_in_i,box_dim,num_ranks);
                                               fprintf(stdout,"\nattempting to create a %d^3 level from %d x %d^3 boxes distributed among %d tasks...\n", box_dim*boxes_in_i,TotalBoxes,box_dim,num_ranks);
    if(domain_boundary_condition==BC_DIRICHLET)fprintf(stdout,"  boundary condition = BC_DIRICHLET\n");
    if(domain_boundary_condition==BC_PERIODIC )fprintf(stdout,"  boundary condition = BC_PERIODIC\n");
    
  }

  int omp_threads = 1;

  #ifdef _OPENMP
  #pragma omp parallel 
  {
    #pragma omp master
    {
      omp_threads = omp_get_num_threads();
    }
  }
  #endif

  if(box_ghosts < stencil_get_radius() ){
    if(my_rank==0)fprintf(stderr,"ghosts(%d) must be >= stencil_get_radius(%d)\n",box_ghosts,stencil_get_radius());
    exit(0);
  }

  level->box_dim        = box_dim;
  level->box_ghosts     = box_ghosts;
  level->numVectors     = 0; // no vectors have been allocated yet
  level->vectors_base   = NULL; // pointer returned by bulk malloc
  level->vectors        = NULL; // pointers to individual vectors
  level->boxes_in.i     = boxes_in_i;
  level->boxes_in.j     = boxes_in_i;
  level->boxes_in.k     = boxes_in_i;
  level->dim.i          = box_dim*level->boxes_in.i;
  level->dim.j          = box_dim*level->boxes_in.j;
  level->dim.k          = box_dim*level->boxes_in.k;
  level->active         = 1;
  level->my_rank        = my_rank;
  level->num_ranks      = num_ranks;
  level->boundary_condition.type = domain_boundary_condition;
  level->must_subtract_mean = -1;
  level->num_threads      = omp_threads;
  level->my_blocks        = NULL;
  level->num_my_blocks    = 0;
  level->allocated_blocks = 0;
  level->tag              = log2(level->dim.i);
  level->fluxes           = NULL;


  // allocate 3D array of integers to hold the MPI rank of the corresponding box and initialize to -1 (unassigned)
     level->rank_of_box = (int*)malloc(level->boxes_in.i*level->boxes_in.j*level->boxes_in.k*sizeof(int));
  if(level->rank_of_box==NULL){fprintf(stderr,"malloc of level->rank_of_box failed\n");exit(0);}
  for(box=0;box<level->boxes_in.i*level->boxes_in.j*level->boxes_in.k;box++){level->rank_of_box[box]=-1;}  // -1 denotes that there is no actual box assigned to this region


  // parallelize the level (i.e. assign a process rank to each box)...
  #ifdef DECOMPOSE_LEX
  // lexicographical ordering... good load balance, potentially high bisection bandwidth requirements, bad surface:volume ratio when #boxes/proc is large
  if(my_rank==0){fprintf(stdout,"  Decomposing level via lexicographical ordering... ");fflush(stdout);}
  decompose_level_lex(level->rank_of_box,level->boxes_in.i,level->boxes_in.j,level->boxes_in.k,num_ranks);
  #elif DECOMPOSE_BISECTION_SPECIAL
  // recursive partitioning by primes
  if(my_rank==0){fprintf(stdout,"  Decomposing level via partitioning by primes... ");fflush(stdout);}
  decompose_level_bisection_special(level->rank_of_box,level->boxes_in.i,level->boxes_in.i*level->boxes_in.j,0,0,0,level->boxes_in.i,level->boxes_in.j,level->boxes_in.k,0,num_ranks);
  #elif DECOMPOSE_BISECTION
  // recursive bisection
  if(my_rank==0){fprintf(stdout,"  Decomposing level via recursive bisection... ");fflush(stdout);}
  decompose_level_bisection(level->rank_of_box,level->boxes_in.i,level->boxes_in.i*level->boxes_in.j,0,0,0,level->boxes_in.i,level->boxes_in.j,level->boxes_in.k,num_ranks,0,level->boxes_in.i*level->boxes_in.j*level->boxes_in.k);
  #else//#elif DECOMPOSE_ZMORT
  if(my_rank==0){fprintf(stdout,"  Decomposing level via Z-mort ordering... ");fflush(stdout);}
  #if 0 // Z-Mort over a power of two bounding box skipping boxes outside the domain
  int idim_padded=1;while(idim_padded<level->boxes_in.i)idim_padded*=2;
  int jdim_padded=1;while(jdim_padded<level->boxes_in.j)jdim_padded*=2;
  int kdim_padded=1;while(kdim_padded<level->boxes_in.k)kdim_padded*=2;
  #else // Z-Mort over the valid domain wtih odd-sized base cases (i.e. zmort on 3x3)
  int idim_padded=level->boxes_in.i;
  int jdim_padded=level->boxes_in.j;
  int kdim_padded=level->boxes_in.k;
  #endif
  decompose_level_zmort(level->rank_of_box,level->boxes_in.i,level->boxes_in.j,level->boxes_in.k,0,0,0,idim_padded,jdim_padded,kdim_padded,num_ranks,0,level->boxes_in.i*level->boxes_in.j*level->boxes_in.k);
  #endif
  if(my_rank==0){fprintf(stdout,"done\n");fflush(stdout);}
//print_decomposition(level);// for debug purposes only


  // calculate how many boxes I own...
  level->num_my_boxes=0;
  for(box=0;box<level->boxes_in.i*level->boxes_in.j*level->boxes_in.k;box++){if(level->rank_of_box[box]==level->my_rank)level->num_my_boxes++;} 
  level->my_boxes = (box_type*)malloc(level->num_my_boxes*sizeof(box_type));
  if((level->num_my_boxes>0)&&(level->my_boxes==NULL)){fprintf(stderr,"malloc failed - create_level/level->my_boxes\n");exit(0);}


  // allocate flattened vector FP data and create pointers...
  if(my_rank==0){fprintf(stdout,"  Allocating vectors... ");fflush(stdout);}
  create_vectors(level,numVectors);
  if(my_rank==0){fprintf(stdout,"done\n");fflush(stdout);}


  // Build and auxilarlly data structure that flattens boxes into blocks...
  for(box=0;box<level->num_my_boxes;box++){
    int blockcopy_i = BLOCKCOPY_TILE_I;
    int blockcopy_j = BLOCKCOPY_TILE_J;
    int blockcopy_k = BLOCKCOPY_TILE_K;

    append_block_to_list(&(level->my_blocks),&(level->allocated_blocks),&(level->num_my_blocks),
      /* dim.i         = */ level->my_boxes[box].dim,
      /* dim.j         = */ level->my_boxes[box].dim,
      /* dim.k         = */ level->my_boxes[box].dim,
      /* read.box      = */ box,
      /* read.ptr      = */ NULL,
      /* read.i        = */ 0,
      /* read.j        = */ 0,
      /* read.k        = */ 0,
      /* read.jStride  = */ level->my_boxes[box].jStride,
      /* read.kStride  = */ level->my_boxes[box].kStride,
      /* read.scale    = */ 1,
      /* write.box     = */ box,
      /* write.ptr     = */ NULL,
      /* write.i       = */ 0,
      /* write.j       = */ 0,
      /* write.k       = */ 0,
      /* write.jStride = */ level->my_boxes[box].jStride,
      /* write.kStride = */ level->my_boxes[box].kStride,
      /* write.scale   = */ 1,
      /* blockcopy_i   = */ blockcopy_i,
      /* blockcopy_j   = */ blockcopy_j,
      /* blockcopy_k   = */ blockcopy_k,
      /* subtype       = */ 0  
    );
  }


  // build an assist structure for Gauss Seidel Red Black that would facilitate unrolling and SIMDization...
  level->RedBlack_base = NULL;
  level->RedBlack_FP   = NULL;
  if(level->num_my_boxes){
    int i,j;
    int kStride = level->my_boxes[0].kStride;
    int jStride = level->my_boxes[0].jStride;
    level->RedBlack_base = (double*)malloc(2*kStride*sizeof(double)+256); // used for free()
    level->RedBlack_FP   = level->RedBlack_base; // aligned version
    // align first *non-ghost* zone element to a 64-Byte boundary...
    while( (uint64_t)(level->RedBlack_FP + level->box_ghosts*(1+level->box_jStride)) & 0x3f ){level->RedBlack_FP++;}
    // initialize RedBlack array...
    for(j=0-level->box_ghosts;j<level->box_dim+level->box_ghosts;j++){
    for(i=0-level->box_ghosts;i<level->box_dim+level->box_ghosts;i++){
      int ij = (i+level->box_ghosts) + (j+level->box_ghosts)*jStride;
      if((i^j^1)&0x1){
        level->RedBlack_FP[ij        ]=1.0;
        level->RedBlack_FP[ij+kStride]=0.0;
      }else{
        level->RedBlack_FP[ij        ]=0.0;
        level->RedBlack_FP[ij+kStride]=1.0;
      }
      // Never update ghost zones
      //if( (i<0) || (i>=level->box_dim) || (j<0) || (j>=level->box_dim) ){
      //  level->RedBlack_FP[ij        ]=0.0;
      //  level->RedBlack_FP[ij+kStride]=0.0;
      //}
    }}
  }


  int shape;
  // create mini program for each stencil shape to perform a ghost zone exchange...
  for(shape=0;shape<STENCIL_MAX_SHAPES;shape++)build_exchange_ghosts(    level,shape);
  // create mini program for each stencil shape to perform a boundary condition...
  for(shape=0;shape<STENCIL_MAX_SHAPES;shape++)build_boundary_conditions(level,shape);


  // duplicate MPI_COMM_WORLD to be the communicator for each level
  #ifdef USE_MPI
  level->MPI_COMM_ALLREDUCE = comm;
  /*
  if(my_rank==0){fprintf(stdout,"  Duplicating MPI_COMM_WORLD... ");fflush(stdout);}
  double time_start = MPI_Wtime();
  MPI_Comm_dup(comm,&level->MPI_COMM_ALLREDUCE);
  double time_end = MPI_Wtime();
  double time_in_comm_dup = 0;
  double time_in_comm_dup_send = time_end-time_start;
  MPI_Allreduce(&time_in_comm_dup_send,&time_in_comm_dup,1,MPI_DOUBLE,MPI_MAX,level->MPI_COMM_ALLREDUCE);
  if(my_rank==0){fprintf(stdout,"done (%0.6f seconds)\n",time_in_comm_dup);fflush(stdout);}
  */
  #endif

  // report on potential load imbalance
  int BoxesPerProcess = level->num_my_boxes;
  #ifdef USE_MPI
  int BoxesPerProcessSend = level->num_my_boxes;
  MPI_Allreduce(&BoxesPerProcessSend,&BoxesPerProcess,1,MPI_INT,MPI_MAX,level->MPI_COMM_ALLREDUCE);
  #endif
  if(my_rank==0){fprintf(stdout,"  Calculating boxes per process... target=%0.3f, max=%d\n",(double)TotalBoxes/(double)num_ranks,BoxesPerProcess);}
}



//---------------------------------------------------------------------------------------------------------------------------------------------------
// zeros are the timers within this level
// useful if one wishes to separate setup(build) timing from solve timing
void reset_level_timers(level_type *level){
  // cycle counters information...
  level->timers.smooth                  = 0;
  level->timers.apply_op                = 0;
  level->timers.residual                = 0;
  level->timers.blas1                   = 0;
  level->timers.blas3                   = 0;
  level->timers.boundary_conditions     = 0;
  level->timers.restriction_total       = 0;
  level->timers.restriction_pack        = 0;
  level->timers.restriction_local       = 0;
  level->timers.restriction_unpack      = 0;
  level->timers.restriction_recv        = 0;
  level->timers.restriction_send        = 0;
  level->timers.restriction_wait        = 0;
  level->timers.interpolation_total     = 0;
  level->timers.interpolation_pack      = 0;
  level->timers.interpolation_local     = 0;
  level->timers.interpolation_unpack    = 0;
  level->timers.interpolation_recv      = 0;
  level->timers.interpolation_send      = 0;
  level->timers.interpolation_wait      = 0;
  level->timers.ghostZone_total         = 0;
  level->timers.ghostZone_pack          = 0;
  level->timers.ghostZone_local         = 0;
  level->timers.ghostZone_unpack        = 0;
  level->timers.ghostZone_recv          = 0;
  level->timers.ghostZone_send          = 0;
  level->timers.ghostZone_wait          = 0;
  level->timers.collectives             = 0;
  level->timers.Total                   = 0;
  // solver events information...
  level->Krylov_iterations              = 0;
  level->CAKrylov_formations_of_G       = 0;
  level->vcycles_from_this_level        = 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------
// free all memory allocated by this level
// n.b. in some cases a malloc was used as the basis for an array of pointers.  As such free(x[0])
void destroy_level(level_type *level){
  int i,j;
  if(level->my_rank==0){fprintf(stdout,"attempting to free the %5d^3 level... ",level->dim.i);fflush(stdout);}

  // box ...
  for(i=0;i<level->num_my_boxes;i++)if(level->my_boxes[i].vectors)free(level->my_boxes[i].vectors);

  // misc ...
  if(level->rank_of_box )free(level->rank_of_box);
  if(level->my_boxes    )free(level->my_boxes);
  if(level->my_blocks   )free(level->my_blocks);
  if(level->RedBlack_base)free(level->RedBlack_base);

  // FP vector data...
  #ifdef VECTOR_MALLOC_BULK
  if(level->vectors_base)free(level->vectors_base);
  if(level->vectors     )free(level->vectors);
  #else
  for(i=0;i<level->numVectors;i++)if(level->vectors[i])free(level->vectors[i]);
  if(level->vectors     )free(level->vectors);
  #endif

  // boundary condition mini program...
  for(i=0;i<STENCIL_MAX_SHAPES;i++){
    if(level->boundary_condition.blocks[i])free(level->boundary_condition.blocks[i]);
  }

  // ghost zone exchange mini programs...
  for(i=0;i<STENCIL_MAX_SHAPES;i++){
    if(level->exchange_ghosts[i].num_recvs>0){
    for(j=0;j<level->exchange_ghosts[i].num_recvs;j++)if(level->exchange_ghosts[i].recv_buffers[j])free(level->exchange_ghosts[i].recv_buffers[j]);
    if(level->exchange_ghosts[i].recv_buffers)free(level->exchange_ghosts[i].recv_buffers);
    if(level->exchange_ghosts[i].recv_ranks  )free(level->exchange_ghosts[i].recv_ranks  );
    if(level->exchange_ghosts[i].recv_sizes  )free(level->exchange_ghosts[i].recv_sizes  );
    }
    if(level->exchange_ghosts[i].num_sends>0){
    for(j=0;j<level->exchange_ghosts[i].num_sends;j++)if(level->exchange_ghosts[i].send_buffers[j])free(level->exchange_ghosts[i].send_buffers[j]);
    if(level->exchange_ghosts[i].send_buffers)free(level->exchange_ghosts[i].send_buffers);
    if(level->exchange_ghosts[i].send_ranks  )free(level->exchange_ghosts[i].send_ranks  );
    if(level->exchange_ghosts[i].send_sizes  )free(level->exchange_ghosts[i].send_sizes  );
    }
    if(level->exchange_ghosts[i].blocks[0]   )free(level->exchange_ghosts[i].blocks[0]   );
    if(level->exchange_ghosts[i].blocks[1]   )free(level->exchange_ghosts[i].blocks[1]   );
    if(level->exchange_ghosts[i].blocks[2]   )free(level->exchange_ghosts[i].blocks[2]   );
    #ifdef USE_MPI
    if(level->exchange_ghosts[i].requests    )free(level->exchange_ghosts[i].requests    );
    if(level->exchange_ghosts[i].status      )free(level->exchange_ghosts[i].status      );
    #endif
  }

  /*
  // MPI subcommunicator
  #ifdef USE_MPI
  #ifdef USE_SUBCOMM
  MPI_Comm_free(&level->MPI_COMM_ALLREDUCE);
  #endif
  #endif
  */

  if(level->my_rank==0){fprintf(stdout,"done\n");}
}
