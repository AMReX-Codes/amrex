//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
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
#include "solvers.h"
#include "mg.h"
//------------------------------------------------------------------------------------------------------------------------------
// structs/routines used to construct the restriction and prolognation lists and ensure a convention on how data is ordered within an MPI buffer
typedef struct {
  int sendRank;
  int sendBoxID;
  int sendBox;
  int recvRank;
  int recvBoxID;
  int recvBox;
  int i,j,k; // offsets used to index into the coarse box
} RP_type;


int qsortRP(const void *a, const void*b){
  RP_type *rpa = (RP_type*)a;
  RP_type *rpb = (RP_type*)b;
  // sort first by sendRank
  if(rpa->sendRank < rpb->sendRank)return(-1);
  if(rpa->sendRank > rpb->sendRank)return( 1);
  // then by sendBoxID
  if(rpa->sendBoxID < rpb->sendBoxID)return(-1);
  if(rpa->sendBoxID > rpb->sendBoxID)return( 1);
  return(0);
}


//----------------------------------------------------------------------------------------------------------------------------------------------------
// print out average time per solve and then decompose by function and level
// note, in FMG, some levels are accessed more frequently.  This routine only prints time per solve in that level
void MGPrintTiming(mg_type *all_grids, int fromLevel){
  if(all_grids->my_rank!=0)return;
  int level,num_levels = all_grids->num_levels;
  #ifdef CALIBRATE_TIMER
  double _timeStart=getTime();sleep(1);double _timeEnd=getTime();
  double SecondsPerCycle = (double)1.0/(double)(_timeEnd-_timeStart);
  #else
  double SecondsPerCycle = 1.0;
  #endif
  double scale = SecondsPerCycle/(double)all_grids->MGSolves_performed; // prints average performance per MGSolve

  double time,total;
          printf("\n\n");
          printf("level                     ");for(level=fromLevel;level<(num_levels  );level++){printf("%12d ",level-fromLevel);}printf("\n");
          printf("level dimension           ");for(level=fromLevel;level<(num_levels  );level++){printf("%10d^3 ",all_grids->levels[level]->dim.i  );}printf("\n");
          printf("box dimension             ");for(level=fromLevel;level<(num_levels  );level++){printf("%10d^3 ",all_grids->levels[level]->box_dim);}printf("       total\n");
  total=0;printf("------------------        ");for(level=fromLevel;level<(num_levels+1);level++){printf("------------ ");}printf("\n");
  total=0;printf("smooth                    ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.smooth;               total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("residual                  ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.residual;             total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("applyOp                   ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.apply_op;             total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("BLAS1                     ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.blas1;                total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("BLAS3                     ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.blas3;                total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("Boundary Conditions       ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.boundary_conditions;  total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("Restriction               ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.restriction_total;    total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  local restriction       ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.restriction_local;    total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #ifdef USE_MPI
  total=0;printf("  pack MPI buffers        ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.restriction_pack;     total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  unpack MPI buffers      ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.restriction_unpack;   total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Isend               ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.restriction_send;     total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Irecv               ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.restriction_recv;     total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Waitall             ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.restriction_wait;     total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #endif
  total=0;printf("Interpolation             ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.interpolation_total;  total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  local interpolation     ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.interpolation_local;  total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #ifdef USE_MPI
  total=0;printf("  pack MPI buffers        ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.interpolation_pack;   total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  unpack MPI buffers      ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.interpolation_unpack; total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Isend               ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.interpolation_send;   total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Irecv               ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.interpolation_recv;   total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Waitall             ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.interpolation_wait;   total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #endif
  total=0;printf("Ghost Zone Exchange       ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.ghostZone_total;      total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  local exchange          ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.ghostZone_local;      total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #ifdef USE_MPI
  total=0;printf("  pack MPI buffers        ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.ghostZone_pack;       total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  unpack MPI buffers      ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.ghostZone_unpack;     total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Isend               ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.ghostZone_send;       total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Irecv               ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.ghostZone_recv;       total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Waitall             ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.ghostZone_wait;       total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #endif
  #ifdef USE_MPI
  total=0;printf("MPI_collectives           ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.collectives;          total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #endif
  total=0;printf("------------------        ");for(level=fromLevel;level<(num_levels+1);level++){printf("------------ ");}printf("\n");
  total=0;printf("Total by level            ");for(level=fromLevel;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->timers.Total;                total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);

  printf("\n");
  printf( "   Total time in MGBuild  %12.6f seconds\n",SecondsPerCycle*(double)all_grids->timers.MGBuild);
  printf( "   Total time in MGSolve  %12.6f seconds\n",scale*(double)all_grids->timers.MGSolve);
  printf( "      number of v-cycles  %12d\n"  ,all_grids->levels[fromLevel]->vcycles_from_this_level/all_grids->MGSolves_performed);
  printf( "Bottom solver iterations  %12d\n"  ,all_grids->levels[num_levels-1]->Krylov_iterations/all_grids->MGSolves_performed);
  #if defined(USE_CABICGSTAB) || defined(USE_CACG)
  printf( "     formations of G[][]  %12d\n"  ,all_grids->levels[num_levels-1]->CAKrylov_formations_of_G/all_grids->MGSolves_performed);
  #endif
  printf("\n\n");fflush(stdout);
}


//----------------------------------------------------------------------------------------------------------------------------------------------------
// zeros all timers within this MG hierarchy
void MGResetTimers(mg_type *all_grids){
  int level;
  for(level=0;level<all_grids->num_levels;level++)reset_level_timers(all_grids->levels[level]);
//all_grids->timers.MGBuild     = 0;
  all_grids->timers.MGSolve     = 0;
  all_grids->MGSolves_performed = 0;
}


//----------------------------------------------------------------------------------------------------------------------------------------------------
// build a list of operations and MPI buffers to affect distributed interpolation
// the three lists constitute
//   - buffer packing (i.e. interpolate a local box (or region of a box) and place the result in an MPI buffer)
//   - local operations (i.e. interpolate a local box (or region of a box) and place the result in another local box)
//   - buffer upacking (i.e. take interpolated data recieved from another process and use it to increment a local box)
void build_interpolation(mg_type *all_grids){
  int level;
  for(level=0;level<all_grids->num_levels;level++){

  // initialize to defaults...
  all_grids->levels[level]->interpolation.num_recvs           = 0;
  all_grids->levels[level]->interpolation.num_sends           = 0;
  all_grids->levels[level]->interpolation.recv_ranks          = NULL;
  all_grids->levels[level]->interpolation.send_ranks          = NULL;
  all_grids->levels[level]->interpolation.recv_sizes          = NULL;
  all_grids->levels[level]->interpolation.send_sizes          = NULL;
  all_grids->levels[level]->interpolation.recv_buffers        = NULL;
  all_grids->levels[level]->interpolation.send_buffers        = NULL;
  all_grids->levels[level]->interpolation.blocks[0]           = NULL;
  all_grids->levels[level]->interpolation.blocks[1]           = NULL;
  all_grids->levels[level]->interpolation.blocks[2]           = NULL;
  all_grids->levels[level]->interpolation.num_blocks[0]       = 0;
  all_grids->levels[level]->interpolation.num_blocks[1]       = 0;
  all_grids->levels[level]->interpolation.num_blocks[2]       = 0;
  all_grids->levels[level]->interpolation.allocated_blocks[0] = 0;
  all_grids->levels[level]->interpolation.allocated_blocks[1] = 0;
  all_grids->levels[level]->interpolation.allocated_blocks[2] = 0;
  #ifdef USE_MPI
  all_grids->levels[level]->interpolation.requests            = NULL;
  all_grids->levels[level]->interpolation.status              = NULL;
  #endif


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // construct pack, send(to level-1), and local...
  if( (level>0) && (all_grids->levels[level]->num_my_boxes>0) ){ // not top  *and*  I have boxes to send
    // construct a list of fine boxes to be coarsened and sent to me...
    int numFineBoxes = (all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i)*
                       (all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j)*
                       (all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k)*
                                                               all_grids->levels[level]->num_my_boxes;
        int *fineRanks = (    int*)malloc(numFineBoxes*sizeof(    int)); // high water mark (assumes every neighboring box is a different process)
    RP_type *fineBoxes = (RP_type*)malloc(numFineBoxes*sizeof(RP_type)); 
        numFineBoxes       = 0;
    int numFineBoxesLocal  = 0;
    int numFineBoxesRemote = 0;
    int coarseBox;
    for(coarseBox=0;coarseBox<all_grids->levels[level]->num_my_boxes;coarseBox++){
      int bi,bj,bk;
      int   coarseBoxID = all_grids->levels[level]->my_boxes[coarseBox].global_box_id;
      int   coarseBox_i = all_grids->levels[level]->my_boxes[coarseBox].low.i / all_grids->levels[level]->box_dim;
      int   coarseBox_j = all_grids->levels[level]->my_boxes[coarseBox].low.j / all_grids->levels[level]->box_dim;
      int   coarseBox_k = all_grids->levels[level]->my_boxes[coarseBox].low.k / all_grids->levels[level]->box_dim;
      for(bk=0;bk<all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k;bk++){
      for(bj=0;bj<all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j;bj++){
      for(bi=0;bi<all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i;bi++){
        int fineBox_i = (all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i)*coarseBox_i + bi;
        int fineBox_j = (all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j)*coarseBox_j + bj;
        int fineBox_k = (all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k)*coarseBox_k + bk;
        int fineBoxID =  fineBox_i + fineBox_j*all_grids->levels[level-1]->boxes_in.i + fineBox_k*all_grids->levels[level-1]->boxes_in.i*all_grids->levels[level-1]->boxes_in.j;
        int fineBox   = -1;int f;for(f=0;f<all_grids->levels[level-1]->num_my_boxes;f++)if( all_grids->levels[level-1]->my_boxes[f].global_box_id == fineBoxID )fineBox=f; // try and find the index of a fineBox global_box_id == fineBoxID
        fineBoxes[numFineBoxes].sendRank  = all_grids->levels[level  ]->rank_of_box[coarseBoxID];
        fineBoxes[numFineBoxes].sendBoxID = coarseBoxID;
        fineBoxes[numFineBoxes].sendBox   = coarseBox;
        fineBoxes[numFineBoxes].recvRank  = all_grids->levels[level-1]->rank_of_box[  fineBoxID];
        fineBoxes[numFineBoxes].recvBoxID = fineBoxID;
        fineBoxes[numFineBoxes].recvBox   = fineBox;
        fineBoxes[numFineBoxes].i         = bi*all_grids->levels[level-1]->box_dim/2;
        fineBoxes[numFineBoxes].j         = bj*all_grids->levels[level-1]->box_dim/2;
        fineBoxes[numFineBoxes].k         = bk*all_grids->levels[level-1]->box_dim/2;
                  numFineBoxes++;
        if(all_grids->levels[level-1]->rank_of_box[fineBoxID] != all_grids->levels[level]->my_rank){
          fineRanks[numFineBoxesRemote++] = all_grids->levels[level-1]->rank_of_box[fineBoxID];
        }else{numFineBoxesLocal++;}
      }}}
    } // my (coarse) boxes
    // sort boxes by sendRank(==my rank) then by sendBoxID... ensures the sends and receive buffers are always sorted by sendBoxID...
    qsort(fineBoxes,numFineBoxes      ,sizeof(RP_type),qsortRP );
    // sort the lists of neighboring ranks and remove duplicates...
    qsort(fineRanks,numFineBoxesRemote,sizeof(    int),qsortInt);
    int numFineRanks=0;
    int _rank=-1;int neighbor=0;
    for(neighbor=0;neighbor<numFineBoxesRemote;neighbor++)if(fineRanks[neighbor] != _rank){_rank=fineRanks[neighbor];fineRanks[numFineRanks++]=fineRanks[neighbor];}

    // allocate structures...
    all_grids->levels[level]->interpolation.num_sends     =                         numFineRanks;
    all_grids->levels[level]->interpolation.send_ranks    =            (int*)malloc(numFineRanks*sizeof(int));
    all_grids->levels[level]->interpolation.send_sizes    =            (int*)malloc(numFineRanks*sizeof(int));
    all_grids->levels[level]->interpolation.send_buffers  =        (double**)malloc(numFineRanks*sizeof(double*));
    if(numFineRanks>0){
    if(all_grids->levels[level]->interpolation.send_ranks  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.send_ranks\n",level);exit(0);}
    if(all_grids->levels[level]->interpolation.send_sizes  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.send_sizes\n",level);exit(0);}
    if(all_grids->levels[level]->interpolation.send_buffers==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.send_buffers\n",level);exit(0);}
    }

    int elementSize = all_grids->levels[level-1]->box_dim*all_grids->levels[level-1]->box_dim*all_grids->levels[level-1]->box_dim;
    double * all_send_buffers = (double*)malloc(numFineBoxesRemote*elementSize*sizeof(double));
          if(numFineBoxesRemote*elementSize>0)
          if(all_send_buffers==NULL){fprintf(stderr,"malloc failed - interpolation/all_send_buffers\n");exit(0);}
                      memset(all_send_buffers,0,numFineBoxesRemote*elementSize*sizeof(double)); // DO NOT DELETE... you must initialize to 0 to avoid getting something like 0.0*NaN and corrupting the solve
    //printf("level=%d, rank=%2d, send_buffers=%6d\n",level,all_grids->my_rank,numFineBoxesRemote*elementSize*sizeof(double));

    // for each neighbor, construct the pack list and allocate the MPI send buffer... 
    for(neighbor=0;neighbor<numFineRanks;neighbor++){
      int fineBox;
      int offset = 0;
      all_grids->levels[level]->interpolation.send_buffers[neighbor] = all_send_buffers;
      for(fineBox=0;fineBox<numFineBoxes;fineBox++)if(fineBoxes[fineBox].recvRank==fineRanks[neighbor]){
        // pack the MPI send buffer...
        append_block_to_list(&(all_grids->levels[level]->interpolation.blocks[0]),&(all_grids->levels[level]->interpolation.allocated_blocks[0]),&(all_grids->levels[level]->interpolation.num_blocks[0]),
          /* dim.i         = */ all_grids->levels[level-1]->box_dim/2,
          /* dim.j         = */ all_grids->levels[level-1]->box_dim/2,
          /* dim.k         = */ all_grids->levels[level-1]->box_dim/2,
          /* read.box      = */ fineBoxes[fineBox].sendBox,
          /* read.ptr      = */ NULL,
          /* read.i        = */ fineBoxes[fineBox].i,
          /* read.j        = */ fineBoxes[fineBox].j,
          /* read.k        = */ fineBoxes[fineBox].k,
          /* read.jStride  = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].sendBox].jStride,
          /* read.kStride  = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].sendBox].kStride,
          /* read.scale    = */ 1,
          /* write.box     = */ -1,
          /* write.ptr     = */ all_grids->levels[level]->interpolation.send_buffers[neighbor],
          /* write.i       = */ offset,
          /* write.j       = */ 0,
          /* write.k       = */ 0,
          /* write.jStride = */ all_grids->levels[level-1]->box_dim,
          /* write.kStride = */ all_grids->levels[level-1]->box_dim*all_grids->levels[level-1]->box_dim,
          /* write.scale   = */ 2,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
        offset+=elementSize;
      }
      all_grids->levels[level]->interpolation.send_ranks[neighbor] = fineRanks[neighbor];
      all_grids->levels[level]->interpolation.send_sizes[neighbor] = offset;
      all_send_buffers+=offset;
    } // neighbor
    {
      int fineBox;
      for(fineBox=0;fineBox<numFineBoxes;fineBox++)if(fineBoxes[fineBox].recvRank==all_grids->my_rank){
        // local interpolations...
        append_block_to_list(&(all_grids->levels[level]->interpolation.blocks[1]),&(all_grids->levels[level]->interpolation.allocated_blocks[1]),&(all_grids->levels[level]->interpolation.num_blocks[1]),
          /* dim.i         = */ all_grids->levels[level-1]->box_dim/2,
          /* dim.j         = */ all_grids->levels[level-1]->box_dim/2,
          /* dim.k         = */ all_grids->levels[level-1]->box_dim/2,
          /* read.box      = */ fineBoxes[fineBox].sendBox,
          /* read.ptr      = */ NULL,
          /* read.i        = */ fineBoxes[fineBox].i,
          /* read.j        = */ fineBoxes[fineBox].j,
          /* read.k        = */ fineBoxes[fineBox].k,
          /* read.jStride  = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].sendBox].jStride,
          /* read.kStride  = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].sendBox].kStride,
          /* read.scale    = */ 1,
          /* write.box     = */ fineBoxes[fineBox].recvBox,
          /* write.ptr     = */ NULL,
          /* write.i       = */ 0,
          /* write.j       = */ 0,
          /* write.k       = */ 0,
          /* write.jStride = */ all_grids->levels[level-1]->my_boxes[fineBoxes[fineBox].recvBox].jStride,
          /* write.kStride = */ all_grids->levels[level-1]->my_boxes[fineBoxes[fineBox].recvBox].kStride,
          /* write.scale   = */ 2,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
      }
    } // local to local interpolation

    // free temporary storage...
    free(fineBoxes);
    free(fineRanks);
  } // pack/send/local


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // construct recv(from level+1) and unpack...
  if( (level<all_grids->num_levels-1) && (all_grids->levels[level]->num_my_boxes>0) ){ // not bottom  *and*  I have boxes to receive

    // construct the list of coarsened boxes and neighboring ranks that will be interpolated and sent to me...
    int numCoarseBoxes = all_grids->levels[level]->num_my_boxes; // I may receive a block for each of my boxes
        int *coarseRanks = (    int*)malloc(numCoarseBoxes*sizeof(    int)); // high water mark (assumes every neighboring box is a different process)
    RP_type *coarseBoxes = (RP_type*)malloc(numCoarseBoxes*sizeof(RP_type)); 
        numCoarseBoxes       = 0;
    int fineBox;
    for(fineBox=0;fineBox<all_grids->levels[level]->num_my_boxes;fineBox++){
      int   fineBoxID = all_grids->levels[level]->my_boxes[fineBox].global_box_id;
      int   fineBox_i = all_grids->levels[level]->my_boxes[fineBox].low.i / all_grids->levels[level]->box_dim;
      int   fineBox_j = all_grids->levels[level]->my_boxes[fineBox].low.j / all_grids->levels[level]->box_dim;
      int   fineBox_k = all_grids->levels[level]->my_boxes[fineBox].low.k / all_grids->levels[level]->box_dim;
      int coarseBox_i = fineBox_i*all_grids->levels[level+1]->boxes_in.i/all_grids->levels[level]->boxes_in.i;
      int coarseBox_j = fineBox_j*all_grids->levels[level+1]->boxes_in.j/all_grids->levels[level]->boxes_in.j;
      int coarseBox_k = fineBox_k*all_grids->levels[level+1]->boxes_in.k/all_grids->levels[level]->boxes_in.k;
      int coarseBoxID =  coarseBox_i + coarseBox_j*all_grids->levels[level+1]->boxes_in.i + coarseBox_k*all_grids->levels[level+1]->boxes_in.i*all_grids->levels[level+1]->boxes_in.j;
      if(all_grids->levels[level]->my_rank != all_grids->levels[level+1]->rank_of_box[coarseBoxID]){
        coarseBoxes[numCoarseBoxes].sendRank  = all_grids->levels[level+1]->rank_of_box[coarseBoxID];
        coarseBoxes[numCoarseBoxes].sendBoxID = coarseBoxID;
        coarseBoxes[numCoarseBoxes].sendBox   = -1; 
        coarseBoxes[numCoarseBoxes].recvRank  = all_grids->levels[level  ]->rank_of_box[  fineBoxID];
        coarseBoxes[numCoarseBoxes].recvBoxID = fineBoxID;
        coarseBoxes[numCoarseBoxes].recvBox   = fineBox;
        coarseRanks[numCoarseBoxes] = all_grids->levels[level+1]->rank_of_box[coarseBoxID];
                    numCoarseBoxes++;
      }
    } // my (fine) boxes

    // sort boxes by sendRank(==my rank) then by sendBoxID... ensures the sends and receive buffers are always sorted by sendBoxID...
    qsort(coarseBoxes,numCoarseBoxes,sizeof(RP_type),qsortRP );
    // sort the lists of neighboring ranks and remove duplicates...
    qsort(coarseRanks,numCoarseBoxes,sizeof(    int),qsortInt);
    int numCoarseRanks=0;
    int _rank=-1;int neighbor=0;
    for(neighbor=0;neighbor<numCoarseBoxes;neighbor++)if(coarseRanks[neighbor] != _rank){_rank=coarseRanks[neighbor];coarseRanks[numCoarseRanks++]=coarseRanks[neighbor];}

    // allocate structures...
    all_grids->levels[level]->interpolation.num_recvs     =                         numCoarseRanks;
    all_grids->levels[level]->interpolation.recv_ranks    =            (int*)malloc(numCoarseRanks*sizeof(int));
    all_grids->levels[level]->interpolation.recv_sizes    =            (int*)malloc(numCoarseRanks*sizeof(int));
    all_grids->levels[level]->interpolation.recv_buffers  =        (double**)malloc(numCoarseRanks*sizeof(double*));
    if(numCoarseRanks>0){
    if(all_grids->levels[level]->interpolation.recv_ranks  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.recv_ranks\n",level);exit(0);}
    if(all_grids->levels[level]->interpolation.recv_sizes  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.recv_sizes\n",level);exit(0);}
    if(all_grids->levels[level]->interpolation.recv_buffers==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.recv_buffers\n",level);exit(0);}
    }

    int elementSize = all_grids->levels[level]->box_dim*all_grids->levels[level]->box_dim*all_grids->levels[level]->box_dim;
    double * all_recv_buffers = (double*)malloc(numCoarseBoxes*elementSize*sizeof(double)); 
          if(numCoarseBoxes*elementSize>0)
          if(all_recv_buffers==NULL){fprintf(stderr,"malloc failed - interpolation/all_recv_buffers\n");exit(0);}
                      memset(all_recv_buffers,0,numCoarseBoxes*elementSize*sizeof(double)); // DO NOT DELETE... you must initialize to 0 to avoid getting something like 0.0*NaN and corrupting the solve
    //printf("level=%d, rank=%2d, recv_buffers=%6d\n",level,all_grids->my_rank,numCoarseBoxes*elementSize*sizeof(double));

    // for each neighbor, construct the unpack list and allocate the MPI recv buffer... 
    for(neighbor=0;neighbor<numCoarseRanks;neighbor++){
      int coarseBox;
      int offset = 0;
      all_grids->levels[level]->interpolation.recv_buffers[neighbor] = all_recv_buffers;
      for(coarseBox=0;coarseBox<numCoarseBoxes;coarseBox++)if(coarseBoxes[coarseBox].sendRank==coarseRanks[neighbor]){
        // unpack MPI recv buffer...
        append_block_to_list(&(all_grids->levels[level]->interpolation.blocks[2]),&(all_grids->levels[level]->interpolation.allocated_blocks[2]),&(all_grids->levels[level]->interpolation.num_blocks[2]),
          /* dim.i         = */ all_grids->levels[level]->box_dim,
          /* dim.j         = */ all_grids->levels[level]->box_dim,
          /* dim.k         = */ all_grids->levels[level]->box_dim,
          /* read.box      = */ -1,
          /* read.ptr      = */ all_grids->levels[level]->interpolation.recv_buffers[neighbor],
          /* read.i        = */ offset,
          /* read.j        = */ 0,
          /* read.k        = */ 0,
          /* read.jStride  = */ all_grids->levels[level]->box_dim,
          /* read.kStride  = */ all_grids->levels[level]->box_dim*all_grids->levels[level]->box_dim,
          /* read.scale    = */ 1,
          /* write.box     = */ coarseBoxes[coarseBox].recvBox,
          /* write.ptr     = */ NULL,
          /* write.i       = */ 0,
          /* write.j       = */ 0,
          /* write.k       = */ 0,
          /* write.jStride = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].recvBox].jStride,
          /* write.kStride = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].recvBox].kStride,
          /* write.scale   = */ 1,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
        offset+=elementSize;
      }
      all_grids->levels[level]->interpolation.recv_ranks[neighbor] = coarseRanks[neighbor];
      all_grids->levels[level]->interpolation.recv_sizes[neighbor] = offset;
      all_recv_buffers+=offset;
    } // neighbor

    // free temporary storage...
    free(coarseBoxes);
    free(coarseRanks);
  } // recv/unpack


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  } // all levels


  #ifdef USE_MPI
  for(level=0;level<all_grids->num_levels;level++){
    all_grids->levels[level]->interpolation.requests = NULL;
    all_grids->levels[level]->interpolation.status   = NULL;
    if(level<all_grids->num_levels-1){  // i.e. bottom never calls interpolation()
    // by convention, level_f allocates a combined array of requests for both level_f recvs and level_c sends...
    int nMessages = all_grids->levels[level+1]->interpolation.num_sends + all_grids->levels[level]->interpolation.num_recvs;
    all_grids->levels[level]->interpolation.requests = (MPI_Request*)malloc(nMessages*sizeof(MPI_Request));
    all_grids->levels[level]->interpolation.status   = (MPI_Status *)malloc(nMessages*sizeof(MPI_Status ));
    }
  }
  #endif
}


//----------------------------------------------------------------------------------------------------------------------------------------------------
// build a list of operations and MPI buffers to affect distributed restriction
// the three lists constitute
//   - buffer packing (i.e. restrict a local box and place the result in an MPI buffer to be sent to a remote coarse grid process)
//   - local operations (i.e. restrict a local box and place the result in another local box or region of another local box)
//   - buffer upacking (i.e. copy restricted data recieved from another process into a local box or region of a local box)
void build_restriction(mg_type *all_grids, int restrictionType){
  int level;
  for(level=0;level<all_grids->num_levels;level++){

  // initialize to defaults...
  all_grids->levels[level]->restriction[restrictionType].num_recvs           = 0;
  all_grids->levels[level]->restriction[restrictionType].num_sends           = 0;
  all_grids->levels[level]->restriction[restrictionType].recv_ranks          = NULL;
  all_grids->levels[level]->restriction[restrictionType].send_ranks          = NULL;
  all_grids->levels[level]->restriction[restrictionType].recv_sizes          = NULL;
  all_grids->levels[level]->restriction[restrictionType].send_sizes          = NULL;
  all_grids->levels[level]->restriction[restrictionType].recv_buffers        = NULL;
  all_grids->levels[level]->restriction[restrictionType].send_buffers        = NULL;
  all_grids->levels[level]->restriction[restrictionType].blocks[0]           = NULL;
  all_grids->levels[level]->restriction[restrictionType].blocks[1]           = NULL;
  all_grids->levels[level]->restriction[restrictionType].blocks[2]           = NULL;
  all_grids->levels[level]->restriction[restrictionType].allocated_blocks[0] = 0;
  all_grids->levels[level]->restriction[restrictionType].allocated_blocks[1] = 0;
  all_grids->levels[level]->restriction[restrictionType].allocated_blocks[2] = 0;
  all_grids->levels[level]->restriction[restrictionType].num_blocks[0]       = 0; // number of unpack/insert operations  = number of boxes on level+1 that I don't own and restrict to 
  all_grids->levels[level]->restriction[restrictionType].num_blocks[1]       = 0; // number of unpack/insert operations  = number of boxes on level+1 that I own and restrict to
  all_grids->levels[level]->restriction[restrictionType].num_blocks[2]       = 0; // number of unpack/insert operations  = number of boxes on level-1 that I don't own that restrict to me
  #ifdef USE_MPI
  all_grids->levels[level]->restriction[restrictionType].requests            = NULL;
  all_grids->levels[level]->restriction[restrictionType].status              = NULL;
  #endif


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // construct pack, send, and local...
  if( (level<all_grids->num_levels-1) && (all_grids->levels[level]->num_my_boxes>0) ){ // not bottom  *and*  I have boxes to send

    // construct the list of coarsened boxes and neighboring ranks...
    int numCoarseBoxes = (all_grids->levels[level]->boxes_in.i/all_grids->levels[level+1]->boxes_in.i)*
                         (all_grids->levels[level]->boxes_in.j/all_grids->levels[level+1]->boxes_in.j)*
                         (all_grids->levels[level]->boxes_in.k/all_grids->levels[level+1]->boxes_in.k)*
                          all_grids->levels[level]->num_my_boxes;
        int *coarseRanks = (    int*)malloc(numCoarseBoxes*sizeof(    int)); // high water mark (assumes every neighboring box is a different process)
    RP_type *coarseBoxes = (RP_type*)malloc(numCoarseBoxes*sizeof(RP_type)); 
        numCoarseBoxes       = 0;
    int numCoarseBoxesLocal  = 0;
    int numCoarseBoxesRemote = 0;
    int fineBox;
    for(fineBox=0;fineBox<all_grids->levels[level]->num_my_boxes;fineBox++){
      int   fineBoxID = all_grids->levels[level]->my_boxes[fineBox].global_box_id;
      int   fineBox_i = all_grids->levels[level]->my_boxes[fineBox].low.i / all_grids->levels[level]->box_dim;
      int   fineBox_j = all_grids->levels[level]->my_boxes[fineBox].low.j / all_grids->levels[level]->box_dim;
      int   fineBox_k = all_grids->levels[level]->my_boxes[fineBox].low.k / all_grids->levels[level]->box_dim;
      int coarseBox_i = fineBox_i*all_grids->levels[level+1]->boxes_in.i/all_grids->levels[level]->boxes_in.i;
      int coarseBox_j = fineBox_j*all_grids->levels[level+1]->boxes_in.j/all_grids->levels[level]->boxes_in.j;
      int coarseBox_k = fineBox_k*all_grids->levels[level+1]->boxes_in.k/all_grids->levels[level]->boxes_in.k;
      int coarseBoxID =  coarseBox_i + coarseBox_j*all_grids->levels[level+1]->boxes_in.i + coarseBox_k*all_grids->levels[level+1]->boxes_in.i*all_grids->levels[level+1]->boxes_in.j;
      int coarseBox   = -1;int c;for(c=0;c<all_grids->levels[level+1]->num_my_boxes;c++)if( all_grids->levels[level+1]->my_boxes[c].global_box_id == coarseBoxID )coarseBox=c; // try and find the coarseBox index of a box with global_box_id == coaseBoxID
      coarseBoxes[numCoarseBoxes].sendRank  = all_grids->levels[level  ]->rank_of_box[  fineBoxID];
      coarseBoxes[numCoarseBoxes].sendBoxID = fineBoxID;
      coarseBoxes[numCoarseBoxes].sendBox   = fineBox;
      coarseBoxes[numCoarseBoxes].recvRank  = all_grids->levels[level+1]->rank_of_box[coarseBoxID];
      coarseBoxes[numCoarseBoxes].recvBoxID = coarseBoxID;
      coarseBoxes[numCoarseBoxes].recvBox   = coarseBox;  // -1 if off-node
      coarseBoxes[numCoarseBoxes].i         = (all_grids->levels[level]->box_dim/2)*( fineBox_i % (all_grids->levels[level]->boxes_in.i/all_grids->levels[level+1]->boxes_in.i) );
      coarseBoxes[numCoarseBoxes].j         = (all_grids->levels[level]->box_dim/2)*( fineBox_j % (all_grids->levels[level]->boxes_in.j/all_grids->levels[level+1]->boxes_in.j) );
      coarseBoxes[numCoarseBoxes].k         = (all_grids->levels[level]->box_dim/2)*( fineBox_k % (all_grids->levels[level]->boxes_in.k/all_grids->levels[level+1]->boxes_in.k) );
                  numCoarseBoxes++;
      if(all_grids->levels[level]->my_rank != all_grids->levels[level+1]->rank_of_box[coarseBoxID]){
        coarseRanks[numCoarseBoxesRemote++] = all_grids->levels[level+1]->rank_of_box[coarseBoxID];
      }else{numCoarseBoxesLocal++;}
    } // my (fine) boxes

    // sort boxes by sendRank(==my rank) then by sendBoxID... ensures the sends and receive buffers are always sorted by sendBoxID...
    qsort(coarseBoxes,numCoarseBoxes      ,sizeof(RP_type),qsortRP );
    // sort the lists of neighboring ranks and remove duplicates...
    qsort(coarseRanks,numCoarseBoxesRemote,sizeof(    int),qsortInt);
    int numCoarseRanks=0;
    int _rank=-1;int neighbor=0;
    for(neighbor=0;neighbor<numCoarseBoxesRemote;neighbor++)if(coarseRanks[neighbor] != _rank){_rank=coarseRanks[neighbor];coarseRanks[numCoarseRanks++]=coarseRanks[neighbor];}

    // allocate structures...
    all_grids->levels[level]->restriction[restrictionType].num_sends     =                         numCoarseRanks;
    all_grids->levels[level]->restriction[restrictionType].send_ranks    =            (int*)malloc(numCoarseRanks*sizeof(int));
    all_grids->levels[level]->restriction[restrictionType].send_sizes    =            (int*)malloc(numCoarseRanks*sizeof(int));
    all_grids->levels[level]->restriction[restrictionType].send_buffers  =        (double**)malloc(numCoarseRanks*sizeof(double*));
    if(numCoarseRanks>0){
    if(all_grids->levels[level]->restriction[restrictionType].send_ranks  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].send_ranks\n",level);exit(0);}
    if(all_grids->levels[level]->restriction[restrictionType].send_sizes  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].send_sizes\n",level);exit(0);}
    if(all_grids->levels[level]->restriction[restrictionType].send_buffers==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].send_buffers\n",level);exit(0);}
    }

    int elementSize;
    int restrict_dim_i=-1;
    int restrict_dim_j=-1;
    int restrict_dim_k=-1;
    switch(restrictionType){
      case RESTRICT_CELL   : restrict_dim_i = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level]->box_dim/2);break;
      case RESTRICT_FACE_I : restrict_dim_i = (1+all_grids->levels[level]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level]->box_dim/2);break;
      case RESTRICT_FACE_J : restrict_dim_i = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_j = (1+all_grids->levels[level]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level]->box_dim/2);break;
      case RESTRICT_FACE_K : restrict_dim_i = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_k = (1+all_grids->levels[level]->box_dim/2);break;
    }
    elementSize = restrict_dim_i*restrict_dim_j*restrict_dim_k;
   
    double * all_send_buffers = (double*)malloc(numCoarseBoxesRemote*elementSize*sizeof(double));
          if(numCoarseBoxesRemote*elementSize>0)
          if(all_send_buffers==NULL){fprintf(stderr,"malloc failed - restriction/all_send_buffers\n");exit(0);}
                      memset(all_send_buffers,0,numCoarseBoxesRemote*elementSize*sizeof(double)); // DO NOT DELETE... you must initialize to 0 to avoid getting something like 0.0*NaN and corrupting the solve

    // for each neighbor, construct the pack list and allocate the MPI send buffer... 
    for(neighbor=0;neighbor<numCoarseRanks;neighbor++){
      int coarseBox;
      int offset = 0;
      all_grids->levels[level]->restriction[restrictionType].send_buffers[neighbor] = all_send_buffers;
      for(coarseBox=0;coarseBox<numCoarseBoxes;coarseBox++)if(coarseBoxes[coarseBox].recvRank==coarseRanks[neighbor]){
        // restrict to MPI send buffer...
        append_block_to_list( &(all_grids->levels[level]->restriction[restrictionType].blocks[0]),
                              &(all_grids->levels[level]->restriction[restrictionType].allocated_blocks[0]),
                              &(all_grids->levels[level]->restriction[restrictionType].num_blocks[0]),
          /* dim.i         = */ restrict_dim_i, 
          /* dim.j         = */ restrict_dim_j, 
          /* dim.k         = */ restrict_dim_k, 
          /* read.box      = */ coarseBoxes[coarseBox].sendBox,
          /* read.ptr      = */ NULL,
          /* read.i        = */ 0,
          /* read.j        = */ 0,
          /* read.k        = */ 0,
          /* read.jStride  = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].sendBox].jStride,
          /* read.kStride  = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].sendBox].kStride,
          /* read.scale    = */ 2,
          /* write.box     = */ -1,
          /* write.ptr     = */ all_grids->levels[level]->restriction[restrictionType].send_buffers[neighbor],
          /* write.i       = */ offset,
          /* write.j       = */ 0,
          /* write.k       = */ 0,
          /* write.jStride = */ restrict_dim_i,
          /* write.kStride = */ restrict_dim_i*restrict_dim_j, 
          /* write.scale   = */ 1,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
        offset+=elementSize;
      }
      all_grids->levels[level]->restriction[restrictionType].send_ranks[neighbor] = coarseRanks[neighbor];
      all_grids->levels[level]->restriction[restrictionType].send_sizes[neighbor] = offset;
      all_send_buffers+=offset;
    }
    // for construct the local restriction list... 
    {
      int coarseBox;
      for(coarseBox=0;coarseBox<numCoarseBoxes;coarseBox++)if(coarseBoxes[coarseBox].recvRank==all_grids->levels[level+1]->my_rank){
        // restrict to local...
        append_block_to_list( &(all_grids->levels[level]->restriction[restrictionType].blocks[1]),
                              &(all_grids->levels[level]->restriction[restrictionType].allocated_blocks[1]),
                              &(all_grids->levels[level]->restriction[restrictionType].num_blocks[1]),
          /* dim.i         = */ restrict_dim_i, 
          /* dim.j         = */ restrict_dim_j, 
          /* dim.k         = */ restrict_dim_k, 
          /* read.box      = */ coarseBoxes[coarseBox].sendBox,
          /* read.ptr      = */ NULL,
          /* read.i        = */ 0, 
          /* read.j        = */ 0,
          /* read.k        = */ 0,
          /* read.jStride  = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].sendBox].jStride,
          /* read.kStride  = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].sendBox].kStride,
          /* read.scale    = */ 2,
          /* write.box     = */ coarseBoxes[coarseBox].recvBox,
          /* write.ptr     = */ NULL,
          /* write.i       = */ coarseBoxes[coarseBox].i,
          /* write.j       = */ coarseBoxes[coarseBox].j,
          /* write.k       = */ coarseBoxes[coarseBox].k,
          /* write.jStride = */ all_grids->levels[level+1]->my_boxes[coarseBoxes[coarseBox].recvBox].jStride,
          /* write.kStride = */ all_grids->levels[level+1]->my_boxes[coarseBoxes[coarseBox].recvBox].kStride,
          /* write.scale   = */ 1,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
      }
    } // local to local

    // free temporary storage...
    free(coarseBoxes);
    free(coarseRanks);
  } // send/pack/local




  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // construct recv and unpack...
  if( (level>0) && (all_grids->levels[level]->num_my_boxes>0) ){ // not top  *and*  I have boxes to receive
    // construct a list of fine boxes to be coarsened and sent to me...
    int numFineBoxesMax = (all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i)*
                          (all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j)*
                          (all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k)*
                                                                  all_grids->levels[level]->num_my_boxes;
        int *fineRanks = (    int*)malloc(numFineBoxesMax*sizeof(    int)); // high water mark (assumes every neighboring box is a different process)
    RP_type *fineBoxes = (RP_type*)malloc(numFineBoxesMax*sizeof(RP_type)); 
    int numFineBoxesRemote = 0;
    int coarseBox;
    for(coarseBox=0;coarseBox<all_grids->levels[level]->num_my_boxes;coarseBox++){
      int bi,bj,bk;
      int   coarseBoxID = all_grids->levels[level]->my_boxes[coarseBox].global_box_id;
      int   coarseBox_i = all_grids->levels[level]->my_boxes[coarseBox].low.i / all_grids->levels[level]->box_dim;
      int   coarseBox_j = all_grids->levels[level]->my_boxes[coarseBox].low.j / all_grids->levels[level]->box_dim;
      int   coarseBox_k = all_grids->levels[level]->my_boxes[coarseBox].low.k / all_grids->levels[level]->box_dim;
      for(bk=0;bk<all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k;bk++){
      for(bj=0;bj<all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j;bj++){
      for(bi=0;bi<all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i;bi++){
        int fineBox_i = (all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i)*coarseBox_i + bi;
        int fineBox_j = (all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j)*coarseBox_j + bj;
        int fineBox_k = (all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k)*coarseBox_k + bk;
        int fineBoxID =  fineBox_i + fineBox_j*all_grids->levels[level-1]->boxes_in.i + fineBox_k*all_grids->levels[level-1]->boxes_in.i*all_grids->levels[level-1]->boxes_in.j;
        if(all_grids->levels[level-1]->rank_of_box[fineBoxID] != all_grids->levels[level]->my_rank){
          fineBoxes[numFineBoxesRemote].sendRank  = all_grids->levels[level-1]->rank_of_box[  fineBoxID];
          fineBoxes[numFineBoxesRemote].sendBoxID = fineBoxID;
          fineBoxes[numFineBoxesRemote].sendBox   = -1; // I don't know the off-node box index
          fineBoxes[numFineBoxesRemote].recvRank  = all_grids->levels[level  ]->rank_of_box[coarseBoxID];
          fineBoxes[numFineBoxesRemote].recvBoxID = coarseBoxID;
          fineBoxes[numFineBoxesRemote].recvBox   = coarseBox;
          fineBoxes[numFineBoxesRemote].i         = bi*all_grids->levels[level-1]->box_dim/2;
          fineBoxes[numFineBoxesRemote].j         = bj*all_grids->levels[level-1]->box_dim/2;
          fineBoxes[numFineBoxesRemote].k         = bk*all_grids->levels[level-1]->box_dim/2;
          fineRanks[numFineBoxesRemote] = all_grids->levels[level-1]->rank_of_box[fineBoxID];
                    numFineBoxesRemote++;
        }
      }}}
    } // my (coarse) boxes
    // sort boxes by sendRank(==my rank) then by sendBoxID... ensures the sends and receive buffers are always sorted by sendBoxID...
    qsort(fineBoxes,numFineBoxesRemote,sizeof(RP_type),qsortRP );
    // sort the lists of neighboring ranks and remove duplicates...
    qsort(fineRanks,numFineBoxesRemote,sizeof(    int),qsortInt);
    int numFineRanks=0;
    int _rank=-1;int neighbor=0;
    for(neighbor=0;neighbor<numFineBoxesRemote;neighbor++)if(fineRanks[neighbor] != _rank){_rank=fineRanks[neighbor];fineRanks[numFineRanks++]=fineRanks[neighbor];}

    // allocate structures...
    all_grids->levels[level]->restriction[restrictionType].num_recvs     =                         numFineRanks;
    all_grids->levels[level]->restriction[restrictionType].recv_ranks    =            (int*)malloc(numFineRanks*sizeof(int));
    all_grids->levels[level]->restriction[restrictionType].recv_sizes    =            (int*)malloc(numFineRanks*sizeof(int));
    all_grids->levels[level]->restriction[restrictionType].recv_buffers  =        (double**)malloc(numFineRanks*sizeof(double*));
    if(numFineRanks>0){
    if(all_grids->levels[level]->restriction[restrictionType].recv_ranks  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].recv_ranks  \n",level);exit(0);}
    if(all_grids->levels[level]->restriction[restrictionType].recv_sizes  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].recv_sizes  \n",level);exit(0);}
    if(all_grids->levels[level]->restriction[restrictionType].recv_buffers==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].recv_buffers\n",level);exit(0);}
    }

    int elementSize;
    int restrict_dim_i=-1;
    int restrict_dim_j=-1;
    int restrict_dim_k=-1;
    switch(restrictionType){
      case RESTRICT_CELL   : restrict_dim_i = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level-1]->box_dim/2);break;
      case RESTRICT_FACE_I : restrict_dim_i = (1+all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level-1]->box_dim/2);break;
      case RESTRICT_FACE_J : restrict_dim_i = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_j = (1+all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level-1]->box_dim/2);break;
      case RESTRICT_FACE_K : restrict_dim_i = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_k = (1+all_grids->levels[level-1]->box_dim/2);break;
    }
    elementSize = restrict_dim_i*restrict_dim_j*restrict_dim_k;

    double * all_recv_buffers = (double*)malloc(numFineBoxesRemote*elementSize*sizeof(double));
          if(numFineBoxesRemote*elementSize>0)
          if(all_recv_buffers==NULL){fprintf(stderr,"malloc failed - restriction/all_recv_buffers\n");exit(0);}
                      memset(all_recv_buffers,0,numFineBoxesRemote*elementSize*sizeof(double)); // DO NOT DELETE... you must initialize to 0 to avoid getting something like 0.0*NaN and corrupting the solve
    //printf("level=%d, rank=%2d, recv_buffers=%6d\n",level,all_grids->my_rank,numFineBoxesRemote*elementSize*sizeof(double));

    // for each neighbor, construct the unpack list and allocate the MPI recv buffer... 
    for(neighbor=0;neighbor<numFineRanks;neighbor++){
      int fineBox;
      int offset = 0;
      all_grids->levels[level]->restriction[restrictionType].recv_buffers[neighbor] = all_recv_buffers;
      for(fineBox=0;fineBox<numFineBoxesRemote;fineBox++)if(fineBoxes[fineBox].sendRank==fineRanks[neighbor]){
        // unpack MPI recv buffer...
        append_block_to_list( &(all_grids->levels[level]->restriction[restrictionType].blocks[2]),
                              &(all_grids->levels[level]->restriction[restrictionType].allocated_blocks[2]),
                              &(all_grids->levels[level]->restriction[restrictionType].num_blocks[2]),
          /* dim.i         = */ restrict_dim_i, 
          /* dim.j         = */ restrict_dim_j, 
          /* dim.k         = */ restrict_dim_k, 
          /* read.box      = */ -1,
          /* read.ptr      = */ all_grids->levels[level]->restriction[restrictionType].recv_buffers[neighbor],
          /* read.i        = */ offset,
          /* read.j        = */ 0,
          /* read.k        = */ 0,
          /* read.jStride  = */ restrict_dim_i,
          /* read.kStride  = */ restrict_dim_i*restrict_dim_j, 
          /* read.scale    = */ 1,
          /* write.box     = */ fineBoxes[fineBox].recvBox,
          /* write.ptr     = */ NULL,
          /* write.i       = */ fineBoxes[fineBox].i,
          /* write.j       = */ fineBoxes[fineBox].j,
          /* write.k       = */ fineBoxes[fineBox].k,
          /* write.jStride = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].recvBox].jStride,
          /* write.kStride = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].recvBox].kStride,
          /* write.scale   = */ 1,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
        offset+=elementSize;
      }
      all_grids->levels[level]->restriction[restrictionType].recv_ranks[neighbor] = fineRanks[neighbor];
      all_grids->levels[level]->restriction[restrictionType].recv_sizes[neighbor] = offset;
      all_recv_buffers+=offset;
    } // neighbor

    // free temporary storage...
    free(fineBoxes);
    free(fineRanks);
  } // recv/unpack



  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  } // level loop


  #ifdef USE_MPI
  for(level=0;level<all_grids->num_levels;level++){
    all_grids->levels[level]->restriction[restrictionType].requests = NULL;
    all_grids->levels[level]->restriction[restrictionType].status   = NULL;
    if(level<all_grids->num_levels-1){ // bottom never calls restriction()
    // by convention, level_f allocates a combined array of requests for both level_f sends and level_c recvs...
    int nMessages = all_grids->levels[level+1]->restriction[restrictionType].num_recvs + all_grids->levels[level]->restriction[restrictionType].num_sends;
    all_grids->levels[level]->restriction[restrictionType].requests = (MPI_Request*)malloc(nMessages*sizeof(MPI_Request));
    all_grids->levels[level]->restriction[restrictionType].status   = (MPI_Status *)malloc(nMessages*sizeof(MPI_Status ));
    }
  }
  #endif
}


//------------------------------------------------------------------------------------------------------------------------------
// given a fine grid input, build a hiearchy of MG levels
// level 0 simply points to fine_grid.  All other levels are created
// rebuild the restriction/interpolation lists for each coarse grid level
// rebuild the operator on each coarse grid level
// add extra vectors to the coarse grid once here instead of on every call to the coarse grid solve
// NOTE, this routine presumes the fine_grid domain is cubical... fine_grid->dim.i==fine_grid->dim.j==fine_grid->dim.k
// NOTE, as this function is not timed, it has not been optimzied for performance
void MGBuild(mg_type *all_grids, level_type *fine_grid, double a, double b, int minCoarseGridDim, const MPI_Comm comm){
  int  maxLevels=100; // i.e. maximum problem size is (2^100)^3
  int     nProcs[100];
  int      dim_i[100];
  int boxes_in_i[100];
  int    box_dim[100];
  int box_ghosts[100];
  all_grids->my_rank = fine_grid->my_rank;
  all_grids->timers.MGBuild = 0;
  double _timeStartMGBuild = getTime();

  // calculate how deep we can make the v-cycle...
  int level=1;
                             int coarse_dim = fine_grid->dim.i;
//if(fine_grid->dim.j<coarse_dim)coarse_dim = fine_grid->dim.j;
//if(fine_grid->dim.k<coarse_dim)coarse_dim = fine_grid->dim.k;
  while( (coarse_dim>=2*minCoarseGridDim) && ((coarse_dim&0x1)==0) ){ // grid dimension is even and big enough...
    level++;
    coarse_dim = coarse_dim / 2;
  }if(level<maxLevels)maxLevels=level;

      nProcs[0] = fine_grid->num_ranks;
       dim_i[0] = fine_grid->dim.i;
  boxes_in_i[0] = fine_grid->boxes_in.i;
     box_dim[0] = fine_grid->box_dim;
  box_ghosts[0] = fine_grid->box_ghosts;

  // build the list of levels...
  all_grids->levels = (level_type**)malloc(maxLevels*sizeof(level_type*));
  if(all_grids->levels == NULL){fprintf(stderr,"malloc failed - MGBuild/all_grids->levels\n");exit(0);}
  all_grids->num_levels=1;
  all_grids->levels[0] = fine_grid;


  // build a table to guide the construction of the v-cycle...
  int doRestrict=1;if(maxLevels<2)doRestrict=0; // i.e. can't restrict if there is only one level !!!
  #ifdef USE_UCYCLES
  while(doRestrict){
    level = all_grids->num_levels;
    doRestrict=0;
    if( (box_dim[level-1] % 2 == 0) ){
          nProcs[level] =     nProcs[level-1];
           dim_i[level] =      dim_i[level-1]/2;
         box_dim[level] =    box_dim[level-1]/2;
      boxes_in_i[level] = boxes_in_i[level-1];
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }
    if(box_dim[level] < box_ghosts[level])doRestrict=0;
    if(dim_i[level]<minCoarseGridDim)doRestrict=0;
    if(doRestrict)all_grids->num_levels++;
  }
  #else // TRUE V-Cycle...
  while(doRestrict){
    level = all_grids->num_levels;
    doRestrict=0;
    int fine_box_dim    =    box_dim[level-1];
    int fine_nProcs     =     nProcs[level-1];
    int fine_dim_i      =      dim_i[level-1];
    int fine_boxes_in_i = boxes_in_i[level-1];
    if( (fine_box_dim % 2 == 0) && (fine_box_dim > MG_AGGLOMERATION_START) && ((fine_box_dim/2)>=stencil_get_radius()) ){ // Boxes are too big to agglomerate
          nProcs[level] = fine_nProcs;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_box_dim/2; // FIX, verify its not less than the stencil radius
      boxes_in_i[level] = fine_boxes_in_i;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }else
    if( (fine_boxes_in_i % 2 == 0) && ((fine_box_dim)>=stencil_get_radius()) ){ // 8:1 box agglomeration
          nProcs[level] = fine_nProcs;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_box_dim;
      boxes_in_i[level] = fine_boxes_in_i/2;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }else
    if( (coarse_dim != 1) && (fine_dim_i == 2*coarse_dim) && ((fine_dim_i/2)>=stencil_get_radius()) ){ // agglomerate everything
          nProcs[level] = 1;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_dim_i/2; // FIX, verify its not less than the stencil radius
      boxes_in_i[level] = 1;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }else
    if( (coarse_dim != 1) && (fine_dim_i == 4*coarse_dim) && ((fine_box_dim/2)>=stencil_get_radius()) ){ // restrict box dimension, and run on fewer ranks
          nProcs[level] = coarse_dim<fine_nProcs ? coarse_dim : fine_nProcs;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_box_dim/2; // FIX, verify its not less than the stencil radius
      boxes_in_i[level] = fine_boxes_in_i;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }else
    if( (coarse_dim != 1) && (fine_dim_i == 8*coarse_dim) && ((fine_box_dim/2)>=stencil_get_radius()) ){ // restrict box dimension, and run on fewer ranks
          nProcs[level] = coarse_dim*coarse_dim<fine_nProcs ? coarse_dim*coarse_dim : fine_nProcs;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_box_dim/2; // FIX, verify its not less than the stencil radius
      boxes_in_i[level] = fine_boxes_in_i;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }else
    if( (fine_box_dim % 2 == 0) && ((fine_box_dim/2)>=stencil_get_radius()) ){ // restrict box dimension, and run on the same number of ranks
          nProcs[level] = fine_nProcs;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_box_dim/2; // FIX, verify its not less than the stencil radius
      boxes_in_i[level] = fine_boxes_in_i;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }
    if(dim_i[level]<minCoarseGridDim)doRestrict=0;
    if(doRestrict)all_grids->num_levels++;
  }
  #endif


  // now build all the coarsened levels...
  for(level=1;level<all_grids->num_levels;level++){
    all_grids->levels[level] = (level_type*)malloc(sizeof(level_type));
    if(all_grids->levels[level] == NULL){fprintf(stderr,"malloc failed - MGBuild/doRestrict\n");exit(0);}
    create_level(all_grids->levels[level],boxes_in_i[level],box_dim[level],box_ghosts[level],all_grids->levels[level-1]->numVectors,all_grids->levels[level-1]->boundary_condition.type,all_grids->levels[level-1]->my_rank,nProcs[level], comm);
    all_grids->levels[level]->h = 2.0*all_grids->levels[level-1]->h;
  }


  // bottom solver (level = all_grids->num_levels-1) gets extra vectors...
  create_vectors(all_grids->levels[all_grids->num_levels-1],all_grids->levels[all_grids->num_levels-1]->numVectors + IterativeSolver_NumVectors() );


  // build the restriction and interpolation communicators...
  if(all_grids->my_rank==0){fprintf(stdout,"\n  Building restriction and interpolation lists... ");fflush(stdout);}
  build_restriction(all_grids,RESTRICT_CELL  ); // cell-centered
  build_restriction(all_grids,RESTRICT_FACE_I); // face-centered, normal to i
  build_restriction(all_grids,RESTRICT_FACE_J); // face-centered, normal to j
  build_restriction(all_grids,RESTRICT_FACE_K); // face-centered, normal to k
  build_interpolation(all_grids);
  if(all_grids->my_rank==0){fprintf(stdout,"done\n");fflush(stdout);}


  // build subcommunicators...
  #ifdef USE_MPI
  #ifdef USE_SUBCOMM
  if(all_grids->my_rank==0){fprintf(stdout,"\n");}
  for(level=1;level<all_grids->num_levels;level++){
    double comm_split_start = MPI_Wtime();
    if(all_grids->my_rank==0){fprintf(stdout,"  Building MPI subcommunicator for level %d... ",level);fflush(stdout);}
    all_grids->levels[level]->active=0;
    int ll;for(ll=level;ll<all_grids->num_levels;ll++)if(all_grids->levels[ll]->num_my_boxes>0)all_grids->levels[level]->active=1;
    MPI_Comm_split(comm, all_grids->levels[level]->active, all_grids->levels[level]->my_rank, &all_grids->levels[level]->MPI_COMM_ALLREDUCE);
    double comm_split_end = MPI_Wtime();
    double comm_split_time_send = comm_split_end-comm_split_start;
    double comm_split_time = 0;
    MPI_Allreduce(&comm_split_time_send,&comm_split_time,1,MPI_DOUBLE,MPI_MAX,all_grids->levels[level]->MPI_COMM_ALLREDUCE);
    if(all_grids->my_rank==0){fprintf(stdout,"done (%0.6f seconds)\n",comm_split_time);fflush(stdout);}
  }
  #endif
  #endif


  // rebuild various coefficients for the operator... must occur after build_restriction !!!
  if(all_grids->my_rank==0){fprintf(stdout,"\n");}
  for(level=1;level<all_grids->num_levels;level++){
    rebuild_operator(all_grids->levels[level],(level>0)?all_grids->levels[level-1]:NULL,a,b);
  }
  if(all_grids->my_rank==0){fprintf(stdout,"\n");}


  // quick tests for Poisson, Neumann, etc...
  for(level=0;level<all_grids->num_levels;level++){
    all_grids->levels[level]->must_subtract_mean = 0;
    int alpha_is_zero = (dot(all_grids->levels[level],VECTOR_ALPHA,VECTOR_ALPHA) == 0.0);
    // For Poisson with Periodic Boundary Conditions, by convention we assume the solution sums to zero.  Eliminate any constants from the solution by subtracting the mean.
    if( (all_grids->levels[level]->boundary_condition.type==BC_PERIODIC) && ((a==0) || (alpha_is_zero==1)) )all_grids->levels[level]->must_subtract_mean = 1;
  }

  
  all_grids->timers.MGBuild += (double)(getTime()-_timeStartMGBuild);
}


//------------------------------------------------------------------------------------------------------------------------------
// deallocate all memory created in the MG hierarchy
// WARNING, this will free the fine_grid level as well (FIX?)
void MGDestroy(mg_type *all_grids){
  int level;
  int i;

  #ifdef USE_MPI
  #ifdef USE_SUBCOMM
  // only MGBuild creates subcommunicators (level_create assigns)
  for(level=all_grids->num_levels-1;level>0;level--){
    if(all_grids->levels[level]->MPI_COMM_ALLREDUCE != MPI_COMM_WORLD)
    MPI_Comm_free(&all_grids->levels[level]->MPI_COMM_ALLREDUCE);
  }
  #endif
  #endif

  if(all_grids->my_rank==0){fprintf(stdout,"attempting to free the restriction and interpolation lists... ");fflush(stdout);}
  for(level=all_grids->num_levels-1;level>=0;level--){
    // destroy restriction mini program created by MGBuild...
    for(i=0;i<4;i++){
      if(all_grids->levels[level]->restriction[i].num_recvs>0){
      //for(j=0;j<all_grids->levels[level]->restriction[i].num_recvs;j++)if(all_grids->levels[level]->restriction[i].recv_buffers[j])free(all_grids->levels[level]->restriction[i].recv_buffers[j]);
      if(all_grids->levels[level]->restriction[i].recv_buffers[0])free(all_grids->levels[level]->restriction[i].recv_buffers[0]); // allocated in bulk
      if(all_grids->levels[level]->restriction[i].recv_buffers   )free(all_grids->levels[level]->restriction[i].recv_buffers   );
      if(all_grids->levels[level]->restriction[i].recv_ranks     )free(all_grids->levels[level]->restriction[i].recv_ranks     );
      if(all_grids->levels[level]->restriction[i].recv_sizes     )free(all_grids->levels[level]->restriction[i].recv_sizes     );
      }
      if(all_grids->levels[level]->restriction[i].num_sends>0){
      //for(j=0;j<all_grids->levels[level]->restriction[i].num_sends;j++)if(all_grids->levels[level]->restriction[i].send_buffers[j])free(all_grids->levels[level]->restriction[i].send_buffers[j]);
      if(all_grids->levels[level]->restriction[i].send_buffers[0])free(all_grids->levels[level]->restriction[i].send_buffers[0]); // allocated in bulk
      if(all_grids->levels[level]->restriction[i].send_buffers   )free(all_grids->levels[level]->restriction[i].send_buffers   );
      if(all_grids->levels[level]->restriction[i].send_ranks     )free(all_grids->levels[level]->restriction[i].send_ranks     );
      if(all_grids->levels[level]->restriction[i].send_sizes     )free(all_grids->levels[level]->restriction[i].send_sizes     );
      }
      if(all_grids->levels[level]->restriction[i].blocks[0]      )free(all_grids->levels[level]->restriction[i].blocks[0]      );
      if(all_grids->levels[level]->restriction[i].blocks[1]      )free(all_grids->levels[level]->restriction[i].blocks[1]      );
      if(all_grids->levels[level]->restriction[i].blocks[2]      )free(all_grids->levels[level]->restriction[i].blocks[2]      );
      #ifdef USE_MPI
      if(all_grids->levels[level]->restriction[i].requests       )free(all_grids->levels[level]->restriction[i].requests       );
      if(all_grids->levels[level]->restriction[i].status         )free(all_grids->levels[level]->restriction[i].status         );
      #endif
    }

    // destroy interpolation mini program created by MGBuild...
    if(all_grids->levels[level]->interpolation.num_recvs>0){
    //for(j=0;j<all_grids->levels[level]->interpolation.num_recvs;j++)if(all_grids->levels[level]->interpolation.recv_buffers[j])free(all_grids->levels[level]->interpolation.recv_buffers[j]);
    if(all_grids->levels[level]->interpolation.recv_buffers[0])free(all_grids->levels[level]->interpolation.recv_buffers[0]); // allocated in bulk
    if(all_grids->levels[level]->interpolation.recv_buffers   )free(all_grids->levels[level]->interpolation.recv_buffers   );
    if(all_grids->levels[level]->interpolation.recv_ranks     )free(all_grids->levels[level]->interpolation.recv_ranks     );
    if(all_grids->levels[level]->interpolation.recv_sizes     )free(all_grids->levels[level]->interpolation.recv_sizes     );
    }
    if(all_grids->levels[level]->interpolation.num_sends>0){
    //for(j=0;j<all_grids->levels[level]->interpolation.num_sends;j++)if(all_grids->levels[level]->interpolation.send_buffers[j])free(all_grids->levels[level]->interpolation.send_buffers[j]);
    if(all_grids->levels[level]->interpolation.send_buffers[0])free(all_grids->levels[level]->interpolation.send_buffers[0]); // allocated in bulk
    if(all_grids->levels[level]->interpolation.send_buffers   )free(all_grids->levels[level]->interpolation.send_buffers   );
    if(all_grids->levels[level]->interpolation.send_ranks     )free(all_grids->levels[level]->interpolation.send_ranks     );
    if(all_grids->levels[level]->interpolation.send_sizes     )free(all_grids->levels[level]->interpolation.send_sizes     );
    }
    if(all_grids->levels[level]->interpolation.blocks[0]      )free(all_grids->levels[level]->interpolation.blocks[0]      );
    if(all_grids->levels[level]->interpolation.blocks[1]      )free(all_grids->levels[level]->interpolation.blocks[1]      );
    if(all_grids->levels[level]->interpolation.blocks[2]      )free(all_grids->levels[level]->interpolation.blocks[2]      );
    #ifdef USE_MPI
    if(all_grids->levels[level]->interpolation.requests       )free(all_grids->levels[level]->interpolation.requests       );
    if(all_grids->levels[level]->interpolation.status         )free(all_grids->levels[level]->interpolation.status         );
    #endif

  }
  if(all_grids->my_rank==0){fprintf(stdout,"done\n");}

  // now destroy the level itself (but don't destroy level 0 as it was not created by MGBuild)
  for(level=all_grids->num_levels-1;level>0;level--){
    destroy_level(all_grids->levels[level]);
  }
  if(all_grids->levels)free(all_grids->levels);
}


//------------------------------------------------------------------------------------------------------------------------------
// perform a richardson error analysis to infer the order of the operator/solver
void richardson_error(mg_type *all_grids, int levelh, int u_id){
  // in FV...
  // +-------+   +---+---+   +-------+   +-------+
  // |       |   | a | b |   |       |   |a+b+c+d|
  // |  u^2h | - +---+---+ = |  u^2h | - |  ---  |
  // |       |   | c | d |   |       |   |   4   |
  // +-------+   +---+---+   +-------+   +-------+
  //
  restriction(all_grids->levels[levelh+1],VECTOR_TEMP,all_grids->levels[levelh  ],u_id,RESTRICT_CELL); // temp^2h = R u^h
  restriction(all_grids->levels[levelh+2],VECTOR_TEMP,all_grids->levels[levelh+1],u_id,RESTRICT_CELL); // temp^4h = R u^2h
  add_vectors(all_grids->levels[levelh+1],VECTOR_TEMP,1.0,u_id,-1.0,VECTOR_TEMP);                      // temp^2h = u^2h - temp^2h = u^2h - R u^h
  add_vectors(all_grids->levels[levelh+2],VECTOR_TEMP,1.0,u_id,-1.0,VECTOR_TEMP);                      // temp^2h = u^4h - temp^4h = u^4h - R u^2h
  double norm_of_u2h_minus_uh  = norm(all_grids->levels[levelh+1],VECTOR_TEMP); // || u^2h - R u^h  ||max
  double norm_of_u4h_minus_u2h = norm(all_grids->levels[levelh+2],VECTOR_TEMP); // || u^4h - R u^2h ||max
  // estimate the error^h using ||u^2h - R u^h||
  if(all_grids->my_rank==0){fprintf(stdout,"  h=%0.15e  ||error||=%0.15e\n",all_grids->levels[levelh]->h,norm_of_u2h_minus_uh);fflush(stdout);}
  // log( ||u^4h - R u^2h|| / ||u^2h - R u^h|| ) / log(2) is an estimate of the order of the method (e.g. 4th order)
  if(all_grids->my_rank==0){fprintf(stdout,"  order=%0.3f\n",log(norm_of_u4h_minus_u2h / norm_of_u2h_minus_uh) / log(2) );fflush(stdout);}
}


//------------------------------------------------------------------------------------------------------------------------------
void MGVCycle(mg_type *all_grids, int e_id, int R_id, double a, double b, int level){
  if(!all_grids->levels[level]->active)return;
  double _LevelStart;

  // bottom solve...
  if(level==all_grids->num_levels-1){
    double _timeBottomStart = getTime();
    IterativeSolver(all_grids->levels[level],e_id,R_id,a,b,MG_DEFAULT_BOTTOM_NORM);
    all_grids->levels[level]->timers.Total += (double)(getTime()-_timeBottomStart);
    return;
  }

  // down...
  _LevelStart = getTime();
       smooth(all_grids->levels[level  ],e_id,R_id,a,b);
     residual(all_grids->levels[level  ],VECTOR_TEMP,e_id,R_id,a,b);
  restriction(all_grids->levels[level+1],R_id,all_grids->levels[level],VECTOR_TEMP,RESTRICT_CELL);
  zero_vector(all_grids->levels[level+1],e_id);
  all_grids->levels[level]->timers.Total += (double)(getTime()-_LevelStart);

  // recursion...
  MGVCycle(all_grids,e_id,R_id,a,b,level+1);

  // up...
  _LevelStart = getTime();
  interpolation_vcycle(all_grids->levels[level  ],e_id,1.0,all_grids->levels[level+1],e_id);
                smooth(all_grids->levels[level  ],e_id,R_id,a,b);

  all_grids->levels[level]->timers.Total += (double)(getTime()-_LevelStart);
}


//------------------------------------------------------------------------------------------------------------------------------
void MGSolve(mg_type *all_grids, int onLevel, int u_id, int F_id, double a, double b, double dtol, double rtol){
  // solves Au=f on level 'onLevel'
  all_grids->MGSolves_performed++;
  if(!all_grids->levels[onLevel]->active)return;
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int e_id = u_id; // __u FIX
  int R_id = VECTOR_R;
  int v;
  int maxVCycles = 20;

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef _OPENMP
  double MG_Start_Time = omp_get_wtime();
  #elif USE_MPI
  double MG_Start_Time = MPI_Wtime();
  #endif
  if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"MGSolve... ");}
  double _timeStartMGSolve = getTime();

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // calculate norm of f for convergence criteria...
  double norm_of_F     = 1.0;
  double norm_of_DinvF = 1.0;
  if(dtol>0){
    mul_vectors(all_grids->levels[onLevel],VECTOR_TEMP,1.0,F_id,VECTOR_DINV); // D^{-1}F
    norm_of_DinvF = norm(all_grids->levels[onLevel],VECTOR_TEMP);		// ||D^{-1}F||
  }
  if(rtol>0)norm_of_F = norm(all_grids->levels[onLevel],F_id);		// ||F||

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // make initial guess for e (=0) and setup the RHS
   zero_vector(all_grids->levels[onLevel],e_id);                  // ee = 0
  scale_vector(all_grids->levels[onLevel],R_id,1.0,F_id);         // R_id = F_id

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // now do v-cycles to calculate the correction...
  for(v=0;v<maxVCycles;v++){   
    int level = onLevel;
    all_grids->levels[level]->vcycles_from_this_level++;

    // do the v-cycle...
    MGVCycle(all_grids,e_id,R_id,a,b,level);

    // now calculate the norm of the residual...
    double _timeStart = getTime();
    if(all_grids->levels[level]->must_subtract_mean == 1){
      double average_value_of_e = mean(all_grids->levels[level],e_id);
      shift_vector(all_grids->levels[level],e_id,e_id,-average_value_of_e);
    }
    residual(all_grids->levels[level],VECTOR_TEMP,e_id,F_id,a,b);
    if(dtol>0)mul_vectors(all_grids->levels[level],VECTOR_TEMP,1.0,VECTOR_TEMP,VECTOR_DINV); //  Using ||D^{-1}(b-Ax)||_{inf} as convergence criteria...
    double norm_of_residual = norm(all_grids->levels[level],VECTOR_TEMP);
    double _timeNorm = getTime();
    all_grids->levels[level]->timers.Total += (double)(_timeNorm-_timeStart);
    if(all_grids->levels[level]->my_rank==0){
      double rel = 0.0;
      if(rtol>0)rel = norm_of_residual/norm_of_F;
           else rel = norm_of_residual/norm_of_DinvF;
      if(   v>0){fprintf(stdout,"\n           v-cycle=%2d  norm=%1.15e  rel=%1.15e  ",v+1,norm_of_residual,rel);}
            else{fprintf(stdout,             "v-cycle=%2d  norm=%1.15e  rel=%1.15e  ",v+1,norm_of_residual,rel);}
    }
    if(norm_of_residual/norm_of_F < rtol)break;
    if(norm_of_residual           < dtol)break;
  } // maxVCycles
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  all_grids->timers.MGSolve += (double)(getTime()-_timeStartMGSolve);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef _OPENMP
  if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"done (%f seconds)\n",omp_get_wtime()-MG_Start_Time);} // used to monitor variability in individual solve times
  #elif USE_MPI
  if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"done (%f seconds)\n",MPI_Wtime()-MG_Start_Time);} // used to monitor variability in individual solve times
  #else
  if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"done\n");}
  #endif
}


//------------------------------------------------------------------------------------------------------------------------------

void FMGSolve(mg_type *all_grids, int onLevel, int u_id, int F_id, double a, double b, double rtol){
 // FMGSolve will iterate on F-cycles then iterate on V-cycles
 // It does this by putting the system into residual correction
 // i.e. calculate the residual (becomes the RHS), solve for a correction, add the correction to the current solution, calculate a new residual, and repeat

 all_grids->MGSolves_performed++;
 if(!all_grids->levels[onLevel]->active)return;
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 #ifdef HPGMG_F_CYCLES
 int maxFCycles=HPGMG_F_CYCLES;
 #else
 int maxFCycles=20;
 #endif

 #ifdef HPGMG_V_CYCLES
 int maxVCycles=HPGMG_V_CYCLES;
 #else
 int maxVCycles= 0;
 #endif

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 int f,v;
 int level;
 int e_id = VECTOR_E;
 int R_id = VECTOR_R;
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 #ifdef _OPENMP
 double FMG_Start_Time = omp_get_wtime();
 #elif USE_MPI
 double FMG_Start_Time = MPI_Wtime();
 #endif
 if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"FMGSolve... ");}
 double _timeStartMGSolve = getTime();


 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // calculate norm of f, calculate R...
 double _LevelStart = getTime();
                         residual(all_grids->levels[onLevel],R_id,u_id,F_id,a,b);  // R_id = F-Au
 double norm_of_residual   = norm(all_grids->levels[onLevel],R_id);                // ||R||
 double norm_of_F =          norm(all_grids->levels[onLevel],F_id);                // ||F||
 all_grids->levels[onLevel]->timers.Total += (double)(getTime()-_LevelStart);


 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"            norm=%1.15e rel=%1.15e\n",norm_of_residual,norm_of_residual/norm_of_F);}
 //if(norm_of_residual > norm_of_F)norm_of_F=norm_of_residual;
 if(norm_of_residual/norm_of_F < rtol){maxFCycles=0;maxVCycles=0;} // initial guess was already converged

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // iterate on f-cycles...
 for(f=0;f<maxFCycles;f++){

   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   // restrict RHS to bottom (coarsest grids)
   for(level=onLevel;level<(all_grids->num_levels-1);level++){
     double _LevelStart = getTime();
     restriction(all_grids->levels[level+1],R_id,all_grids->levels[level],R_id,RESTRICT_CELL);
     all_grids->levels[level]->timers.Total += (double)(getTime()-_LevelStart);
   }


   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   // solve coarsest grid...
     double _timeBottomStart = getTime();
     level = all_grids->num_levels-1;
     if(level>onLevel)zero_vector(all_grids->levels[level],e_id);//else use whatever was the initial guess
     IterativeSolver(all_grids->levels[level],e_id,R_id,a,b,MG_DEFAULT_BOTTOM_NORM);  // -1 == exact solution
     all_grids->levels[level]->timers.Total += (double)(getTime()-_timeBottomStart);


   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   // now do the F-cycle proper...
   for(level=all_grids->num_levels-2;level>=onLevel;level--){
     // high-order interpolation
     _LevelStart = getTime();
     interpolation_fcycle(all_grids->levels[level],e_id,0.0,all_grids->levels[level+1],e_id);
     all_grids->levels[level]->timers.Total += (double)(getTime()-_LevelStart);

     // v-cycle
     all_grids->levels[level]->vcycles_from_this_level++;
     MGVCycle(all_grids,e_id,R_id,a,b,level);
   }

   // correct current solution and calculate residual (new RHS)...
   _LevelStart = getTime();
   add_vectors(all_grids->levels[onLevel],u_id,1.0,u_id,1.0,e_id); // u = u+e
   if(all_grids->levels[onLevel]->must_subtract_mean == 1){
     double average_value_of_u = mean(all_grids->levels[onLevel],u_id);
     shift_vector(all_grids->levels[onLevel],u_id,u_id,-average_value_of_u);
   }
   residual(all_grids->levels[onLevel],R_id,u_id,F_id,a,b);
   double norm_of_residual = norm(all_grids->levels[onLevel],R_id);
   all_grids->levels[onLevel]->timers.Total += (double)(getTime()-_LevelStart);

   // test convergence...
   if(all_grids->levels[onLevel]->my_rank==0){
     double rel = norm_of_residual/norm_of_F;
     fprintf(stdout,"            f-cycle=%2d  norm=%1.15e  rel=%1.15e\n",f+1,norm_of_residual,rel);
   }
   if(norm_of_residual/norm_of_F < rtol){maxVCycles=0;break;}

 } // F-cycle


 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 // iterate on v-cycles...
 for(v=0;v<maxVCycles;v++){

   int level = onLevel;
   all_grids->levels[level]->vcycles_from_this_level++;

   // do the v-cycle...
   zero_vector(all_grids->levels[onLevel],e_id);
   MGVCycle(all_grids,e_id,R_id,a,b,level);

   // correct current solution and calculate residual (new RHS)...
   _LevelStart = getTime();
   add_vectors(all_grids->levels[onLevel],u_id,1.0,u_id,1.0,e_id); // u = u+e
   if(all_grids->levels[onLevel]->must_subtract_mean == 1){
     double average_value_of_u = mean(all_grids->levels[onLevel],u_id);
     shift_vector(all_grids->levels[onLevel],u_id,u_id,-average_value_of_u);
   }
   residual(all_grids->levels[onLevel],R_id,u_id,F_id,a,b);
   double norm_of_residual = norm(all_grids->levels[onLevel],R_id);
   all_grids->levels[onLevel]->timers.Total += (double)(getTime()-_LevelStart);

   // test convergence...
   if(all_grids->levels[onLevel]->my_rank==0){
     double rel = norm_of_residual/norm_of_F;
     fprintf(stdout,"            v-cycle=%2d  norm=%1.15e  rel=%1.15e\n",v+1,norm_of_residual,rel);
   }
   if(norm_of_residual/norm_of_F < rtol)break;

 } // V-cycles

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 all_grids->timers.MGSolve += (double)(getTime()-_timeStartMGSolve);
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 #ifdef _OPENMP
 if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"            done (%f seconds)\n\n",omp_get_wtime()-FMG_Start_Time);} // used to monitor variability in individual solve times
 #elif USE_MPI
 if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"            done (%f seconds)\n\n",MPI_Wtime()-FMG_Start_Time);} // used to monitor variability in individual solve times
 #else
 if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"            done\n\n");}
 #endif

}


//------------------------------------------------------------------------------------------------------------------------------
void MGPCG(mg_type *all_grids, int onLevel, int x_id, int F_id, double a, double b, double dtol, double rtol){
  // Algorithm 9.1 in Iterative Methods for Sparse Linear Systems(Yousef Saad) using a MG V-Cycle as M^{-1}
  level_type * level = all_grids->levels[onLevel];
  if(!level->active)return;

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // CG with a MG preconditioner, every level needs 3 extra vectors (p, Ap, z)
  int l;
  for(l=0;l<all_grids->num_levels;l++){
    create_vectors(all_grids->levels[l],VECTORS_RESERVED+3);
  }

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // Test for Poisson with Periodic BCs
  for(l=0;l<all_grids->num_levels;l++){
    if(all_grids->levels[l]->must_subtract_mean==-1){
      all_grids->levels[l]->must_subtract_mean=0;
      int alpha_is_zero = (dot(all_grids->levels[l],VECTOR_ALPHA,VECTOR_ALPHA) == 0.0);
      if( (all_grids->levels[l]->boundary_condition.type==BC_PERIODIC) && ((a==0) || (alpha_is_zero)) )all_grids->levels[l]->must_subtract_mean = 1;
    }
  }

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int   r_id = VECTOR_R;
  int   p_id = VECTORS_RESERVED+0;
  int  Ap_id = VECTORS_RESERVED+1;
  int   z_id = VECTORS_RESERVED+2;

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef _OPENMP
  double MGPCG_Start_Time = omp_get_wtime();
  #elif USE_MPI
  double MGPCG_Start_Time = MPI_Wtime();
  #endif
  if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"MGPCG...  ");}
  double _timeStartMGSolve = getTime();
  all_grids->MGSolves_performed++;
  int jMax=20;
  int j=0;
  int CGFailed    = 0;
  int CGConverged = 0;

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  zero_vector(level,x_id);                                                      // x[] = 0
  residual(level,r_id,x_id,F_id,a,b);                                           // r[] = F_id[] - A(x_id)
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(level->must_subtract_mean == 1){
    double mean_of_r = mean(level,r_id);
    shift_vector(level,r_id,r_id,-mean_of_r);
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  double norm_of_r0 = norm(level,r_id);                                         // the norm of the initial residual...
  if(norm_of_r0 == 0.0){CGConverged=1;}                                         // entered CG with exact solution
  level->vcycles_from_this_level++;                                             //
  zero_vector(level,z_id);                                                      // z[] = 0
  MGVCycle(all_grids,z_id,r_id,a,b,onLevel);                                    // z[] = M^{-1}r[]
  scale_vector(level,p_id,1.0,z_id);                                            // p[] = z[]
  double r_dot_z = dot(level,r_id,z_id);                                        // r_dot_z = dot(r,z)
  while( (j<jMax) && (!CGFailed) && (!CGConverged) ){                           // while(not done){
    j++;level->Krylov_iterations++;                                             //
    apply_op(level,Ap_id,p_id,a,b);                                             //   Ap[] = A(p)
    double Ap_dot_p = dot(level,Ap_id,p_id);                                    //   Ap_dot_p = dot(Ap,p)
    if(Ap_dot_p == 0.0){CGFailed=1;break;}                                      //   pivot breakdown ???
    double alpha = r_dot_z / Ap_dot_p;                                          //   alpha = r_dot_z / Ap_dot_p
    if(isinf(alpha)){CGFailed=1;break;}                                         //   ???
    add_vectors(level,x_id,1.0,x_id, alpha,p_id );                              //   x_id[] = x_id[] + alpha*p[]
    add_vectors(level,r_id,1.0,r_id,-alpha,Ap_id);                              //   r[]    = r[]    - alpha*Ap[]   (intermediate residual?)
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(level->must_subtract_mean == 1){
      double mean_of_r = mean(level,r_id);
      shift_vector(level,r_id,r_id,-mean_of_r);
    }
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //double norm_of_r = norm(level,r_id);                                        //   norm of intermediate residual (delusional convergence)
    residual(level,VECTOR_TEMP,x_id,F_id,a,b);                                  //   true residual
    double norm_of_r = norm(level,VECTOR_TEMP);                                 //   norm of true residual (true convergence test)
    if(norm_of_r == 0.0){CGConverged=1;break;}                                  //
    if(level->my_rank==0){
      if(   j>1){fprintf(stdout,"\n          ");}
      if(rtol>0){fprintf(stdout,"iter=%3d  norm=%1.15e  rel=%1.15e  ",j,norm_of_r,norm_of_r/norm_of_r0    );}
    }
    if(norm_of_r/norm_of_r0 < rtol)break;                                       //   norm if true residual is small enough
    level->vcycles_from_this_level++;                                           //
    zero_vector(level,z_id);                                                    //   z[] = 0
    MGVCycle(all_grids,z_id,r_id,a,b,onLevel);                                  //   z[] = M^{-1}r[]
    double r_dot_z_new = dot(level,r_id,z_id);                                  //   r_dot_z_new = dot(r_{j+1},z_{j+1})
    if(r_dot_z_new == 0.0){CGFailed=1;break;}                                   //   Lanczos breakdown ???
    double beta = (r_dot_z_new/r_dot_z);                                        //   beta = (r_dot_z_new/r_dot_z)
    if(isinf(beta)){CGFailed=1;break;}                                          //   ???
    add_vectors(level,p_id,1.0,z_id,beta,p_id );                                //   p[] = z[] + beta*p[]
    r_dot_z = r_dot_z_new;                                                      //   r_dot_r = r_dot_r_new   (save old r_dot_r)
    // FIX... need to test for stalled convergence...
  }                                                                             // }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  all_grids->timers.MGSolve += (double)(getTime()-_timeStartMGSolve);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef _OPENMP
  if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"done (%f seconds)\n",omp_get_wtime()-MGPCG_Start_Time);} // used to monitor variability in individual solve times
  #elif USE_MPI
  if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"done (%f seconds)\n",MPI_Wtime()-MGPCG_Start_Time);} // used to monitor variability in individual solve times
  #else
  if(all_grids->levels[onLevel]->my_rank==0){fprintf(stdout,"done\n");}
  #endif
}
//------------------------------------------------------------------------------------------------------------------------------
