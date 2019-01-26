//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#ifndef LEVEL_H
#define LEVEL_H
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
#ifdef USE_MPI
#include <mpi.h>
#endif
//------------------------------------------------------------------------------------------------------------------------------
// supported boundary conditions
#define BC_PERIODIC  0
#define BC_DIRICHLET 1
//------------------------------------------------------------------------------------------------------------------------------
// regiment communication by defining a series of stencil shapes...
#define STENCIL_SHAPE_BOX         0	// faces, edges, and corners
#define STENCIL_SHAPE_STAR        1	// just faces
#define STENCIL_SHAPE_NO_CORNERS  2	// faces and edges, but no corners
#define STENCIL_MAX_SHAPES        3
//------------------------------------------------------------------------------------------------------------------------------
// regiment threading around the 'block' or 'tile' concepts.  Define default tilings...
#ifndef BLOCKCOPY_TILE_I
#define BLOCKCOPY_TILE_I 10000
#else
#warning By overriding BLOCKCOPY_TILE_I, you are tiling in the unit stride.  I hope you know what you are doing.
#endif
#ifndef BLOCKCOPY_TILE_J
#define BLOCKCOPY_TILE_J 8
#endif
#ifndef BLOCKCOPY_TILE_K
#define BLOCKCOPY_TILE_K 8
#endif
//------------------------------------------------------------------------------------------------------------------------------
// FP data for a vector within a box is padded to ensure alignment
#ifndef BOX_ALIGN_JSTRIDE
#define BOX_ALIGN_JSTRIDE   4  // j-stride(unit stride dimension including ghosts and padding) is a multiple of BOX_ALIGN_JSTRIDE... useful for SIMD in j+/-1
#endif
#ifndef BOX_ALIGN_KSTRIDE
#define BOX_ALIGN_KSTRIDE   8  // k-stride is a multiple of BOX_ALIGN_KSTRIDE ... useful for SIMD in k+/-1
#endif
#ifndef BOX_ALIGN_VOLUME
#define BOX_ALIGN_VOLUME    8  // box volumes are a multiple of BOX_ALIGN_VOLUME ... useful for SIMD on different vectors
#endif
//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  int subtype;			// e.g. used to calculate normal to domain for BC's
  struct {int i, j, k;}dim;	// dimensions of the block to copy
  struct {int box, i, j, k, jStride, kStride;double * __restrict__ ptr;}read,write;
  // coordinates in the read grid to extract data, 
  // coordinates in the write grid to insert data
  // if read/write.box<0, then use write/read.ptr, otherwise use boxes[box].vectors[id]
  // Thus, you can do grid->grid, grid->buf, buf->grid, or buf->buf
} __attribute__((aligned(64))) blockCopy_type;


//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
    int                           num_recvs;	//   number of neighbors by type
    int                           num_sends;	//   number of neighbors by type
    int     * __restrict__       recv_ranks;	//   MPI rank of each neighbor...          recv_ranks[neighbor]
    int     * __restrict__       send_ranks;	//   MPI rank of each neighbor...          send_ranks[neighbor]
    int     * __restrict__       recv_sizes;	//   size of each MPI recv buffer...       recv_sizes[neighbor]
    int     * __restrict__       send_sizes;	//   size of each MPI send buffer...       send_sizes[neighbor]
    double ** __restrict__     recv_buffers;	//   MPI recv buffer for each neighbor...  recv_buffers[neighbor][ recv_sizes[neighbor] ]
    double ** __restrict__     send_buffers;	//   MPI send buffer for each neighbor...  send_buffers[neighbor][ send_sizes[neighbor] ]
    int                 allocated_blocks[3];	//   number of blocks allocated (not necessarily used) each list...
    int                       num_blocks[3];	//   number of blocks in each list...        num_blocks[pack,local,unpack]
    blockCopy_type *              blocks[3];	//   list of block copies...                     blocks[pack,local,unpack]
    #ifdef USE_MPI
    MPI_Request * __restrict__     requests;
    MPI_Status  * __restrict__       status;
    #endif
} communicator_type;


//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  int                         global_box_id;	// used to inded into level->rank_of_box
  struct {int i, j, k;}low;			// global coordinates of the first (non-ghost) element of subdomain
  int                                   dim;	// dimension of this box's core (owned)
  int                                ghosts;	// ghost zone depth
  int                jStride,kStride,volume;	// useful for offsets
  int                            numVectors;	//
  double   ** __restrict__          vectors;	// vectors[c] = pointer to 3D array for vector c for one box
} box_type;


//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  double h;					// grid spacing at this level
  int active;					// I am an active process (I have work to do on this or subsequent levels)
  int num_ranks;				// total number of MPI ranks
  int my_rank;					// my MPI rank
  int box_dim;					// dimension of each cubical box (not counting ghost zones)
  int box_ghosts;				// ghost zone depth for each box
  int box_jStride,box_kStride,box_volume;	// useful for offsets
  int numVectors;				// number of vectors stored in each box
  int tag;					// tag each level uniquely... FIX... replace with sub commuicator
  struct {int i, j, k;}boxes_in;		// total number of boxes in i,j,k across this level
  struct {int i, j, k;}dim;			// global dimensions at this level (NOTE: dim.i == boxes_in.i * box_dim)

  int * rank_of_box;				// 3D array containing rank of each box.  i-major ordering
  int    num_my_boxes;				//           number of boxes owned by this rank
  box_type * my_boxes;				// pointer to array of boxes owned by this rank

  // create flattened FP data... useful for CUDA/OpenMP4/OpenACC when you want to copy an entire vector to/from an accelerator
  double   ** __restrict__          vectors;	// vectors[v][box][k][j][i] = pointer to 5D array for vector v encompasing all boxes on this process... 
  double    * __restrict__     vectors_base;    // pointer used for malloc/free.  vectors[v] are shifted from this for alignment

  int       allocated_blocks;			//       number of blocks allocated by this rank (note, this represents a flattening of the box/cell hierarchy to facilitate threading)
  int          num_my_blocks;			//       number of blocks     owned by this rank (note, this represents a flattening of the box/cell hierarchy to facilitate threading)
  blockCopy_type * my_blocks;			// pointer to array of blocks owned by this rank (note, this represents a flattening of the box/cell hierarchy to facilitate threading)

  struct {
    int                type;			// BC_PERIODIC or BC_DIRICHLET
    int    allocated_blocks[STENCIL_MAX_SHAPES];// number of blocks allocated (not necessarily used) for boundary conditions on this level for [shape]
    int          num_blocks[STENCIL_MAX_SHAPES];// number of blocks used for boundary conditions on this level for [shape]
    blockCopy_type * blocks[STENCIL_MAX_SHAPES];// pointer to array of blocks used for boundary conditions on this level for [shape]
  } boundary_condition;				// boundary conditions on this level

  communicator_type exchange_ghosts[STENCIL_MAX_SHAPES];// mini program that performs a neighbor ghost zone exchange for [shape]
  communicator_type restriction[4];			// mini program that performs restriction and agglomeration for [0=cell centered, 1=i-face, 2=j-face, 3-k-face]
  communicator_type interpolation;			// mini program that performs interpolation and dissemination...
  #ifdef USE_MPI
  MPI_Comm MPI_COMM_ALLREDUCE;			// MPI sub communicator for just the ranks that have boxes on this level or any subsequent level... 
  #endif
  double dominant_eigenvalue_of_DinvA;		// estimate on the dominate eigenvalue of D^{-1}A
  int must_subtract_mean;			// e.g. Poisson with Periodic BC's
  double    * __restrict__ RedBlack_base;       // allocated pointer... will be aligned for the first non ghost zone element
  double    * __restrict__ RedBlack_FP;	        // Red/Black Mask (i.e. 0.0 or 1.0) for even/odd planes (2*kStride).  

  int num_threads;
  double    * __restrict__ fluxes;		// temporary array used to hold the flux values used by FV operators

  // statistics information...
  struct {
    double              smooth;
    double            apply_op;
    double            residual;
    double               blas1;
    double               blas3;
    double boundary_conditions;
    // Distributed Restriction
    double   restriction_total;
    double   restriction_pack;
    double   restriction_local;
    double   restriction_unpack;
    double   restriction_recv;
    double   restriction_send;
    double   restriction_wait;
    // Distributed interpolation
    double interpolation_total;
    double interpolation_pack;
    double interpolation_local;
    double interpolation_unpack;
    double interpolation_recv;
    double interpolation_send;
    double interpolation_wait;
    // Ghost Zone Exchanges...
    double     ghostZone_total;
    double     ghostZone_pack;
    double     ghostZone_local;
    double     ghostZone_unpack;
    double     ghostZone_recv;
    double     ghostZone_send;
    double     ghostZone_wait;
    // Collectives...
    double   collectives;
    double         Total;
  }timers;
  int Krylov_iterations;        // total number of bottom solver iterations
  int CAKrylov_formations_of_G; // i.e. [G,g] = [P,R]^T[P,R,rt]
  int vcycles_from_this_level;  // number of vcycles performed that were initiated from this level
} level_type;


//------------------------------------------------------------------------------------------------------------------------------
void create_level(level_type *level, int boxes_in_i, int box_dim, int box_ghosts, int numVectors, int domain_boundary_condition, int my_rank, int num_ranks, const MPI_Comm comm);
void destroy_level(level_type *level);
void create_vectors(level_type *level, int numVectors);
void reset_level_timers(level_type *level);
int qsortInt(const void *a, const void *b);
void append_block_to_list(blockCopy_type ** blocks, int *allocated_blocks, int *num_blocks,
                          int dim_i, int dim_j, int dim_k,
                          int  read_box, double*  read_ptr, int  read_i, int  read_j, int  read_k, int  read_jStride, int  read_kStride, int  read_scale,
                          int write_box, double* write_ptr, int write_i, int write_j, int write_k, int write_jStride, int write_kStride, int write_scale,
                          int my_blockcopy_tile_i, int my_blockcopy_tile_j, int my_blockcopy_tile_k,
                          int subtype
                         );
//------------------------------------------------------------------------------------------------------------------------------
#endif
