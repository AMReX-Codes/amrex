#ifndef BL_HPGMG_H
#define BL_HPGMG_H

#include <AMReX_MultiFab.H>

#ifdef USEHPGMG

#define STENCIL_MAX_SHAPES 3
#define VECTOR_ALPHA       5  // cell centered coefficient
#define BC_PERIODIC        0
#define VECTOR_E           1  // error used in residual correction FMG
#define VECTOR_F           3  // original right-hand side (Au=f), cell centered
#define VECTOR_U           4  // numerical solution
#define BC_DIRICHLET       1
#define VECTOR_BETA_I      6 // face centered coefficient (n.b. element 0 is the left face of the ghost zone element)
#define VECTOR_BETA_J      7 // face centered coefficient (n.b. element 0 is the back face of the ghost zone element)
#define VECTOR_BETA_K      8 // face centered coefficient (n.b. element 0 is the bottom face of the ghost zone element)
#define RESTRICT_CELL      0

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

typedef struct {
  int subtype;			// e.g. used to calculate normal to domain for BC's
  struct {int i, j, k;}dim;	// dimensions of the block to copy
  struct {int box, i, j, k, jStride, kStride;double * __restrict__ ptr;}read,write;
  // coordinates in the read grid to extract data,
  // coordinates in the write grid to insert data
  // if read/write.box<0, then use write/read.ptr, otherwise use boxes[box].vectors[id]
  // Thus, you can do grid->grid, grid->buf, buf->grid, or buf->buf
} __attribute__((aligned(64))) blockCopy_type;


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
    #ifdef BL_USE_MPI
    MPI_Request * __restrict__     requests;
    MPI_Status  * __restrict__       status;
    #endif
} communicator_type;

typedef struct {
  int                         global_box_id;	// used to inded into level->rank_of_box
  struct {int i, j, k;}low;			// global coordinates of the first (non-ghost) element of subdomain
  int                                   dim;	// dimension of this box's core (owned)
  int                                ghosts;	// ghost zone depth
  int                jStride,kStride,volume;	// useful for offsets
  int                            numVectors;	//
  double   ** __restrict__          vectors;	// vectors[c] = pointer to 3D array for vector c for one box
} box_type;

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
  #ifdef BL_USE_MPI
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
    uint64_t              smooth;
    uint64_t            apply_op;
    uint64_t            residual;
    uint64_t               blas1;
    uint64_t               blas3;
    uint64_t boundary_conditions;
    // Distributed Restriction
    uint64_t   restriction_total;
    uint64_t   restriction_pack;
    uint64_t   restriction_local;
    uint64_t   restriction_unpack;
    uint64_t   restriction_recv;
    uint64_t   restriction_send;
    uint64_t   restriction_wait;
    // Distributed interpolation
    uint64_t interpolation_total;
    uint64_t interpolation_pack;
    uint64_t interpolation_local;
    uint64_t interpolation_unpack;
    uint64_t interpolation_recv;
    uint64_t interpolation_send;
    uint64_t interpolation_wait;
    // Ghost Zone Exchanges...
    uint64_t     ghostZone_total;
    uint64_t     ghostZone_pack;
    uint64_t     ghostZone_local;
    uint64_t     ghostZone_unpack;
    uint64_t     ghostZone_recv;
    uint64_t     ghostZone_send;
    uint64_t     ghostZone_wait;
    // Collectives...
    uint64_t   collectives;
    uint64_t         Total;
  }cycles;
  int Krylov_iterations;        // total number of bottom solver iterations
  int CAKrylov_formations_of_G; // i.e. [G,g] = [P,R]^T[P,R,rt]
  int vcycles_from_this_level;  // number of vcycles performed that were initiated from this level
} level_type;

typedef struct {
  int num_ranks;	// total number of MPI ranks for MPI_COMM_WORLD
  int my_rank;		// my MPI rank for MPI_COMM_WORLD
  int       num_levels;	// depth of the v-cycle
  level_type ** levels;	// array of pointers to levels

  struct {
    uint64_t MGBuild; // total time spent building the coefficients...
    uint64_t MGSolve; // total time spent in MGSolve
  }cycles;
  int MGSolves_performed;
} mg_type;

extern "C" {
    void create_vectors(level_type* level, int numVectors);
    void append_block_to_list (blockCopy_type ** blocks, int *allocated_blocks,
                               int *num_blocks, int dim_i, int dim_j, int dim_k,
                               int read_box, double *read_ptr, int read_i,
                               int read_j, int read_k, int read_jStride,
                               int read_kStride, int read_scale, int write_box,
                               double *write_ptr, int write_i, int write_j,
                               int write_k, int write_jStride, int write_kStride,
                               int write_scale, int my_blockcopy_tile_i,
                               int my_blockcopy_tile_j, int my_blockcopy_tile_k,
                               int subtype);
    void build_exchange_ghosts (level_type * level, int shape);
    void build_boundary_conditions (level_type * level, int shape);
    void initialize_problem (level_type * level, double hLevel, double a, double b);
    double dot (level_type * level, int id_a, int id_b);
    double mean (level_type * level, int id_a);
    void shift_vector (level_type * level, int id_c, int id_a, double shift_a);
    void rebuild_operator (level_type * level, level_type * fromLevel, double a, double b);
    void MGBuild (mg_type * all_grids, level_type * fine_grid, double a, double b, int minCoarseGridDim, const MPI_Comm comm);
    void MGResetTimers (mg_type * all_grids);
    void zero_vector (level_type * level, int component_id);
    void MGSolve (mg_type *all_grids, int onLevel, int u_id, int F_id, double a, double b, double dtol, double rtol);
    void FMGSolve(mg_type *all_grids, int onLevel, int u_id, int F_id, double a, double b, double rtol);
    void MGPrintTiming (mg_type * all_grids, int fromLevel);
    double error (level_type * level, int id_a, int id_b);
    void destroy_level (level_type * level);
    int stencil_get_radius();
    double norm(level_type * level, int component_id);
    void restriction(level_type * level_c, int id_c, level_type *level_f, int id_f, int restrictionType);
    void richardson_error(mg_type *all_grids, int levelh, int u_id);
    void MGDestroy(mg_type *all_grids);
}

// If we want to use the multigrid solver from HPGMG then we must convert our
// MultiFabs to its own level data structures. This function essentially
// replaces the create_level() function in HPGMG.
void CreateHPGMGLevel(level_type* level,
                      const amrex::MultiFab& mf,
                      const int n_cell,
                      const int max_grid_size,
                      const int my_rank,
                      const int num_ranks,
                      const int domain_boundary_condition,
                      const int numVectors,
                      const double h0);

void SetupHPGMGCoefficients(const double a,
                            const double b,
                            const amrex::MultiFab& alpha,
                            const amrex::MultiFab& beta_cc,
                            level_type* level);

void ConvertToHPGMGLevel(const amrex::MultiFab& mf,
                         const int n_cell,
                         const int max_grid_size,
                         level_type* level,
                         const int component_id);

void ConvertFromHPGMGLevel(amrex::MultiFab& mf,
                           const level_type* level,
                           const int component_id);

#endif /* USEHPGMG */

#endif /* BL_HPGMG_H */
