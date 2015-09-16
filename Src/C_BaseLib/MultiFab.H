
#ifndef BL_MULTIFAB_H
#define BL_MULTIFAB_H

#include <winstd.H>

#include <stdint.h>

#include <BLassert.H>
#include <FArrayBox.H>
#include <FabArray.H>



#ifdef USEHPGMG
// ---------------------------- BEGIN HPGMG STUFF ----------------------------
#define STENCIL_MAX_SHAPES        3
#define  VECTOR_ALPHA        5  // cell centered coefficient
#define BC_PERIODIC  0
#define  VECTOR_UTRUE        1  // exact solution used to generate f
#define  VECTOR_F            3  // original right-hand side (Au=f), cell centered
#define  VECTOR_U            4  // numerical solution
#define BC_DIRICHLET 1
#define  VECTOR_BETA_I       6 // face centered coefficient (n.b. element 0 is the left face of the ghost zone element)
#define  VECTOR_BETA_J       7 // face centered coefficient (n.b. element 0 is the back face of the ghost zone element)
#define  VECTOR_BETA_K       8 // face centered coefficient (n.b. element 0 is the bottom face of the ghost zone element)
#define RESTRICT_CELL   0

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
  double    * __restrict__ RedBlack_FP;	        // Red/Black Mask (i.e. 0.0 or 1.0) for even/odd planes (2*kStride).

  int num_threads;

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
    void initialize_valid_region (level_type * level);
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
    void FMGSolve(mg_type *all_grids, int onLevel, int u_id, int F_id, double a, double b, double dtol, double rtol);
    void MGPrintTiming (mg_type * all_grids, int fromLevel);
    double error (level_type * level, int id_a, int id_b);
    void destroy_level (level_type * level);
    int stencil_get_radius();
    double norm(level_type * level, int component_id);
    void restriction(level_type * level_c, int id_c, level_type *level_f, int id_f, int restrictionType);
    void richardson_error(mg_type *all_grids, int levelh, int u_id);
    void MGDestroy(mg_type *all_grids);
}
// ---------------------------- END HPGMG STUFF ----------------------------
#endif /* USEHPGMG */


//
// Forward declaration.
//
typedef FabArrayId MultiFabId;

class MultiFabCopyDescriptor;

class MultiFab;

namespace BoxLib
{
    void InterpAddBox (MultiFabCopyDescriptor& fabCopyDesc,
		       BoxList*                returnUnfilledBoxes,
		       Array<FillBoxId>&       returnedFillBoxIds,
		       const Box&              subbox,
		       MultiFabId              faid1,
		       MultiFabId              faid2,
		       Real                    t1,
		       Real                    t2,
		       Real                    t,
		       int                     src_comp,
		       int                     dest_comp,
		       int                     num_comp,
		       bool                    extrap);

    void InterpFillFab (MultiFabCopyDescriptor& fabCopyDesc,
			const Array<FillBoxId>& fillBoxIds,
			MultiFabId              faid1,
			MultiFabId              faid2,
			FArrayBox&              dest,
			Real                    t1,
			Real                    t2,
			Real                    t,
			int                     src_comp,
			int                     dest_comp,
			int                     num_comp,
			bool                    extrap);
}
//
// A Collection of FArrayBoxes
//
// The MultiFab class is publically derived from the
// FabArray<Real,FArrayBox> class.  It is a collection (stored as an array) of
// FArrayBoxes useful for storing floating point data on a domain defined by
// a union of rectangular regions embedded in a uniform index space.  The
// MultiFab class extends the function of the underlying FabArray class just
// as the FArrayBox class extends the funtion of BaseFab<Real>.  Additional
// member functions are defined for I/O and simple arithmetic operations on
// these aggregate objects.
//
// This class does NOT provide a copy constructor or assignment operator.
//
class MultiFab
    :
    public FabArray<FArrayBox>
{
public:
    //
    // Constructs an empty MultiFab.  Data can be defined at a later
    // time using the define member functions inherited
    // from FabArray.
    //
    MultiFab ();
    //
    // Constructs a MultiFab with a valid region defined by bxs and
    // a region of definition defined by the grow factor ngrow.
    // If mem_mode is defined to be Fab_allocate then FArrayBoxes are
    // allocated for each Box in the BoxArray.  The size of the
    // FArrayBox is given by the Box grown by ngrow and the number of
    // components is given by ncomp.  If mem_mode is defined to be
    // Fab_noallocate then no FArrayBoxes are allocated at this time
    // but can be defined later.
    //
    MultiFab (const BoxArray& bs,
              int             ncomp,
              int             ngrow,
              FabAlloc        mem_mode = Fab_allocate,
	      const IntVect&  nodal = IntVect::TheZeroVector());

    MultiFab (const BoxArray&            bs,
              int                        ncomp,
              int                        ngrow,
              const DistributionMapping& dm,
              FabAlloc                   mem_mode = Fab_allocate,
	      const IntVect&             nodal = IntVect::TheZeroVector());

    void operator= (const Real& r);
    //
    // Returns the minimum value contained in component comp of the
    // MultiFab.  The parameter nghost determines the number of
    // boundary cells to search for the minimum.  The default is to
    // search only the valid regions of the FArrayBoxes.
    //
    Real min (int comp,
              int nghost = 0,
	      bool local = false) const;
    //
    // Identical to the previous min() function, but confines its
    // search to intersection of Box b and the MultiFab.
    //
    Real min (const Box& b,
              int        comp,
              int        nghost = 0,
	      bool       local = false) const;
    //
    // Returns the maximum value contained in component comp of the
    // MultiFab.  The parameter nghost determines the number of
    // boundary cells to search for the maximum.  The default is to
    // search only the valid regions of the FArrayBoxes.
    //
    Real max (int comp,
              int nghost = 0,
	      bool local = false) const;
    //
    // Identical to the previous max() function, but confines its
    // search to intersection of Box b and the MultiFab.
    //
    Real max (const Box& b,
              int        comp,
              int        nghost = 0,
	      bool       local = false) const;
    //
    // Returns the maximum *absolute* value contained in 
    // component comp of the MultiFab.  No ghost cells are used.
    //
    Real norm0 (int comp = 0, bool local = false) const;
    //
    // Returns the maximum *absolute* value contained in 
    // component comp of the MultiFab in the intersection of the BoxArray
    // and the valid regions of the MultiFab.  No ghost cells are used.
    //
    Real norm0 (int comp, const BoxArray& ba, bool local = false) const;
    //
    // Returns the maximum *absolute* values contained in 
    // each component of "comps" of the MultiFab.  No ghost cells are used.
    //
    Array<Real> norm0 (const Array<int>& comps, bool local = false) const;
    //
    // Returns the L1 norm of component "comp" over the MultiFab.
    // ngrow ghost cells are used.
    //
    Real norm1 (int comp = 0, int ngrow = 0, bool local = false) const;
    //
    // Returns the L1 norm of each component of "comps" over the MultiFab.
    // ngrow ghost cells are used.
    //
    Array<Real> norm1 (const Array<int>& comps, int ngrow = 0, bool local = false) const;
    //
    // Returns the L2 norm of component "comp" over the MultiFab.
    // No ghost cells are used.
    //
    Real norm2 (int comp = 0) const;
    //
    // Returns the L2 norm of each component of "comps" over the MultiFab.
    // No ghost cells are used.
    //
    Array<Real> norm2 (const Array<int>& comps) const;
    //
    // Adds the scalar value val to the value of each cell in the
    // specified subregion of the MultiFab.  The subregion consists
    // of the num_comp components starting at component comp.
    // The value of nghost specifies the number of cells in the
    // boundary region of each FArrayBox in the subregion that should
    // be modified.  
    //
    void plus (Real val,
               int  comp,
               int  num_comp,
               int  nghost = 0);
    //
    // Identical to the previous version of plus(), with the
    // restriction that the subregion is further constrained to
    // the intersection with Box region.
    //
    void plus (Real       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);
    //
    // Adds the scalar value val to the value of each cell in the
    // valid region of each component of the MultiFab.  The value
    // of nghost specifies the number of cells in the boundary
    // region that should be modified.
    //
    void plus (Real val,
               int  nghost);
    //
    // Adds the scalar value val to the value of each cell in the
    // valid region of each component of the MultiFab, that also
    // intersects the Box region.  The value of nghost specifies the
    // number of cells in the boundary region of each FArrayBox in
    // the subregion that should be modified.
    //
    void plus (Real       val,
               const Box& region,
               int        nghost);
    //
    // Scales the value of each cell in the specified subregion of the
    // MultiFab by the scalar val (a[i] <- a[i]*val). The subregion
    // consists of the num_comp components starting at component comp.
    // The value of nghost specifies the number of cells in the
    // boundary region of each FArrayBox in the subregion that should
    // be modified.  
    //
    void mult (Real val,
               int  comp,
               int  num_comp,
               int  nghost = 0);
    //
    // Identical to the previous version of mult(), with the
    // restriction that the subregion is further constrained to the
    // intersection with Box region.  The value of nghost specifies the
    // number of cells in the boundary region of each FArrayBox in
    // the subregion that should be modified.
    //
    void mult (Real       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);
    //
    // Scales the value of each cell in the valid region of each
    // component of the MultiFab by the scalar val (a[i] <- a[i]*val).
    // The value of nghost specifies the number of cells in the
    // boundary region that should be modified.
    //
    void mult (Real val,
               int  nghost = 0);
    //
    // Scales the value of each cell in the valid region of each
    // component of the MultiFab by the scalar val (a[i] <- a[i]*val),
    // that also intersects the Box region.  The value of nghost
    // specifies the number of cells in the boundary region of each
    // FArrayBox in the subregion that should be modified.
    //
    void mult (Real       val,
               const Box& region,
               int        nghost = 0);
    //
    // Replaces the value of each cell in the specified subregion of
    // the MultiFab with its reciprocal multiplied by the value of
    // numerator. The subregion consists of the num_comp components
    // starting at component comp.  The value of nghost specifies the
    // number of cells in the boundary region of each FArrayBox in the
    // subregion that should be modified.
    //
    void invert (Real numerator,
                 int  comp,
                 int  num_comp,
                 int  nghost = 0);
    //
    // Identical to the previous version of invert(), with the
    // restriction that the subregion is further constrained to the
    // intersection with Box region.  The value of nghost specifies the
    // number of cells in the boundary region of each FArrayBox in the
    // subregion that should be modified.
    //
    void invert (Real       numerator,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    //
    // Replaces the value of each cell in the specified subregion of
    // the MultiFab with its reciprocal multiplied by the value of
    // numerator.  The value of nghost specifies the number of cells
    // in the boundary region that should be modified.
    //
    void invert (Real numerator,
                 int  nghost);
    //
    // Replaces the value of each cell in the specified subregion of
    // the MultiFab, that also intersects the Box region, with its
    // reciprocal multiplied by the value of numerator.  The value
    // of nghost specifies the number of cells in the boundary region
    // of each FArrayBox in the subregion that should be modified.
    //
    void invert (Real       numerator,
                 const Box& region,
                 int        nghost);
    //
    // Negates the value of each cell in the specified subregion of
    // the MultiFab.  The subregion consists of the num_comp
    // components starting at component comp.  The value of nghost
    // specifies the number of cells in the boundary region of each
    // FArrayBox in the subregion that should be modified.  
    //
    void negate (int comp,
                 int num_comp,
                 int nghost = 0);
    //
    // Identical to the previous version of negate(), with the
    // restriction that the subregion is further constrained to
    // the intersection with Box region.
    //
    void negate (const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    //
    // Negates the value of each cell in the valid region of
    // the MultiFab.  The value of nghost specifies the number of
    // cells in the boundary region that should be modified.  
    //
    void negate (int nghost = 0);
    //
    // Negates the value of each cell in the valid region of
    // the MultiFab that also intersects the Box region.  The value
    // of nghost specifies the number of cells in the boundary region
    // that should be modified.  
    //
    void negate (const Box& region,
                 int        nghost = 0);

    IntVect minIndex (int comp,
                      int nghost = 0) const;

    IntVect maxIndex (int comp,
                      int nghost = 0) const;
    //
    // This function adds the values of the cells in mf to the corresponding
    // cells of this MultiFab.  mf is required to have the same BoxArray or
    // "valid region" as this MultiFab.  The addition is done only to num_comp
    // components, starting with component number strt_comp.  The parameter
    // nghost specifies the number of boundary cells that will be modified.
    // If nghost == 0, only the valid region of each FArrayBox will be
    // modified.
    //
    void plus (const MultiFab& mf,
               int             strt_comp,
               int             num_comp,
               int             nghost);
    //
    //
    // This function subtracts the values of the cells in mf from the
    // corresponding cells of this MultiFab.  mf is required to have the
    // same BoxArray or "valid region" as this MultiFab.  The subtraction is
    // done only to num_comp components, starting with component number
    // strt_comp.  The parameter nghost specifies the number of boundary
    // cells that will be modified.  If nghost == 0, only the valid region of
    // each FArrayBox will be modified.
    //
    void minus (const MultiFab& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost);
    //
    //
    // This function divides the values of the cells in mf from the
    // corresponding cells of this MultiFab.  mf is required to have the
    // same BoxArray or "valid region" as this MultiFab.  The division is
    // done only to num_comp components, starting with component number
    // strt_comp.  The parameter nghost specifies the number of boundary
    // cells that will be modified.  If nghost == 0, only the valid region of
    // each FArrayBox will be modified.  Note, nothing is done to protect 
    // against divide by zero.
    //
    void divide (const MultiFab& mf,
		 int             strt_comp,
		 int             num_comp,
		 int             nghost);
    //
    // Add src to dst including nghost ghost cells.
    // The two MultiFabs MUST have the same underlying BoxArray.
    //
    static void Add (MultiFab&       dst,
		     const MultiFab& src,
		     int             srccomp,
		     int             dstcomp,
		     int             numcomp,
		     int             nghost);
    //
    // Copy from src to dst including nghost ghost cells.
    // The two MultiFabs MUST have the same underlying BoxArray.
    // The copy is local.  The parallel copy function is in FabArray.H
    //
    static void Copy (MultiFab&       dst,
                      const MultiFab& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);
    //
    // Subtract src from dst including nghost ghost cells.
    // The two MultiFabs MUST have the same underlying BoxArray.
    //
    static void Subtract (MultiFab&       dst,
			  const MultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);
    //
    // Multiply dst by src including nghost ghost cells.
    // The two MultiFabs MUST have the same underlying BoxArray.
    //
    static void Multiply (MultiFab&       dst,
			  const MultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);
    //
    // Divide dst by src including nghost ghost cells.
    // The two MultiFabs MUST have the same underlying BoxArray.
    //
    static void Divide (MultiFab&       dst,
			const MultiFab& src,
			int             srccomp,
			int             dstcomp,
			int             numcomp,
			int             nghost);

    void define (const BoxArray& bxs,
                 int             nvar,
                 int             ngrow,
                 FabAlloc        mem_mode,
		 const IntVect&  nodal = IntVect::TheZeroVector());

    void define (const BoxArray&            bxs,
                 int                        nvar,
                 int                        ngrow,
		 const DistributionMapping& dm,
                 FabAlloc                   mem_mode,
		 const IntVect&             nodal = IntVect::TheZeroVector());
    //
    // Are there any NaNs in the MF?
    // This may return false, even if the MF contains NaNs, if the machine
    // doesn't support the appropriate NaN testing functions.
    //
    // This version checks all components and grow cells.
    //
    bool contains_nan (bool local=false) const;

    bool contains_nan (int scomp, int ncomp, int ngrow = 0, bool local=false) const;
    //
    // Are there any Infs in the MF?
    // This may return false, even if the MF contains Infs, if the machine
    // doesn't support the appropriate Inf testing functions.
    //
    // This version checks all components and grow cells.
    //
    bool contains_inf (bool local=false) const;

    bool contains_inf (int scomp, int ncomp, int ngrow = 0, bool local=false) const;

    const FArrayBox& operator[] (int K) const;

    FArrayBox& operator[] (int K);

    const FArrayBox& operator[] (const MFIter& mfi) const { 
	return this->FabArray<FArrayBox>::get(mfi); } 

    const FArrayBox& get (const MFIter& mfi) const { return operator[](mfi); }

    FArrayBox& operator[] (const MFIter& mfi) { 
	return this->FabArray<FArrayBox>::get(mfi); } 

    FArrayBox& get (const MFIter& mfi) { return operator[](mfi); }
    //
    // Set ghost cells that are covered by valid cells to the value of valid cells.
    //
    void FillBoundary (bool local = false, bool cross = false);

    void FillBoundary (int scomp, int ncomp, bool local = false, bool cross = false);
    //
    // Sum ghost cells that are covered by valid cells into the valid cells.
    //
    void SumBoundary ();

    void SumBoundary (int scomp, int ncomp);

    static void Initialize ();
    static void Finalize ();
    //
    // Send a MultiFab from Compute MPI group to Sidecar MPI group.
    // The pointer to the MultiFab on the sidecars will be null on the compute
    // processes, and vice versa. So be mindful of which MPI processes you're
    // on when you access these pointers.
    //
    // TODO: allow it to send only selected components
    //
#ifdef BL_USE_MPI
    static void SendMultiFabToSidecars(MultiFab *mf);
#endif

#ifdef USEHPGMG
    // If we want to use the multigrid solver from HPGMG then we must convert our
    // MultiFabs to its own level data structures. This function essentially
    // replaces the create_level() function in HPGMG.
    static void CreateHPGMGLevel(level_type* level,
                                 const MultiFab& mf,
                                 const int n_cell,
                                 const int max_grid_size,
                                 const int my_rank,
                                 const int num_ranks,
                                 const int domain_boundary_condition,
                                 const int numVectors,
                                 const double h0);

    static void SetupHPGMGCoefficients(const double a,
                                       const double b,
                                       const MultiFab& alpha,
                                       const MultiFab& beta_cc,
                                       level_type* level);

    static void ConvertToHPGMGLevel(const MultiFab& mf,
                                    const int n_cell,
                                    const int max_grid_size,
                                    level_type* level,
                                    const int component_id);

    static void ConvertFromHPGMGLevel(MultiFab& mf,
                                      const level_type* level,
                                      const int component_id);
#endif

private:
    //
    // These are disabled.
    //
    MultiFab (const MultiFab& rhs);
    MultiFab& operator= (const MultiFab& rhs);

    static bool check_for_nan;
    static bool check_for_inf;

    //
    // Some useful typedefs.
    //
    typedef FabArrayBase::CopyComTagsContainer CopyComTagsContainer;
    typedef FabArrayBase::MapOfCopyComTagContainers MapOfCopyComTagContainers;
};

class MultiFabCopyDescriptor
    :
    public FabArrayCopyDescriptor<FArrayBox>
{
  public:

    MultiFabCopyDescriptor ();

    ~MultiFabCopyDescriptor ();

    MultiFabId RegisterMultiFab (MultiFab* mf) { return RegisterFabArray(mf); }

  private:
    //
    // These are disallowed.
    //
    MultiFabCopyDescriptor (const MultiFabCopyDescriptor&);
    MultiFabCopyDescriptor& operator= (const MultiFabCopyDescriptor&);
};

#endif /*BL_MULTIFAB_H*/
