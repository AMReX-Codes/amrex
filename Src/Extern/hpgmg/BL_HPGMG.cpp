#include <BL_HPGMG.H>

using namespace amrex;

// If we want to use the multigrid solver from HPGMG then we must convert our
// MultiFabs to HPGMG's level data structures. This function essentially
// replaces the create_level() function in HPGMG.
#ifdef USEHPGMG
void CreateHPGMGLevel (level_type* level,
                       const MultiFab& mf,
                       const int n_cell,
                       const int max_grid_size,
                       const int my_rank,
                       const int num_ranks,
                       const int domain_boundary_condition,
                       const int numVectors,
                       const double h0)
{
    int box;
    const int boxes_in_i = n_cell / max_grid_size;
    int TotalBoxes = boxes_in_i * boxes_in_i * boxes_in_i;

    // HPGMG requires perfect cubes for all boxes
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        if (!bx.isSquare()) {
             amrex::Error("All boxes must be square in HPGMG");
        }
    }

    // HPGMG also requires all boxes to be the same size, so we iterate over
    // all boxes and make sure they're the same.
    for (MFIter mfi1(mf); mfi1.isValid(); ++mfi1)
    {
        const Box& bx1 = mfi1.validbox();
        for (MFIter mfi2(mf); mfi2.isValid(); ++mfi2)
        {
            const Box& bx2 = mfi2.validbox();
            if (!(bx1.sameSize(bx2)))
            {
                amrex::Error("All boxes must be identical in HPGMG!");
            }
        }
    }

    // All the boxes have identical size and shape, so we just pick one of them
    // as a representative to fill in all the level data for HPGMG.
    MFIter mfi(mf);
    while (!mfi.isValid()) ++mfi;

    const Box& bx = mfi.validbox();
    const int box_dim = bx.length(0); /* Since we've already checked that all boxes are the same size, we can just use the size from one of them here. */

    if (TotalBoxes / num_ranks == 0)
      amrex::Error("Must have at least one box per MPI task when using HPGMG");

    if (ParallelDescriptor::IOProcessor())
    {
      std::cout << std::endl << "attempting to create a " << box_dim*boxes_in_i << "^3 level from " << TotalBoxes << " x " << box_dim << "^3 boxes distributed among " << num_ranks << " tasks..." << std::endl;
      if (domain_boundary_condition==BC_DIRICHLET)
      {
        std::cout << "boundary condition = BC_DIRICHLET" << std::endl;
      }
      else if (domain_boundary_condition==BC_PERIODIC)
      {
        std::cout << "boundary condition = BC_PERIODIC" << std::endl;
      }
      else
      {
        amrex::Error("Unknown boundary condition supplied");
      }
    }

    int omp_threads = 1;

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
      {
        omp_threads = omp_get_num_threads ();
      }
    }
#endif

    int box_ghosts = stencil_get_radius();

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
    level->h                = h0;
    level->fluxes           = NULL;

    // allocate 3D array of integers to hold the MPI rank of the corresponding box and initialize to -1 (unassigned)
    level->rank_of_box = (int*)malloc(level->boxes_in.i*level->boxes_in.j*level->boxes_in.k*sizeof(int));
    if(level->rank_of_box==NULL)
        amrex::Error("malloc of level->rank_of_box failed");
    for(box=0;box<level->boxes_in.i*level->boxes_in.j*level->boxes_in.k;box++){level->rank_of_box[box]=-1;}  // -1 denotes that there is no actual box assigned to this region


    // Now convert our rank distribution of boxes to HPGMG's rank_of_box array.
    // This is convoluted because HPGMG first assigns boxes to ranks, and then
    // lexicographically assigns the coordinates of each box. This
    // lexicographical ordering of box coordinates is *required* in order for
    // the MPI communication patterns in HPGMG to function correctly, via the
    // global_box_id variable. In other words, HPGMG anticipates the geometric
    // relationship between boxes based on their respective values of
    // global_box_id, and routes MPI traffic accordingly. However, in BoxLib
    // the box ranks and indices are not necessarily in this order, so we have
    // to "fake" the box ordering in HPGMG here (even though the coordinates
    // aren't actually assigned until we call create_vectors()) in order to
    // match the box ranks between BoxLib and HPGMG. This whole method is dumb
    // and deserves a better solution, but I don't know a better way to do it.

    int num_local_boxes = 0;
    int i,j,k;
    for(k=0;k<level->boxes_in.k;k++){
    for(j=0;j<level->boxes_in.j;j++){
    for(i=0;i<level->boxes_in.i;i++){
      int jStride = level->boxes_in.i;
      int kStride = level->boxes_in.i*level->boxes_in.j;
      int b=i + j*jStride + k*kStride;

      // These will be the coordinates of a box in HPGMG. These are also the
      // coordinates of a box already created in BoxLib. Now we iterate through
      // every rank's local boxes until we find the matching one, and assign
      // the rank of the HPGMG box to the same rank in BoxLib.

      const int low_i      = i*level->box_dim;
      const int low_j      = j*level->box_dim;
      const int low_k      = k*level->box_dim;

      bool found = false;
      for (MFIter mfi(mf); mfi.isValid(); ++mfi)
      {
        const Box &bx = mfi.validbox();
        const int *loVect = bx.loVect();

        // Found the matching box!
        if ((low_i == loVect[0]) &&
            (low_j == loVect[1]) &&
            (low_k == loVect[2]))
        {
            found = true;
            num_local_boxes++;
            break;
        }
      }
      if (found)
      {
        level->rank_of_box[b] = my_rank;
      }
    }}}

    // Now tell all the ranks what each other's box ranks are.
    const int tot_num_boxes = level->boxes_in.i * level->boxes_in.j * level->boxes_in.k;
    int all_box_ranks[tot_num_boxes];
    std::fill_n(all_box_ranks, tot_num_boxes, 1);
    MPI_Allreduce(level->rank_of_box, all_box_ranks, tot_num_boxes, MPI_INT, MPI_PROD, ParallelDescriptor::Communicator());
    for (unsigned int i = 0; i < tot_num_boxes; ++i)
    {
        level->rank_of_box[i] = std::abs(all_box_ranks[i]);
    }

    std::vector<int> box_ranks(level->rank_of_box, level->rank_of_box + tot_num_boxes);

    // calculate how many boxes I own...
    level->num_my_boxes=0;
    for(box=0;box<level->boxes_in.i*level->boxes_in.j*level->boxes_in.k;box++){if(level->rank_of_box[box]==level->my_rank)level->num_my_boxes++;}
    level->my_boxes = (box_type*)malloc(level->num_my_boxes*sizeof(box_type));
    if((level->num_my_boxes>0)&&(level->my_boxes==NULL))
        amrex::Error("malloc failed - create_level/level->my_boxes");

    // allocate flattened vector FP data and create pointers...
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Allocating vectors... ";
    create_vectors (level, numVectors);
    if (ParallelDescriptor::IOProcessor())
        std::cout << "done." << std::endl;

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
    level->RedBlack_FP = NULL;
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
      }}
    }

    int shape;
    // create mini program for each stencil shape to perform a ghost zone exchange...
    for(shape=0;shape<STENCIL_MAX_SHAPES;shape++)build_exchange_ghosts(    level,shape);
    // create mini program for each stencil shape to perform a boundary condition...
    for(shape=0;shape<STENCIL_MAX_SHAPES;shape++)build_boundary_conditions(level,shape);


    // duplicate the parent communicator to be the communicator for each level
    #ifdef BL_USE_MPI
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Duplicating MPI communicator... ";
    double time_start = MPI_Wtime();
    MPI_Comm_dup(ParallelDescriptor::Communicator(),&level->MPI_COMM_ALLREDUCE);
    double time_end = MPI_Wtime();
    double time_in_comm_dup = 0;
    double time_in_comm_dup_send = time_end-time_start;
    MPI_Allreduce(&time_in_comm_dup_send,&time_in_comm_dup,1,MPI_DOUBLE,MPI_MAX,ParallelDescriptor::Communicator());
    if (ParallelDescriptor::IOProcessor())
      std::cout << "done (" << time_in_comm_dup << " seconds)" << std::endl;
    #endif /* BL_USE_MPI */

    // report on potential load imbalance
    int BoxesPerProcess = level->num_my_boxes;
    #ifdef BL_USE_MPI
    int BoxesPerProcessSend = level->num_my_boxes;
    MPI_Allreduce(&BoxesPerProcessSend,&BoxesPerProcess,1,MPI_INT,MPI_MAX,ParallelDescriptor::Communicator());
    #endif /* BL_USE_MPI */
    if (ParallelDescriptor::IOProcessor())
      std::cout << "Calculating boxes per process... target=" << (double)TotalBoxes/(double)num_ranks << ", max=" << BoxesPerProcess << std::endl;
}


void SetupHPGMGCoefficients(const double a,
                            const double b,
                            const MultiFab& alpha,
                            const MultiFab& beta_cc,
                            level_type* level)
{

    // First set the alphas (cell-centered).
    bool found = false;
    for (MFIter mfi(alpha); mfi.isValid(); ++mfi) {

      const Box &bx = mfi.validbox();

      const int *loVect = bx.loVect();
      unsigned int box;
      for (box = 0; box < level->num_my_boxes; ++box)
      {
        if ((level->my_boxes[box].low.i == loVect[0]) &&
            (level->my_boxes[box].low.j == loVect[1]) &&
            (level->my_boxes[box].low.k == loVect[2]))
        {
          found = true;
          break;
        }
      }
      if (!found)
      {
        amrex::Error("Could not find matching boxes between HPGMG and BoxLib");
      }

      const Box &fabbox = mfi.fabbox();
      const double *alpha_data_ptr = alpha[mfi].dataPtr();
      int i,j,k;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int  ghosts = level->my_boxes[box].ghosts;
      const int   dim_i = level->my_boxes[box].dim;
      const int   dim_j = level->my_boxes[box].dim;
      const int   dim_k = level->my_boxes[box].dim;

      const int BL_jStride = fabbox.length(0);
      const int BL_kStride = fabbox.length(0) * fabbox.length(1);
      const int BoxLib_ghosts = alpha.nGrow();
      #ifdef _OPENMP
      #pragma omp parallel for private(k,j,i) collapse(3)
      #endif
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){
        int ijk_HPGMG = (i+ghosts) + (j+ghosts)*jStride + (k+ghosts)*kStride;
        const int ijk_BoxLib = (i+BoxLib_ghosts) + (j+BoxLib_ghosts)*BL_jStride + (k+BoxLib_ghosts)*BL_kStride;
        level->my_boxes[box].vectors[VECTOR_ALPHA][ijk_HPGMG] = alpha_data_ptr[ijk_BoxLib];
      }}}
    }


    // Now convert the cell-centered beta to faces.
    found = false;
    for (MFIter mfi(beta_cc); mfi.isValid(); ++mfi) {

      const Box &bx = mfi.validbox();

      const int *loVect = bx.loVect();
      unsigned int box;
      for (box = 0; box < level->num_my_boxes; ++box)
      {
        if ((level->my_boxes[box].low.i == loVect[0]) &&
            (level->my_boxes[box].low.j == loVect[1]) &&
            (level->my_boxes[box].low.k == loVect[2]))
        {
          found = true;
          break;
        }
      }
      if (!found)
      {
        amrex::Error("Could not find matching boxes between HPGMG and BoxLib");
      }

      const Box &fabbox = mfi.fabbox();

      const double *beta_data_ptr = beta_cc[mfi].dataPtr();
      int i,j,k;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int  ghosts = level->my_boxes[box].ghosts;
      const int   dim_i = level->my_boxes[box].dim;
      const int   dim_j = level->my_boxes[box].dim;
      const int   dim_k = level->my_boxes[box].dim;
      const int BL_jStride = fabbox.length(0);
      const int BL_kStride = fabbox.length(0) * fabbox.length(1);
      const int BoxLib_ghosts = beta_cc.nGrow();

      #ifdef _OPENMP
      #pragma omp parallel for private(k,j,i) collapse(3)
      #endif
      for(k=0;k<=dim_k;k++){ // include high face
      for(j=0;j<=dim_j;j++){ // include high face
      for(i=0;i<=dim_i;i++){ // include high face
        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        int ijk_HPGMG = (i+ghosts) + (j+ghosts)*jStride + (k+ghosts)*kStride;
        const int ijk_BoxLib   = (i  +BoxLib_ghosts) + (j  +BoxLib_ghosts)*BL_jStride + (k  +BoxLib_ghosts)*BL_kStride;
        const int im1jk_BoxLib = (i-1+BoxLib_ghosts) + (j  +BoxLib_ghosts)*BL_jStride + (k  +BoxLib_ghosts)*BL_kStride;
        const int ijm1k_BoxLib = (i  +BoxLib_ghosts) + (j-1+BoxLib_ghosts)*BL_jStride + (k  +BoxLib_ghosts)*BL_kStride;
        const int ijkm1_BoxLib = (i  +BoxLib_ghosts) + (j  +BoxLib_ghosts)*BL_jStride + (k-1+BoxLib_ghosts)*BL_kStride;
        level->my_boxes[box].vectors[VECTOR_BETA_I][ijk_HPGMG] = 0.5 * (beta_data_ptr[ijk_BoxLib] + beta_data_ptr[im1jk_BoxLib]);
        level->my_boxes[box].vectors[VECTOR_BETA_J][ijk_HPGMG] = 0.5 * (beta_data_ptr[ijk_BoxLib] + beta_data_ptr[ijm1k_BoxLib]);
        level->my_boxes[box].vectors[VECTOR_BETA_K][ijk_HPGMG] = 0.5 * (beta_data_ptr[ijk_BoxLib] + beta_data_ptr[ijkm1_BoxLib]);
      }}}
    }
}


void ConvertToHPGMGLevel (const MultiFab& mf,
                     const int n_cell,
                     const int max_grid_size,
                     level_type* level,
                     const int component_id)
{
    bool found = false;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

      const Box &bx = mfi.validbox();

      // The local box indices are ordered differently in HPGMG and BoxLib. So
      // as a simple (but SLOW) hack, we just find the boxes with matching
      // lower indices.
      // TODO: make this box matching less hacky
      const int *loVect = bx.loVect();
      unsigned int box;
      for (box = 0; box < level->num_my_boxes; ++box)
      {
        if ((level->my_boxes[box].low.i == loVect[0]) &&
            (level->my_boxes[box].low.j == loVect[1]) &&
            (level->my_boxes[box].low.k == loVect[2]))
        {
          found = true;
          break;
        }
      }

      if (!found)
      {
        amrex::Error("Could not find matching boxes between HPGMG and BoxLib");
      }

      const Box &fabbox = mfi.fabbox();
      const int BL_jStride = fabbox.length(0);
      const int BL_kStride = fabbox.length(0) * fabbox.length(1);

      const double *fab_data = mf[mfi].dataPtr();
      int i,j,k;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int  ghosts = level->my_boxes[box].ghosts;
      const int   dim_i = level->my_boxes[box].dim;
      const int   dim_j = level->my_boxes[box].dim;
      const int   dim_k = level->my_boxes[box].dim;
      const int BoxLib_ghosts = mf.nGrow();

      #ifdef _OPENMP
      #pragma omp parallel for private(k,j,i) collapse(3)
      #endif
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){

	// The HPGMG strides are padded to align memory and encourage
	// SIMD-ization, so they are different than the BoxLib strides.

        const int ijk_HPGMG = (i+ghosts) + (j+ghosts)*jStride + (k+ghosts)*kStride;
        const int ijk_BoxLib = (i+BoxLib_ghosts) + (j+BoxLib_ghosts)*BL_jStride + (k+BoxLib_ghosts)*BL_kStride;

        level->my_boxes[box].vectors[component_id][ijk_HPGMG] = fab_data[ijk_BoxLib];

      }}}

    }
}

void ConvertFromHPGMGLevel(MultiFab& mf,
                           const level_type* level,
                           const int component_id)
{
  for (MFIter mfi(mf); mfi.isValid(); ++mfi)
  {
      const Box &bx = mfi.validbox();
      double *fab_data = mf[mfi].dataPtr();

      // First find the HPGMG box corresponding to this BoxLib box.
      const int *loVect = bx.loVect();
      int box;
      for (box = 0; box < level->num_my_boxes; ++box)
      {
        if ((level->my_boxes[box].low.i == loVect[0]) &&
            (level->my_boxes[box].low.j == loVect[1]) &&
            (level->my_boxes[box].low.k == loVect[2]))
          break;
      }

      const Box &fabbox = mfi.fabbox();

      // Found the matching boxes, now fill the data.
      const int dim_i = level->my_boxes[box].dim;
      const int dim_j = level->my_boxes[box].dim;
      const int dim_k = level->my_boxes[box].dim;
      const int ghosts = level->my_boxes[box].ghosts;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int BoxLib_ghosts = mf.nGrow();

      int i, j, k;
      #ifdef _OPENMP
      #pragma omp parallel for private(k,j,i) collapse(3)
      #endif
      for(k=0;k<dim_k;k++){
      for(j=0;j<dim_j;j++){
      for(i=0;i<dim_i;i++){

        const int ijk_HPGMG = (i+ghosts) + (j+ghosts)*jStride + (k+ghosts)*kStride;

        // WARNING: this indexing stride works for FABs *ONLY* if we have ONE
        // component in the FAB. If we have more than one we have to stride
        // over the components in the outermost loop (outside of k).
        const int BL_jStride = fabbox.length(0);
        const int BL_kStride = fabbox.length(0) * fabbox.length(1);
        const int ijk_BoxLib = (i+BoxLib_ghosts) + (j+BoxLib_ghosts)*BL_jStride + (k+BoxLib_ghosts)*BL_kStride;

        fab_data[ijk_BoxLib] = level->my_boxes[box].vectors[VECTOR_U][ijk_HPGMG];
      }}}
  }
}
#endif /* USEHPGMG */
