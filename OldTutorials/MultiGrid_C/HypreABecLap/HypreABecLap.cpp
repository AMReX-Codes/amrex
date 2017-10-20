
#include <HypreABecLap.H>
#include <HypreABec_F.H>

#include <_hypre_struct_mv.h>

HypreABecLap::HypreABecLap(const BoxArray& grids, const Geometry& geom)
{
  ParmParse pp("hypre");

  // solver_flag = 0 for SMG
  // solver_flag = 1 for PFMG
  solver_flag = 1;
  pp.query("solver_flag", solver_flag);

  int num_procs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs );
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );

  for (int i = 0; i < BL_SPACEDIM; i++) {
    dx[i] = geom.CellSize(i);
  }

  HYPRE_StructGridCreate(MPI_COMM_WORLD, BL_SPACEDIM, &grid);

  for (int i = 0; i < BL_SPACEDIM; i++) {
    is_periodic[i] = 0;
    if (geom.isPeriodic(i)) {
      is_periodic[i] = geom.period(i);
      BL_ASSERT(ispow2(is_periodic[i]));
      BL_ASSERT(geom.Domain().smallEnd(i) == 0);
    }
  }
  if (geom.isAnyPeriodic()) {
    HYPRE_StructGridSetPeriodic(grid, is_periodic);
  }

  if (num_procs != 1) {
    // parallel section:
    BL_ASSERT(ParallelDescriptor::NProcs() == num_procs);
    BL_ASSERT(ParallelDescriptor::MyProc() == myid);
    DistributionMapping distributionMap(grids, num_procs);

    for (int i = 0; i < grids.size(); i++) {
      if (distributionMap[i] == myid) {
	HYPRE_StructGridSetExtents(grid, loV(grids[i]), hiV(grids[i]));
      }
    }
  }
  else {
    for (int i = 0; i < grids.size(); i++) {
      HYPRE_StructGridSetExtents(grid, loV(grids[i]), hiV(grids[i]));
    }
  }

  HYPRE_StructGridAssemble(grid);

#if (BL_SPACEDIM == 2)
  int offsets[5][2] = {{-1,  0},      // 0
		       { 0, -1},      // 1
		       { 1,  0},      // 2
		       { 0,  1},      // 3
		       { 0,  0}};     // 4
#elif (BL_SPACEDIM == 3)
  int offsets[7][3] = {{-1,  0,  0},  // 0
		       { 0, -1,  0},  // 1
		       { 0,  0, -1},  // 2
		       { 1,  0,  0},  // 3
		       { 0,  1,  0},  // 4
		       { 0,  0,  1},  // 5
		       { 0,  0,  0}}; // 6
#endif

  HYPRE_StructStencil stencil;

  HYPRE_StructStencilCreate(BL_SPACEDIM, 2 * BL_SPACEDIM + 1, &stencil);

  for (int i = 0; i < 2 * BL_SPACEDIM + 1; i++) {
    HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
  }

  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
  HYPRE_StructMatrixInitialize(A);

  HYPRE_StructStencilDestroy(stencil); 

  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);

  HYPRE_StructVectorInitialize(b);
  HYPRE_StructVectorInitialize(x);

  int ncomp=1;
  int ngrow=0;
  acoefs = new MultiFab(grids, ncomp, ngrow);
  acoefs->setVal(0.0);
 
  for (int i = 0; i < BL_SPACEDIM; i++) {
    BoxArray edge_boxes(grids);
    edge_boxes.surroundingNodes(i);
    bcoefs[i] = new MultiFab(edge_boxes, ncomp, ngrow);
  }

  pfmg_rap_type = -1;
  pp.query("pfmg_rap_type", pfmg_rap_type);

  pfmg_relax_type = -1;
  pp.query("pfmg_relax_type", pfmg_relax_type);

  num_pre_relax = -1;
  pp.query("num_pre_relax", num_pre_relax);
  num_post_relax = -1;
  pp.query("num_post_relax", num_post_relax);

  skip_relax = -1;
  pp.query("skip_relax", skip_relax);

  print_level = -1;
  pp.query("print_level", print_level);
}

HypreABecLap::~HypreABecLap()
{
  HYPRE_StructGridDestroy(grid);
  HYPRE_StructMatrixDestroy(A);
  HYPRE_StructVectorDestroy(b);
  HYPRE_StructVectorDestroy(x);

  delete acoefs;
  for (int i = 0; i < BL_SPACEDIM; i++) {
    delete bcoefs[i];
  }
}

void HypreABecLap::setScalars(Real sa, Real sb)
{
  scalar_a = sa;
  scalar_b = sb;
}

void HypreABecLap::setACoeffs(const MultiFab& alpha)
{
  BL_ASSERT( alpha.ok() );
  BL_ASSERT( alpha.boxArray() == acoefs->boxArray() );

  acoefs->copy(alpha);
}

void HypreABecLap::setBCoeffs(const MultiFab beta[])
{
  for (int idim=0; idim<BL_SPACEDIM; idim++) {
    BL_ASSERT( beta[idim].ok() );
    BL_ASSERT( beta[idim].boxArray() == bcoefs[idim]->boxArray() );

    bcoefs[idim]->copy(beta[idim]);
  }
}

void HypreABecLap::setBCoeffs(const PArray<MultiFab>& beta)
{
  for (int idim=0; idim<BL_SPACEDIM; idim++) {
    BL_ASSERT( beta[idim].ok() );
    BL_ASSERT( beta[idim].boxArray() == bcoefs[idim]->boxArray() );

    bcoefs[idim]->copy(beta[idim]);
  }
}

void HypreABecLap::setVerbose(int _verbose)
{
  verbose = _verbose;
}

void HypreABecLap::solve(MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol,
			 int max_iter, const BndryData& bndry)
{
  reltol = rel_tol;
  abstol = abs_tol;
  maxiter = max_iter;

  // set up solver
  const BoxArray& grids = acoefs->boxArray();

  const int size = 2 * BL_SPACEDIM + 1;
  int stencil_indices[size];
  for (int i = 0; i < size; i++) {
    stencil_indices[i] = i;
  }

  soln.setVal(0.0);

  Real *mat;
  Real *vec;
  for (MFIter mfi(soln); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[i];

    FArrayBox *f;
    if (soln.nGrow() == 0) { // need a temporary if soln is the wrong size
      f = &soln[mfi];
    }
    else {
      f = new FArrayBox(reg);
      f->copy(soln[mfi], 0, 0, 1);
    }
    vec = f->dataPtr(); // sharing space, soln will be overwritten below

    HYPRE_StructVectorSetBoxValues(x, loV(reg), hiV(reg), vec);

    f->copy(rhs[mfi], 0, 0, 1); 
    // vec now contains rhs, but we need to add bc's before SetBoxValues

    int volume = reg.numPts();
    mat = hypre_CTAlloc(double, size*volume);

    // build matrix interior
    const int* alo = (*acoefs)[mfi].loVect();
    const int* ahi = (*acoefs)[mfi].hiVect();
    FORT_HPACOEF(mat, (*acoefs)[mfi].dataPtr(), ARLIM(alo), ARLIM(ahi),
		 reg.loVect(), reg.hiVect(), scalar_a);

    for (int idim = 0; idim < BL_SPACEDIM; idim++) {
      const int* blo = (*bcoefs[idim])[mfi].loVect();
      const int* bhi = (*bcoefs[idim])[mfi].hiVect();
      FORT_HPBCOEF(mat, (*bcoefs[idim])[mfi].dataPtr(), ARLIM(blo), ARLIM(bhi),
		   reg.loVect(), reg.hiVect(), scalar_b, dx, idim);
    }

    const Vector< Vector<BoundCond> > & bcs_i = bndry.bndryConds(i);
    const BndryData::RealTuple      & bcl_i = bndry.bndryLocs(i);
    const BndryData::MaskTuple      & msk_i = bndry.bndryMasks(i);

    // add b.c.'s for A matrix and b vector
    const Box& domain = bndry.getDomain();
    for (OrientationIter oitr; oitr; oitr++) {
      int cdir(oitr());
      int idim = oitr().coordDir();
      const BoundCond &bct = bcs_i[cdir][0];
      const Real      &bcl = bcl_i[cdir];
      const Mask      &msk = *msk_i[cdir];
      const FArrayBox &bcv = bndry.bndryValues(oitr())[mfi];
      const int* blo = (*bcoefs[idim])[mfi].loVect();
      const int* bhi = (*bcoefs[idim])[mfi].hiVect();
      const int* mlo = msk.loVect();
      const int* mhi = msk.hiVect();
      const int* bvlo = bcv.loVect();
      const int* bvhi = bcv.hiVect();
      
      if (reg[oitr()] == domain[oitr()] && is_periodic[idim] == 0) {
	int bctype = bct;
	FORT_HPBVEC3(vec, (*bcoefs[idim])[mfi].dataPtr(), ARLIM(blo), ARLIM(bhi),
		     reg.loVect(), reg.hiVect(), scalar_b, dx, cdir, bctype, bcl, 
		     msk.dataPtr(), ARLIM(mlo), ARLIM(mhi),
		     bcv.dataPtr(), ARLIM(bvlo), ARLIM(bvhi));
	FORT_HPBMAT3(mat, (*bcoefs[idim])[mfi].dataPtr(), ARLIM(blo), ARLIM(bhi),
		     reg.loVect(), reg.hiVect(), scalar_b, dx, cdir, bctype, bcl, 
		     msk.dataPtr(), ARLIM(mlo), ARLIM(mhi));
      }
    }

    // initialize matrix & vectors
    HYPRE_StructMatrixSetBoxValues(A, loV(reg), hiV(reg),
				   size, stencil_indices, mat);
    HYPRE_StructVectorSetBoxValues(b, loV(reg), hiV(reg), vec);

    hypre_TFree(mat);

    if (soln.nGrow() != 0) {
      delete f;       // contains vec
    }    
  }

  HYPRE_StructMatrixAssemble(A);
  HYPRE_StructVectorAssemble(x); 
  HYPRE_StructVectorAssemble(b); 

  if (solver_flag == 0) {
    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructSMGSetMemoryUse(solver, 0);
    HYPRE_StructSMGSetMaxIter(solver, maxiter);
    HYPRE_StructSMGSetRelChange(solver, 0);
    HYPRE_StructSMGSetTol(solver, reltol);
    HYPRE_StructSMGSetNumPreRelax(solver, 1);
    HYPRE_StructSMGSetNumPostRelax(solver, 1);
    HYPRE_StructSMGSetLogging(solver, 1);
    HYPRE_StructSMGSetup(solver, A, b, x);
  }
  else if (solver_flag == 1) {
    HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &solver);
    if (skip_relax >= 0) {
      HYPRE_StructPFMGSetSkipRelax(solver, skip_relax); // default: 1
    }
    HYPRE_StructPFMGSetMaxIter(solver, maxiter);
    //    HYPRE_StructPFMGSetRelChange(solver, 0);    // default: 0
    HYPRE_StructPFMGSetTol(solver, reltol);           // default: 1.e-6
    if (pfmg_rap_type >= 0) {
      HYPRE_StructPFMGSetRAPType(solver, pfmg_rap_type);  // default: 1
    }
    if (pfmg_relax_type >= 0) {
      HYPRE_StructPFMGSetRelaxType(solver, pfmg_relax_type);  // default: 1
    }
    if (num_pre_relax >= 0) {
      HYPRE_StructPFMGSetNumPreRelax(solver, num_pre_relax);  // default: 1
    }
    if (num_post_relax <= 0) {
      HYPRE_StructPFMGSetNumPostRelax(solver, num_post_relax);// default: 1
    }
    HYPRE_StructPFMGSetLogging(solver, 1);                  // default: 0
    if (print_level >= 0) {
      HYPRE_StructPFMGSetPrintLevel(solver, print_level);     // default: 0
    }
    HYPRE_StructPFMGSetup(solver, A, b, x);
  }

  if (abstol > 0.0) {
    Real bnorm;
    bnorm = hypre_StructInnerProd((hypre_StructVector *) b,
				  (hypre_StructVector *) b);
    bnorm = sqrt(bnorm);

    const BoxArray& grids = acoefs->boxArray();
    Real volume = 0.0;
    for (int i = 0; i < grids.size(); i++) {
      volume += grids[i].numPts();
    }

    Real reltol_new = (bnorm > 0.0
		       ? abstol / bnorm * sqrt(volume)
		       : reltol);

    if (reltol_new > reltol) {
      if (solver_flag == 0) {
	HYPRE_StructSMGSetTol(solver, reltol_new);
      }
      else if(solver_flag == 1) {
	HYPRE_StructPFMGSetTol(solver, reltol_new);
      }
    }
  }

  if (solver_flag == 0) {
    HYPRE_StructSMGSolve(solver, A, b, x);
  }
  else if (solver_flag == 1) {
    HYPRE_StructPFMGSolve(solver, A, b, x);
  }

  for (MFIter mfi(soln); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[i];

    FArrayBox *f;
    if (soln.nGrow() == 0) { // need a temporary if soln is the wrong size
      f = &soln[mfi];
    }
    else {
      f = new FArrayBox(reg);
    }

    vec = f->dataPtr();
    HYPRE_StructVectorGetBoxValues(x, loV(reg), hiV(reg), vec);

    if (soln.nGrow() != 0) {
      soln[mfi].copy(*f, 0, 0, 1);
      delete f;
    }
  }

  if (verbose >= 2 && ParallelDescriptor::IOProcessor()) {
    int num_iterations;
    Real res;
    if (solver_flag == 0) {
      HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
      HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &res);
    }
    else if(solver_flag == 1) {
      HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
      HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &res);
    }

    int oldprec = std::cout.precision(20);
    std::cout << num_iterations
	      << " Hypre Multigrid Iterations, Relative Residual "
	      << res << std::endl;
    std::cout.precision(oldprec);
  }

  // claer solver
  if (solver_flag == 0) {
    HYPRE_StructSMGDestroy(solver);
  }
  else if (solver_flag == 1) {
    HYPRE_StructPFMGDestroy(solver);
  }
}

