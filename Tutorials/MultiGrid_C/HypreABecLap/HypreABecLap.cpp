
#include "HypreABecLap.H"
#include "HypreABec_F.H"

#include "_hypre_struct_mv.h"

static int ispow2(int i)
{
  return (i == 1) ? 1 : (((i <= 0) || (i & 1)) ? 0 : ispow2(i / 2));
}

static int* loV(const Box& b) {
#if (BL_SPACEDIM == 1)
  vl[0] = b.smallEnd(0);
  return vl;
#else
  return (int*) b.loVect();
#endif
}

static int* hiV(const Box& b) {
#if (BL_SPACEDIM == 1)
  vh[0] = b.bigEnd(0);
  return vh;
#else
  return (int*) b.hiVect();
#endif
}

HypreABecLap::HypreABecLap(const BoxArray& grids, const Geometry& geom,
			   int _solver_flag)
  : solver_flag(_solver_flag)
{
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

#if (BL_SPACEDIM == 1)
  // if we were really 1D:
/*
  int offsets[2][1] = {{-1},
		       { 0}};
*/
  // fake 1D as a 2D problem:
  int offsets[2][2] = {{-1,  0},
		       { 0,  0}};
#elif (BL_SPACEDIM == 2)
  int offsets[3][2] = {{-1,  0},
		       { 0, -1},
		       { 0,  0}};
#elif (BL_SPACEDIM == 3)
  int offsets[4][3] = {{-1,  0,  0},
		       { 0, -1,  0},
		       { 0,  0, -1},
		       { 0,  0,  0}};
#endif

#if   (BL_SPACEDIM == 1)
  int A_num_ghost[6] = { 1, 1, 0, 0, 0, 0 };
#elif (BL_SPACEDIM == 2)
  //int A_num_ghost[4] = { 1, 1, 1, 1 };
  int A_num_ghost[6] = { 1, 1, 1, 1, 0, 0 };
#elif (BL_SPACEDIM == 3)
  int A_num_ghost[6] = { 1, 1, 1, 1, 1, 1 };
#endif

  HYPRE_StructStencil stencil;

#if (BL_SPACEDIM == 1)
  HYPRE_StructStencilCreate(2, 2, &stencil);
#else
  HYPRE_StructStencilCreate(BL_SPACEDIM, BL_SPACEDIM + 1, &stencil);
#endif

  for (int i = 0; i < BL_SPACEDIM + 1; i++) {
    HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
  }

  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
  HYPRE_StructMatrixSetSymmetric(A, 1);
  HYPRE_StructMatrixSetNumGhost(A, A_num_ghost);
  HYPRE_StructMatrixInitialize(A);

  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);

  HYPRE_StructStencilDestroy(stencil); // no longer needed

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
  scalar_b  = sb;
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

void HypreABecLap::setVerbose(int _verbose)
{
  verbose = _verbose;
}

void HypreABecLap::solve(MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol,
			 int max_iter, const BndryData& bndry)
{
  solver_flag = 1;

  reltol = rel_tol;
  abstol = abs_tol;
  maxiter = max_iter;

  // set up solver
  const BoxArray& grids = acoefs->boxArray();

  const int size = BL_SPACEDIM + 1;
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
      f = &soln[i];
    }
    else {
      f = new FArrayBox(reg);
      f->copy(soln[i], 0, 0, 1);
    }
    vec = f->dataPtr(); // sharing space, soln will be overwritten below

    HYPRE_StructVectorSetBoxValues(x, loV(reg), hiV(reg), vec);

    f->copy(rhs[i], 0, 0, 1); 
    // vec now contains rhs, but we nned to add bc's before SetBoxValues

    int volume = reg.numPts();
    mat = hypre_CTAlloc(double, size*volume);

    // build matrix interior
    const int* alo = (*acoefs)[i].loVect();
    const int* ahi = (*acoefs)[i].hiVect();
    FORT_HPACOEF(mat, (*acoefs)[i].dataPtr(), ARLIM(alo), ARLIM(ahi),
		 reg.loVect(), reg.hiVect(), scalar_a);

    for (int idim = 0; idim < BL_SPACEDIM; idim++) {
      const int* blo = (*bcoefs[idim])[i].loVect();
      const int* bhi = (*bcoefs[idim])[i].hiVect();
      FORT_HPBCOEF(mat, (*bcoefs[idim])[i].dataPtr(), ARLIM(blo), ARLIM(bhi),
		   reg.loVect(), reg.hiVect(), scalar_b, dx, idim);
    }

    // add b.c.'s for A matrix and b vector
    const Box& domain = bndry.getDomain();
    for (OrientationIter oitr; oitr; oitr++) {
      int cdir(oitr());
      int idim = oitr().coordDir();
      if (is_periodic[idim] == 0) {
	const BoundCond &bct = bndry.bndryConds(oitr())[i][0];
	const Real      &bcl = bndry.bndryLocs(oitr())[i];
	const FArrayBox &bcv = bndry.bndryValues(oitr())[i];
	const Mask      &msk = bndry.bndryMasks(oitr())[i];
	const int* blo = (*bcoefs[idim])[i].loVect();
	const int* bhi = (*bcoefs[idim])[i].hiVect();
	const int* mlo = msk.loVect();
	const int* mhi = msk.hiVect();
	const int* bvlo = bcv.loVect();
	const int* bvhi = bcv.hiVect();
	
	if (reg[oitr()] == domain[oitr()]) {
	  int bctype = bct;
	  FORT_HPBVEC3(vec, (*bcoefs[idim])[i].dataPtr(), ARLIM(blo), ARLIM(bhi),
		       reg.loVect(), reg.hiVect(), scalar_b, dx, cdir, bctype, bcl, 
		       msk.dataPtr(), ARLIM(mlo), ARLIM(mhi),
		       bcv.dataPtr(), ARLIM(bvlo), ARLIM(bvhi));
	  FORT_HPBMAT3(mat, (*bcoefs[idim])[i].dataPtr(), ARLIM(blo), ARLIM(bhi),
		       reg.loVect(), reg.hiVect(), scalar_b, dx, cdir, bctype, bcl, 
		       msk.dataPtr(), ARLIM(mlo), ARLIM(mhi));
	}
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

  if (solver_flag == 1) {
    HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &solver);
    //HYPRE_StructPFMGSetMemoryUse(solver, 0);
    HYPRE_StructPFMGSetSkipRelax(solver, 0);
    HYPRE_StructPFMGSetMaxIter(solver, maxiter);
    HYPRE_StructPFMGSetRelChange(solver, 0);
    HYPRE_StructPFMGSetTol(solver, reltol);
    // in following line, Falgout says use 1 as relax type, not 2 (rbp, 9/27/05)
    // weighted Jacobi = 1; red-black GS = 2
    int pfmg_relax_type = 1;
    HYPRE_StructPFMGSetRelaxType(solver, pfmg_relax_type);
    HYPRE_StructPFMGSetNumPreRelax(solver, 1);
    HYPRE_StructPFMGSetNumPostRelax(solver, 1);
    HYPRE_StructPFMGSetLogging(solver, 1);
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
      else if(solver_flag == 2) {
	// nothing for this option
      }
      else if(solver_flag == 3 || solver_flag == 4) {
	HYPRE_StructPCGSetTol(solver, reltol_new);
      }
    }
  }

  if (solver_flag == 0) {
    HYPRE_StructSMGSolve(solver, A, b, x);
  }
  else if (solver_flag == 1) {
    HYPRE_StructPFMGSolve(solver, A, b, x);
  }
  else if (solver_flag == 2) {
    HYPRE_StructJacobiSolve(solver, A, b, x);
  }
  else if (solver_flag == 3 || solver_flag == 4) {
    HYPRE_StructPCGSolve(solver, A, b, x);
  }
  else if (solver_flag == 5 || solver_flag == 6) {
    HYPRE_StructHybridSolve(solver, A, b, x);
  }

  for (MFIter mfi(soln); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[i];

    FArrayBox *f;
    if (soln.nGrow() == 0) { // need a temporary if soln is the wrong size
      f = &soln[i];
    }
    else {
      f = new FArrayBox(reg);
    }

    vec = f->dataPtr();
    HYPRE_StructVectorGetBoxValues(x, loV(reg), hiV(reg), vec);

    if (soln.nGrow() != 0) {
      soln[i].copy(*f, 0, 0, 1);
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
    else if(solver_flag == 2) {
      HYPRE_StructJacobiGetNumIterations(solver, &num_iterations);
      HYPRE_StructJacobiGetFinalRelativeResidualNorm(solver, &res);
    }
    else if(solver_flag == 3 || solver_flag == 4) {
      HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
      HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &res);
    }
    else if(solver_flag == 5 || solver_flag == 6) {
      HYPRE_StructHybridGetNumIterations(solver, &num_iterations);
      HYPRE_StructHybridGetFinalRelativeResidualNorm(solver, &res);
    }

    int oldprec = std::cout.precision(20);
    std::cout << num_iterations
	      << " Hypre Multigrid Iterations, Relative Residual "
	      << res << std::endl;
    std::cout.precision(oldprec);
  }

  // claer solver
  if (solver_flag == 1) {
    HYPRE_StructPFMGDestroy(solver);
  }
}

