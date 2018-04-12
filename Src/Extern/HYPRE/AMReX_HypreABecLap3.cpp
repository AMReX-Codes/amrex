#include <AMReX_HypreABecLap3.H>
#include <AMReX_HypreABec_F.H>
#include <cmath>

namespace {
static int ispow2(int i) {
  return (i == 1) ? 1 : (((i <= 0) || (i & 1)) ? 0 : ispow2(i / 2));
}
}

namespace amrex {

HypreABecLap3::HypreABecLap3(const BoxArray& grids,
                             const DistributionMapping& dmap,
                             const Geometry& geom_,
                             MPI_Comm comm_)
    : comm(comm_),
      geom(geom_),
      verbose(1),
      A(NULL), b(NULL), x(NULL) {

  const int ncomp = 1;
  int ngrow = 0;
  acoefs.define(grids, dmap, ncomp, ngrow);
  acoefs.setVal(0.0);

  for (int i = 0; i < BL_SPACEDIM; ++i) {
    BoxArray edge_boxes(grids);
    edge_boxes.surroundingNodes(i);
    bcoefs[i].define(edge_boxes, dmap, ncomp, ngrow);
    bcoefs[i].setVal(0.0);
  }

  bd.define(grids, dmap, ncomp, geom);

  int num_procs, myid;
  MPI_Comm_size(comm, &num_procs);
  MPI_Comm_rank(comm, &myid);

  // initialize with ngrow = 1
  ngrow = 1;

  // Store the global integer index
  GbInd.define(grids, dmap, ncomp, ngrow);
  GbInd.setVal(0);

  // Arrays needed to build the global indices for all procs
  CellsGIndex.resize(grids.size());
  numCellsProc.resize(num_procs+1);

  // These arrays store the starting global indices of all the boxes involved
  StartIndex.resize(grids.size());

  // Check if a face is at the physical boundary or interface between boxes
  FaceBcOffset.resize(grids.size()*BL_SPACEDIM*2);

  // Initialize all the arrays to 0
  std::fill(CellsGIndex.begin(), CellsGIndex.end(), 0);
  std::fill(numCellsProc.begin(), numCellsProc.end(), 0);
  std::fill(StartIndex.begin(), StartIndex.end(), 0);
  std::fill(FaceBcOffset.begin(), FaceBcOffset.end(), 0);

  // Fill the numpoints in the box and in the proc where this box resides
  for (int i = 0; i < grids.size(); i++) {
    const Box& bx = grids[i];
    CellsGIndex[i] = numCellsProc[ dmap[i]+1 ];
    numCellsProc[ dmap[i]+1 ] = numCellsProc[ dmap[i]+1 ] + bx.numPts();
  }

  // making sure that counting starts from 0 for procId = 0
  for (int i = 1; i< (num_procs+1); i++) {
    if (numCellsProc[i] == 0) {
      amrex::Abort("No rows in this core - Check domain decomposition");
    }
    numCellsProc[i] += numCellsProc[i-1];
  }

  // Box starting index starting from the proc offset
  for (int i = 0; i < grids.size(); i++) {
    CellsGIndex[i] = CellsGIndex[i] + numCellsProc[dmap[i]];
  }

  // Starting and ending global index for each box
  for (int i = 0; i < (grids.size()-1); i++) {
    StartIndex[i] = CellsGIndex[i];
  }
  StartIndex[grids.size()-1] = CellsGIndex[grids.size()-1];

  // Fill up the imultifab with global indices
  for (MFIter mfi(GbInd); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = mfi.validbox();

    // build the global index for all cells on this level
    amrex_BuildGlobalIndex(BL_TO_FORTRAN(GbInd[mfi]), ARLIM(reg.loVect()),
                           ARLIM(reg.hiVect()), CellsGIndex[i]);
  }

  const int nghost = 0;
  bool local = true;
  int ilower = GbInd.min(0, nghost, local);
  int iupper = GbInd.max(0, nghost, local);

  // Fill the global indices in the ghost cells along the grid edges
  GbInd.FillBoundary();

  // Create the HYPRE matrix object
  HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &A);

  // Parallel csr format storage
  HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

  // Initialize before setting coefficients
  HYPRE_IJMatrixInitialize(A);

  // Create the RHS and solution vector object
  HYPRE_IJVectorCreate(comm, ilower, iupper, &b);
  HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(b);

  HYPRE_IJVectorCreate(comm, ilower, iupper, &x);
  HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(x);
}

HypreABecLap3::~HypreABecLap3() {
  HYPRE_IJMatrixDestroy(A);
  HYPRE_IJVectorDestroy(b);
  HYPRE_IJVectorDestroy(x);
}

void HypreABecLap3::setScalars(Real sa, Real sb) {
  scalar_a = sa;
  scalar_b = sb;
}

void
HypreABecLap3::setACoeffs(const MultiFab& alpha) {
  MultiFab::Copy(acoefs, alpha, 0, 0, 1, 0);
}

void
HypreABecLap3::setBCoeffs(const std::array<const MultiFab*,
                          BL_SPACEDIM>& beta) {
  for (int idim=0; idim < BL_SPACEDIM; idim++) {
    MultiFab::Copy(bcoefs[idim], *beta[idim], 0, 0, 1, 0);
  }
}

void HypreABecLap3::setVerbose(int _verbose) {
  verbose = _verbose;
}

void
HypreABecLap3::solve(MultiFab& soln, const MultiFab& rhs,
                     Real rel_tol_, Real abs_tol_,
                     int max_iter_, LinOpBCType bc_type, Real bc_value) {
  loadBndryData(bc_type, bc_value);
  loadMatrix();
  finalizeMatrix();
  loadVectors(soln, rhs);
  finalizeVectors();
  setupSolver(rel_tol_, abs_tol_, max_iter_);
  solveDoIt();
  getSolution(soln);
  clearSolver();
}

void
HypreABecLap3::solve(MultiFab& soln, const MultiFab& rhs,
                     Real rel_tol_, Real abs_tol_,
                     int max_iter_, const BndryData& _bndry) {
  bd = _bndry;

  loadMatrix();
  finalizeMatrix();
  loadVectors(soln, rhs);
  finalizeVectors();
  setupSolver(rel_tol_, abs_tol_, max_iter_);
  solveDoIt();
  getSolution(soln);
  clearSolver();
}

void HypreABecLap3::loadBndryData(LinOpBCType bc_type, Real bc_value) {
  const int comp = 0;
  const Real* dx = geom.CellSize();
  for (int n=0; n< BL_SPACEDIM; ++n) {
    for (MFIter mfi(acoefs); mfi.isValid(); ++mfi) {
      int i = mfi.index();

      const Box& bx = mfi.validbox();

      // Our default will be that the face of this grid is
      // either touching another grid
      // across an interior boundary or a periodic boundary.
      // We will test for the other
      // cases below.
      {
        // Define the type of boundary conditions to be Dirichlet
        // (even for periodic)
        bd.setBoundCond(Orientation(n, Orientation::low), i,
                        comp, LO_DIRICHLET);
        bd.setBoundCond(Orientation(n, Orientation::high), i,
                        comp, LO_DIRICHLET);

        // Set the boundary conditions to the
        // cell centers outside the domain
        bd.setBoundLoc(Orientation(n, Orientation::low), i, 0.5*dx[n]);
        bd.setBoundLoc(Orientation(n, Orientation::high), i, 0.5*dx[n]);
      }

      // Now test to see if we should override
      // the above with Dirichlet or Neumann physical bc's
      if(bc_type != LinOpBCType::interior) {
        int ibnd = static_cast<int>(bc_type);
        // either LO_DIRICHLET or LO_NEUMANN

        // We are on the low side of the
        // domain in coordinate direction n

        if (bx.smallEnd(n) == geom.Domain().smallEnd(n)) {
          // Set the boundary conditions to
          // live exactly on the faces of the domain
          bd.setBoundLoc(Orientation(n, Orientation::low), i, 0.0);

          // Set the Dirichlet/Neumann boundary values
          bd.setValue(Orientation(n, Orientation::low), i, bc_value);

          // Define the type of boundary conditions
          bd.setBoundCond(Orientation(n, Orientation::low),
                          i, comp, ibnd);
        }

        // We are on the high side of the
        // domain in coordinate direction n
        if (bx.bigEnd(n) == geom.Domain().bigEnd(n)) {
          // Set the boundary conditions to
          // live exactly on the faces of the domain
          bd.setBoundLoc(Orientation(n, Orientation::high), i, 0.0);

          // Set the Dirichlet/Neumann boundary values
          bd.setValue(Orientation(n, Orientation::high), i, bc_value);

          // Define the type of boundary conditions
          bd.setBoundCond(Orientation(n, Orientation::high),
                          i, comp, ibnd);
        }
      }
    }
  }
}

void HypreABecLap3::loadMatrix() {
  static_assert(BL_SPACEDIM > 1, "HypreABecLap2: 1D not supported");

  const int size = 2 * BL_SPACEDIM + 1;
  const int bho = 0;
  const Real* dx = geom.CellSize();

  for (MFIter mfi(acoefs); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = mfi.validbox();
    const Box& domain = geom.Domain();

    for (OrientationIter oitr; oitr; oitr++) {
      int cdir(oitr());
      FaceBcOffset[i*BL_SPACEDIM*2+cdir] = 0;
      if (reg[oitr()] == domain[oitr()]) {
        FaceBcOffset[i*BL_SPACEDIM*2+cdir] = 1;
      }
    }
  }


  for (MFIter mfi(acoefs); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = mfi.validbox();

    int volume = reg.numPts();

    // build matrix interior

    amrex_hmac_ij(BL_TO_FORTRAN(acoefs[mfi]), ARLIM(reg.loVect()),
                  ARLIM(reg.hiVect()), scalar_a,
                  A, BL_TO_FORTRAN(GbInd[mfi]));

    int* BCface = FaceBcOffset.data();

    for (int idim = 0; idim < BL_SPACEDIM; idim++) {
      amrex_hmbc_ij(BL_TO_FORTRAN(bcoefs[idim][mfi]),
                    ARLIM(reg.loVect()), ARLIM(reg.hiVect()), scalar_b,
                    geom.CellSize(), idim, A, BL_TO_FORTRAN(GbInd[mfi]));
    }

    // add b.c.'s to matrix diagonal, and
    // zero out offdiag values at domain boundaries

    const Vector< Vector<BoundCond> > & bcs_i = bd.bndryConds(i);
    const BndryData::RealTuple      & bcl_i = bd.bndryLocs(i);

    const Box& domain = geom.Domain();
    for (OrientationIter oitr; oitr; oitr++) {
      int cdir(oitr());
      int idim = oitr().coordDir();
      const int bctype = bcs_i[cdir][0];
      const Real &bcl  = bcl_i[cdir];
      const Mask &msk  = bd.bndryMasks(oitr())[mfi];

      // Treat an exposed grid edge here as a boundary condition
      // for the linear solver:

      if (reg[oitr()] == domain[oitr()]) {
        amrex_hmmat3_ij(ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
                        cdir, bctype, bho, bcl,
                        BL_TO_FORTRAN(msk),
                        BL_TO_FORTRAN(bcoefs[idim][mfi]),
                        scalar_b, dx, A, BL_TO_FORTRAN(GbInd[mfi]));
      } else {
        amrex_hmmat_ij(ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
                       cdir, bctype, bho, bcl,
                       BL_TO_FORTRAN(msk),
                       BL_TO_FORTRAN(bcoefs[idim][mfi]),
                       scalar_b, dx, A, BL_TO_FORTRAN(GbInd[mfi]));
      }
    }
  }
}

void HypreABecLap3::finalizeMatrix() {
  // Assemble after setting the coefficients
  HYPRE_IJMatrixAssemble(A);
}

void HypreABecLap3::loadVectors(MultiFab& soln, const MultiFab& rhs) {
  const int bho = 0;
  const Real* dx = geom.CellSize();

  FArrayBox fnew;

  soln.setVal(0.0);

  for (MFIter mfi(soln); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = mfi.validbox();

    // initialize soln, since we will reuse the space to set up rhs below:

    FArrayBox *f;
    if (soln.nGrow() == 0) {  // need a temporary if soln is the wrong size
      f = &soln[mfi];
    } else {
      f = &fnew;
      f->resize(reg);
      f->copy(soln[mfi], 0, 0, 1);
    }
    // sharing space, soln will be overwritten below
    Real* vec = f->dataPtr();

    // Convert the 3D vec array to 1D array with indices
    // corresponding to row indices in the matrix
    amrex_conv_Vec_Local_Global(x, vec, reg.numPts(),
                                ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
                                BL_TO_FORTRAN(GbInd[mfi]));

    // Copy the rhs vector
    f->copy(rhs[mfi], 0, 0, 1);

    // add b.c.'s to rhs
    const Vector< Vector<BoundCond> > & bcs_i = bd.bndryConds(i);
    const BndryData::RealTuple      & bcl_i = bd.bndryLocs(i);

    const Box& domain = geom.Domain();
    for (OrientationIter oitr; oitr; oitr++) {
      int cdir(oitr());
      int idim = oitr().coordDir();
      const int bctype = bcs_i[cdir][0];
      const Real &bcl  = bcl_i[cdir];

      const Mask &msk  = bd.bndryMasks(oitr())[mfi];
      const FArrayBox &fs = bd.bndryValues(oitr())[mfi];

      // Treat an exposed grid edge here as a boundary condition
      // for the linear solver:
      if (reg[oitr()] == domain[oitr()]) {
        amrex_hbvec3(vec, ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
                     cdir, bctype, bho, bcl,
                     BL_TO_FORTRAN(fs),
                     BL_TO_FORTRAN(msk),
                     BL_TO_FORTRAN(bcoefs[idim][mfi]),
                     scalar_b, dx);
      } else {
        amrex_hbvec(vec, ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
                    cdir, bctype, bho, bcl,
                    BL_TO_FORTRAN(fs),
                    BL_TO_FORTRAN(msk),
                    BL_TO_FORTRAN(bcoefs[idim][mfi]),
                    scalar_b, dx);
      }
    }

    // initialize rhs
    amrex_conv_Vec_Local_Global(b, vec, reg.numPts(),
                                ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
                                BL_TO_FORTRAN(GbInd[mfi]));
  }
}

void HypreABecLap3::finalizeVectors() {
  HYPRE_IJVectorAssemble(x);
  HYPRE_IJVectorAssemble(b);
}

void HypreABecLap3::setupSolver(Real rel_tol_, Real abs_tol_, int max_iter_) {
  rel_tol = rel_tol_;
  abs_tol = abs_tol_;
  max_iter = max_iter_;

  HYPRE_BoomerAMGCreate(&solver);

  HYPRE_BoomerAMGSetMinIter(solver, 1);
  HYPRE_BoomerAMGSetMaxIter(solver, max_iter);
  HYPRE_BoomerAMGSetTol(solver, rel_tol);

  // Set some parameters (See Reference Manual for more parameters)
  // Falgout coarsening with modified classical interpolation
  HYPRE_BoomerAMGSetOldDefault(solver);
  HYPRE_BoomerAMGSetCoarsenType(solver, 6);
  HYPRE_BoomerAMGSetCycleType(solver, 1);
  HYPRE_BoomerAMGSetRelaxType(solver, 6);   /* G-S/Jacobi hybrid relaxation */
  HYPRE_BoomerAMGSetRelaxOrder(solver, 1);   /* uses C/F relaxation */
  HYPRE_BoomerAMGSetNumSweeps(solver, 2);   /* Sweeeps on each level */
  HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
  HYPRE_BoomerAMGSetStrongThreshold(solver, 0.6);

  // Get the parcsr matrix object to use

  HYPRE_IJMatrixGetObject(A, (void**)  &par_A);
  HYPRE_IJVectorGetObject(b, (void **) &par_b);
  HYPRE_IJVectorGetObject(x, (void **) &par_x);

  HYPRE_BoomerAMGSetup(solver, par_A, par_b, par_x);
}

void HypreABecLap3::clearSolver() {
  HYPRE_BoomerAMGDestroy(solver);
}

void HypreABecLap3::solveDoIt() {
  if (abs_tol > 0.0) {
    Real bnorm;
    bnorm = hypre_ParVectorInnerProd(par_b, par_b);
    bnorm = std::sqrt(bnorm);

    const BoxArray& grids = acoefs.boxArray();
    Real volume = grids.numPts();
    Real rel_tol_new = bnorm > 0.0 ? abs_tol / bnorm *
        std::sqrt(volume) : rel_tol;

    if (rel_tol_new > rel_tol) {
      HYPRE_BoomerAMGSetTol(solver, rel_tol_new);
    }
  }

  HYPRE_BoomerAMGSolve(solver, par_A, par_b, par_x);
  if (verbose >= 2 && ParallelDescriptor::IOProcessor()) {
    int num_iterations;
    Real res;
    HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &res);

    int oldprec = std::cout.precision(17);
    std::cout <<"\n" <<  num_iterations
              << " Hypre Multigrid Iterations, Relative Residual "
              << res << std::endl;
    std::cout.precision(oldprec);
  }
}

void HypreABecLap3::getSolution(MultiFab& soln) {
  const int part = 0;

  std::vector<int> VecIndices;

  FArrayBox fnew;
  for (MFIter mfi(soln); mfi.isValid(); ++mfi) {
    const Box &reg = mfi.validbox();

    FArrayBox *f;
    if (soln.nGrow() == 0) {
      // need a temporary if soln is the wrong size
      f = &soln[mfi];
    } else {
      f = &fnew;
      f->resize(reg);
    }

    int i = mfi.index();

    // Storage for the solution vector returned by HYPRE
    Real *VecGB = hypre_CTAlloc(double, reg.numPts(), HYPRE_MEMORY_HOST);

    // Generate indices corresponding to all the boxes
    VecIndices.resize(reg.numPts());
    std::iota(VecIndices.begin(), VecIndices.end(), StartIndex[i]);
    int* RowIndices = VecIndices.data();

    // Get the solution from HYPRE
    HYPRE_IJVectorGetValues(x, reg.numPts(), RowIndices, VecGB);

    // Move the solution vector back to the soln multifab
    amrex_conv_Vec_Global_Local(BL_TO_FORTRAN(*f), VecGB,
                                reg.numPts(), ARLIM(reg.loVect()),
                                ARLIM(reg.hiVect()));

    if (soln.nGrow() != 0) {
      soln[mfi].copy(*f, 0, 0, 1);
    }

    hypre_TFree(VecGB, HYPRE_MEMORY_HOST);
  }
}

}  // namespace amrex
