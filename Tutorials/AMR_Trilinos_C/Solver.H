#ifndef _SOLVER_H
#define _SOLVER_H

#include "ml_include.h"
#include <Utility.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <MultiFab.H>
#include "Epetra_MpiComm.h"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Operator.h"
#include "EpetraExt_RowMatrixOut.h"
#include <Epetra_Import.h>

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosStatusTestGenResNorm.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include <Teuchos_ParameterList.hpp>

#include "ml_MultiLevelPreconditioner.h"
#include "ml_MultiLevelOperator.h"
#include "ml_epetra_utils.h"

#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>

#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <deque>
#include <fstream>
#include <sstream>
#include <iostream>

//using namespace Teuchos;
using Teuchos::RCP;
using Teuchos::rcp;
using namespace ML_Epetra;
using namespace Isorropia;

class Solver {

public:
    /** Constructor for the SAAMGPCG solver.
     * \param hr the mesh-spacings in each direction.
     * \param Comm the communicator object used in the solver routines.
     * \param verbose flag if output should be verbose.
     * \sa ~Solver(), SetupMLList() and SetupProblem()
     */
    Solver(const Box& domain_in, const std::vector<double>& hr_in, 
           Epetra_MpiComm& Comm, bool verbose, 
           double tol, int maxIterations, int numBlocks, int recycleBlocks, 
           int maxOldLHS, MultiFab& rhs, MultiFab& soln) : 

          Comm_m(Comm)
          {
            Diag = 0;

            domain = domain_in;
            dom_lo = domain.loVect();
            dom_hi = domain.hiVect();

            nx = domain.hiVect()[0] - domain.loVect()[0] + 1;
            ny = domain.hiVect()[1] - domain.loVect()[1] + 1;
#if (BL_SPACEDIM == 3)
            nz = domain.hiVect()[2] - domain.loVect()[2] + 1;
#endif

            hr = hr_in;
            verbose_m = verbose;

            isReusingHierarchy_m = false;
            isReusingPreconditioner_m = false;
            //XXX: currently set reusing to hierarchy
            isReusingHierarchy_m = true;

            numBlocks_m = numBlocks;
            recycleBlocks_m = recycleBlocks;
            nLHS_m = maxOldLHS;

            // setup boundary treatment
            InitializeDomainGeometry();

            // setup ml preconditioner parameters
            SetupMLList();

            SetupProblem(domain, rhs, soln);

            // setup extrapolation
	    if(nLHS_m > 0)
		    P = rcp(new Epetra_MultiVector(*Map, nLHS_m, false));

	    MLPrec = Teuchos::null;

	    // setup Belos parameters
	    belosList.set( "Maximum Iterations", maxIterations );  // Maximum number of iterations allowed
	    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
	    if(numBlocks_m != 0 && recycleBlocks_m != 0){          // only set if solver==RCGSolMgr
		    belosList.set("Num Blocks", numBlocks_m);          // Maximum number of blocks in Krylov space
		    belosList.set("Num Recycled Blocks", recycleBlocks_m); // Number of vectors in recycle space
	    }
	    if(verbose_m) {
		    belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails);
		    belosList.set("Output Frequency", 1);
	    } else
		    belosList.set("Verbosity", Belos::Errors + Belos::Warnings);

	    // setup Belos solver
	    if(numBlocks_m == 0 || recycleBlocks_m == 0)
		    solver = rcp( new Belos::BlockCGSolMgr<double,MV,OP>() );
	    else
		    solver = rcp( new Belos::RCGSolMgr<double,MV,OP>() );
	    convStatusTest = rcp( new Belos::StatusTestGenResNorm<ST,MV,OP> (tol) );
	    convStatusTest->defineScaleForm(Belos::NormOfRHS, Belos::TwoNorm);
#ifdef UserConv
	    solver->setUserConvStatusTest(convStatusTest);
#endif
    }
      /** Destructor for the SAAMGPCG solver.
     */
    ~Solver() {
	    delete Diag;
	    delete Map;
    }

    /** Actual computation of self field.
     * \param repetitions the number of calls to the solver.
     * \param tol the convergence criterion for the PCG solver.
     * \param maxIterations the maximal number of iterations performed by the PCG solver.
     */
    void Compute();

    /** Performs setup of data distribution and problem (system matrix, RHS, LHS).
     */
    void SetupProblem(const Box& domain, MultiFab& rhs, MultiFab& soln);

    /** Copy LHS->Values() into soln
     */
    void CopySolution(Box& domain, MultiFab& soln);

    int getNumIters();

private:

    /// Map corresponding to data distribution
    Epetra_Map* Map;

    /// preconditioner object
    RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec;//now RCP TODO

    /// parameter list used for the preconditioner (filled in SetupMLList)
    Teuchos::ParameterList MLList_m;

    /// system matrix
    RCP<Epetra_CrsMatrix> A;

    /// diagonal matrix of A
    Epetra_CrsMatrix* Diag;

    /// right hand side of our problem
    RCP<Epetra_Vector> RHS;

    /// left hand side
    RCP<Epetra_Vector> LHS;
    RCP<Epetra_MultiVector> lhssol;

    /// last N LHS's for extrapolating the new LHS as starting vector
    //uint nLHS_m;
    int nLHS_m;
    RCP<Epetra_MultiVector> P;
    std::deque< Epetra_Vector > OldLHS;

    /// Solver (Belos BlockCG or RCG)
    /// if numBlocks_m OR recycleBlocks_m == 0 BlockCG is used, else RCG
    /// maximum number of blocks in Krylov space
    int numBlocks_m;
    /// number of vectors in recycle space
    int recycleBlocks_m;

    typedef double                          ST;
    typedef Epetra_Operator                 OP;
    typedef Epetra_MultiVector              MV;
    typedef Belos::OperatorTraits<ST,MV,OP> OPT;// TODO
    typedef Belos::MultiVecTraits<ST,MV>    MVT;// TODO

    Teuchos::ParameterList belosList;
    Belos::LinearProblem<ST,MV,OP> problem;
    RCP< Belos::EpetraPrecOp > prec;//TODO just wrapper for MLPrec->LinearProblem
    RCP< Belos::StatusTestGenResNorm< ST, MV, OP > > convStatusTest;
	RCP< Belos::SolverManager<ST,MV,OP> > solver;

    /// Domain
    Box domain;
    const int* dom_lo;
    const int* dom_hi;

    int nx, ny, nz;

    /// mapping (i,j) or (i,j,k) -> idx
    std::map<int, int> IdxMap;
    /// mapping idx -> (i,j) or (i,j,k)
    std::map<int, int> CoordMap;

    /// mesh spacings
    std::vector<double> hr;

    /// communicator
    Epetra_MpiComm Comm_m;

    /// verbosity of the code
    bool verbose_m;

    /// flag specifying if hierarchy is reused
    bool isReusingHierarchy_m;

    /// flag specifying if whole preconditioner is reued
    bool isReusingPreconditioner_m;

    /// rectangular data distribution
    int my_slide_size[3], my_start_offset[3], my_end_offset[3];

    /** Performs a redistribution of the data with RCB to avoid idle processors.
     */
    void redistributeWithRCB();

    /** Computing the stencil in A and rescaling RHS with help of the geometry object.
     */
    inline void ComputeStencil();

    void printLoadBalanceStats();

    void extrapolateLHS();

    void InitializeDomainGeometry();

#if (BL_SPACEDIM == 2)

    inline bool isInside(int i, int j) {
        return (i > dom_lo[0] && i < dom_hi[0] && j > dom_lo[1] && j < dom_hi[1]);
    }

    inline int toCoordIdx(int i, int j) { return j*nx + i; }

    void ConstantInterpolation(int i, int j, 
                               double& W, double& E, double& S, double& N, 
                               double& C);

    void getBoundaryStencil(int x, int y, 
                            double& W, double& E, double& S, double& N, 
                            double& C);

    void getBoundaryStencil(int idx, 
                            double& W, double& E, double& S, double& N, 
                            double& C);

    void getNeighbours(int x, int y, 
                       double& W, double& E, double& S, double& N);

    void getNeighbours(int idx, 
                       double& W, double& E, double& S, double& N);

    void getCoord(int idx, int& i, int& j)
    {
        int ij = CoordMap[idx];
        j = ij / nx;
        i = ij % nx;  // (This is the same as idx - j*nx)
    }

    inline int getIdx(int i, int j) {
        if(isInside(i,j))
            return IdxMap[toCoordIdx(i,j)];
        else 
            return -1;
    }

#elif (BL_SPACEDIM == 3)

    inline bool isInside(int i, int j, int k) {
        return (i > dom_lo[0] && i < dom_hi[0] && j > dom_lo[1] && j < dom_hi[1] &&
                k > dom_lo[2] && k < dom_hi[2]);
    }

    inline int toCoordIdx(int i, int j, int k) { return k*nx*ny + j*nx + i; }

    void ConstantInterpolation(int i, int j, int k, 
                               double& W, double& E, double& S, double& N, 
                               double& F, double& B, double& C);

    void getBoundaryStencil(int x, int y, int z, 
                            double& W, double& E, double& S, double& N, 
                            double& F, double& B, double& C);

    void getBoundaryStencil(int idx, 
                            double& W, double& E, double& S, double& N, 
                            double& F, double& B, double& C);

    void getNeighbours(int x, int y, int z, 
                       double& W, double& E, double& S, double& N, double& F, double& B);

    void getNeighbours(int idx, 
                       double& W, double& E, double& S, double& N, double& F, double& B);

    void getCoord(int idx, int& i, int& j, int& k) 
    {
        int ijk = CoordMap[idx];
        k = ijk / (nx*ny);
        j = (ijk - k*nx*ny) / ny;
        i = (ijk - k*ny*nz) - j*nx;
    }

    inline int getIdx(int i, int j, int k) {
        if(isInside(i,j,k))
            return IdxMap[toCoordIdx(i,j,k)];
        else 
            return -1;
    }
#endif

protected:

    /** Setup the parameters for the SAAMG preconditioner.
     */
    inline void SetupMLList() {
        ML_Epetra::SetDefaults("SA", MLList_m);
        MLList_m.set("max levels", 8);
        MLList_m.set("increasing or decreasing", "increasing");

        // we use a V-cycle
        MLList_m.set("prec type", "MGV");

        // uncoupled aggregation is used (every processor aggregates
        // only local data)
        MLList_m.set("aggregation: type", "Uncoupled");

        // smoother related parameters
        MLList_m.set("smoother: type","Chebyshev");
        MLList_m.set("smoother: sweeps", 3);
        MLList_m.set("smoother: pre or post", "both");

        // on the coarsest level we solve with  Tim Davis' implementation of
        // Gilbert-Peierl's left-looking sparse partial pivoting algorithm,
        // with Eisenstat & Liu's symmetric pruning. Gilbert's version appears
        // as \c [L,U,P]=lu(A) in MATLAB. It doesn't exploit dense matrix
        // kernels, but it is the only sparse LU factorization algorithm known to be
        // asymptotically optimal, in the sense that it takes time proportional to the
        // number of floating-point operations.
        MLList_m.set("coarse: type", "Amesos-KLU");

        //XXX: or use Chebyshev coarse level solver
        // SEE PAPER FOR EVALUATION KLU vs. Chebyshev
        //MLList.set("coarse: sweeps", 10);
        //MLList.set("coarse: type", "Chebyshev");

        // turn on all output
        if(verbose_m)
            MLList_m.set("ML output", 101);
        else
            MLList_m.set("ML output", 10);

        // heuristic for max coarse size depending on number of processors
        int coarsest_size = std::max(Comm_m.NumProc() * 10, 1024);
        MLList_m.set("coarse: max size", coarsest_size);
    }


public:

    void setHr(const std::vector<double>& hr) { this->hr = hr; }
};

#endif
