
#ifndef AMREX_MLCGSOLVER_H_
#define AMREX_MLCGSOLVER_H_

#include <cmath>

#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLLinOp.H>

namespace amrex {

class MLCGSolver
{
public:

    enum class Solver { BiCGStab, CG };

    MLCGSolver (MLLinOp& _lp, Solver _solver = Solver::BiCGStab);
    ~MLCGSolver ();

    MLCGSolver (const MLCGSolver& rhs) = delete;
    MLCGSolver& operator= (const MLCGSolver& rhs) = delete;

    //
    // solve the system, Lp(solnL)=rhsL to relative err, tolerance
    // RETURNS AN INT!!!! indicating success or failure.  
    // 0 means success
    // 1 means failed for loss of precision
    // 2 means iterations exceeded
    //
    int solve (MultiFab&       solnL,
	       const MultiFab& rhsL,
	       Real            eps_rel,
	       Real            eps_abs);

    void setVerbose (int _verbose) { verbose = _verbose; }
    int getVerbose () const { return verbose; }

    void setMaxIter (int _maxiter) { maxiter = _maxiter; }
    int getMaxIter () const { return maxiter; }

    void setUnstableCriterion (int i) { unstable_criterion = i; }
    int getUnstableCriterion () const { return unstable_criterion; }

protected:

    int solve_cg (MultiFab&       solnL,
		  const MultiFab& rhsL,
		  Real            eps_rel,
		  Real            eps_abs);

    int solve_bicgstab (MultiFab&       solnL,
                        const MultiFab& rhsL,
                        Real            eps_rel,
                        Real            eps_abs);

    int solve_cabicgstab (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          Real            eps_rel,
                          Real            eps_abs);

    int solve_cabicgstab_quad (MultiFab&       solnL,
                               const MultiFab& rhsL,
                               Real            eps_rel,
                               Real            eps_abs);

private:

    MLLinOp& Lp;
    const int amrlev;
    const int mglev;
    Solver cg_solver;
    int    verbose   = 0;
    int    maxiter   = 100;
    int    unstable_criterion = 10;
};

}

#endif /*_CGSOLVER_H_*/
