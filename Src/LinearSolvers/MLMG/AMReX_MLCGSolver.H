
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

    MLCGSolver (MLLinOp& _lp);
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

private:

    MLLinOp& Lp;
    const int amrlev;
    const int mglev;
    int    verbose   = 0;
    int    maxiter   = 100;

    Real dotxy (const MultiFab& r, const MultiFab& z, bool local = false);
    Real norm_inf (const MultiFab& res, bool local = false);
};

}

#endif /*_CGSOLVER_H_*/
