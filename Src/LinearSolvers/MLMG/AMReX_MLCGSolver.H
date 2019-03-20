
#ifndef AMREX_MLCGSOLVER_H_
#define AMREX_MLCGSOLVER_H_

#include <cmath>

#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLLinOp.H>

namespace amrex {

class MLMG;

class MLCGSolver
{
public:

    enum struct Type { BiCGStab, CG };

    MLCGSolver (MLMG* a_mlmg, MLLinOp& _lp, Type _typ = Type::BiCGStab);
    ~MLCGSolver ();

    MLCGSolver (const MLCGSolver& rhs) = delete;
    MLCGSolver& operator= (const MLCGSolver& rhs) = delete;

    void setSolver (Type _typ) noexcept { solver_type = _typ; }

    /**
    * solve the system, Lp(solnL)=rhsL to relative err, tolerance
    * RETURNS AN INT!!!! indicating success or failure.
    * 0 means success
    * 1 means failed for loss of precision
    * 2 means iterations exceeded
    */
    int solve (MultiFab&       solnL,
	       const MultiFab& rhsL,
	       Real            eps_rel,
	       Real            eps_abs);

    void setVerbose (int _verbose) { verbose = _verbose; }
    int getVerbose () const { return verbose; }

    void setMaxIter (int _maxiter) { maxiter = _maxiter; }
    int getMaxIter () const { return maxiter; }

    Real dotxy (const MultiFab& r, const MultiFab& z, bool local = false);
    Real norm_inf (const MultiFab& res, bool local = false);
    int solve_bicgstab (MultiFab&       solnL,
                        const MultiFab& rhsL,
                        Real            eps_rel,
                        Real            eps_abs);
    int solve_cg (MultiFab&       solnL,
                  const MultiFab& rhsL,
                  Real            eps_rel,
                  Real            eps_abs);

private:

    MLMG* mlmg;
    MLLinOp& Lp;
    Type solver_type;
    const int amrlev;
    const int mglev;
    int    verbose   = 0;
    int    maxiter   = 100;
};

}

#endif /*_CGSOLVER_H_*/
