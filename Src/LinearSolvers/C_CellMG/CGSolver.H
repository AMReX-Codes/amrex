//BL_COPYRIGHT_NOTICE

#ifndef _CGSOLVER_H_
#define _CGSOLVER_H_

//
// $Id: CGSolver.H,v 1.8 2000-08-30 19:19:13 car Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cmath>
#else
#include <math.h>
#endif

#include <Array.H>
#include <Pointers.H>
#include <MultiFab.H>
#include <LinOp.H>

class MultiGrid;

//@Man:
/*@Memo:
  A CGSolver solves the linear equation, L(phi)=rhs, for a LinOp L and
  MultiFabs rhs and phi, using the conjugate gradient algorithm, either
  vanilla, or preconditioned via a single V-cycle of multigrid.
*/        
/*@Doc:
        A CGSolver object solves the linear equation, L(phi)=rhs for a LinOp
        L, and MultiFabs phi and rhs.  A CGSolver is constructed with a
        fully initialized 2D or 3D LinOp, and responds to "solve" requests of
        various signatures.  The solve request (implicitly) includes a flag
        as to whether the system is to be solved with homogeneous boundary
        conditions or no.  If homogeneous BCs are requested, all boundary
        information within the LinOp is used, except that the values of
        boundary FabSets are zeroed.  The algorithm follows closely that
        described of p15 of the SIAM book, "Templates for the Solution of
        Linear Systems".  Before solving a system, a small number of control
        parameters may be modified (maximum number of allowed iterations,
        reporting verbosity, etc).  Also, either on construction or via
        member access functions, the user may choose to use a V-cycle
        1-step MultiGrid preconditioner within the CG loop by setting
        the appropriate bool flag (see members/ctrs below).  The MG
        preconditioner used is just an instantiation of a MultiGrid class
        object (the MultiGrid class is documented separately).

        Implementation Notes:
        \begin{itemize}
        \item This algorithm implements solution of equations in a way that
        requires linearity of the operator.  In "residual correction form",
        equations are solved only for the "correction" to the initial guess
        that is required to solve the desired system.  In particular,
        boundary conditions are assumed to be satisfied after a single
        application of the linear operator (therefore, the correction is
        homogeneous at the boundary).  As a result, after putting the
        problem in residual correction form, the entire system CG system
        has homigeneous boundary conditions.  This solver therefore cannot
        incorporate fully nonlinear systems.

        \item In general, a LinOp has the ability to apply high-order
        interpolation methods for computing boundary values.  However, if
        the operator used in this CGSolver is non-symmetric, cg performance
        suffers tremendously.  As a result, BC's should be applied with
        interpolants of order less than or equal to 2.  We do not enforce
        this condition rigourosly, since it would interfere with the user's
        definition of the operator.  If cg performance is bad however,
        setting the LinOp's maxorder=2 is a good place to start.

	\item CGSolver can fail to converge, or more likely diverge, if you
	run into a solvability problem (likely the accuracy cannot be 
	achieved).  CGSolver can be put into a different mode by using 
	setExpert, where instead of bombing out, solve() returns an error
	code.  0 means good, 1 means loss-of-accuracy, 2 means iteration
	count exceeded.

        \end{itemize}
        
        Default settings:
        The user has access to the following options in the cg algorithm:
        \begin{itemize}
        \item maxiter(40) Maximum number of cg cycles allowed to solve the
        system before abandoning
        \item verbose(0) Verbosity (1-results, 2-progress, 3-detailed progress)
        \item use\_mg\_precond(false) Whether to use the V-cycle multigrid
        solver for the preconditioner system
	\item unstable\_criterion(10) if norm of residual grows by more than 
	this factor, it is taken as signal that you've run into a solvability
	problem.
        \end{itemize}
        
        This class does NOT provide a copy constructor or assignment operator.
        
*/

class CGSolver
{
public:
    enum Solver { CG, BiCGStab, CG_Alt, Automatic};

  //
    //@ManDoc: constructor
    //
    CGSolver (LinOp& _Lp,
              bool   _use_mg_precond = false,
              int    _lev = 0);
    //
    //@ManDoc: destructor
    //
    virtual ~CGSolver ();
    //
    //@ManDoc: solve the system, Lp(solnL)=rhsL to relative err, tolerance
    // RETURNS AN INT!!!! indicating success or failure.  
    // 0 means success
    // 1 means failed for loss of precision
    // 2 means iterations exceeded
    //
    int solve (MultiFab&       solnL,
	       const MultiFab& rhsL,
	       Real            eps_rel = -1.0,
	       Real            eps_abs = -1.0,
	       LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC,
	       Solver cg_solver_ = Automatic);
    //
    //@ManDoc: set maximum allowed number of CG iterations
    //
    void setMaxIter (int _maxiter);
    //
    //@ManDoc: get maximum allowed number of CG iterations
    //
    int getMaxIter () const;
    //
    //@ManDoc: set expert mode
    void setExpert( bool flag ) { isExpert = flag; }
    //
    //@ManDoc: get expert flag
    bool getExpert() { return isExpert; }
    //
    //@ManDoc: set flag determining whether MG preconditioning is used
    //
    void setUseMGPrecond (bool _use_mg_precond);
    //
    //@ManDoc: get  flag determining whether MG preconditioning is used
    //
    bool getUseMGPrecond () const;
    //
    //@ManDoc: set the verbosity value
    //
    void setVerbose (int _verbose);
    //
    //@ManDoc: return the verbosity value
    //
    int getVerbose ();

    void setCGSolver(Solver cg_solver_);
    Solver getCGSolver() const;

protected:

    int solve_00 (MultiFab&       solnL,
	       const MultiFab& rhsL,
	       Real            eps_rel = -1.0,
	       Real            eps_abs = -1.0,
	       LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);
    int solve_cg (MultiFab&       solnL,
		  const MultiFab& rhsL,
		  Real            eps_rel = -1.0,
		  Real            eps_abs = -1.0,
		  LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);
    int solve_bicgstab (MultiFab&       solnL,
		  const MultiFab& rhsL,
		  Real            eps_rel = -1.0,
		  Real            eps_abs = -1.0,
		  LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);
    //
    //@ManDoc: construct work space, initialize parameters
    //
    static void initialize ();
    //
    //@ManDoc: if  (use\_mg\_precond == 1) then define the MultiGrid * mg\_precond
    //
    virtual void set_mg_precond ();
    //
    //@ManDoc: compute p = z  +  beta p in the CG algorithm
    //
    void advance (MultiFab&       p,
                  Real            beta,
                  const MultiFab& z);
    //
    //@ManDoc: compute x =+ alpha p  and  r -= alpha w in the CG algorithm
    //
    void update (MultiFab&       sol,
                 Real            alpha,
                 MultiFab&       r,
                 const MultiFab& p,
                 const MultiFab& w);
    //
    //@ManDoc: compute w = A.p, and return Transpose(p).w in the CG algorithm
    //
    Real axp (MultiFab&      w,
              MultiFab&      p,
              LinOp::BC_Mode bc_mode);
    //
    //@ManDoc: Initialized yet?
    //
    static int initialized;
    //
    //@ManDoc: Default maximum number of allowed iterations.
    //
    static int def_maxiter;
    //
    //@ManDoc: default unstable_criterion used to test for loss of accuracy in
    //
    static double def_unstable_criterion;
    //
    //@ManDoc: flag determining action when error conditions are hit
    //
    bool isExpert;
    //
    //@ManDoc: Default verbosity flag.
    //
    static int def_verbose;
    static Solver def_cg_solver;
    Solver cg_solver;
    //
    // The data.
    //
    LinOp&     Lp;             // Operator for linear system to be solved.
    MultiGrid* mg_precond;     // MultiGrid solver to be used as preconditioner
    int        maxiter;        // Current maximum number of allowed iterations.
    int        verbose;        // Current verbosity level.
    int        lev;            // Level of the linear operator to use
    bool       use_mg_precond; // Use multigrid as a preconditioner.
    
private:
    //
    // Disable copy constructor and assignment operator.
    //
    CGSolver (const CGSolver&);
    CGSolver& operator= (const CGSolver&);
};

inline
void
CGSolver::setMaxIter (int _maxiter)
{
    maxiter = _maxiter;
}

inline
int
CGSolver::getMaxIter () const
{
    return maxiter;
}

inline
void
CGSolver::setUseMGPrecond (bool _use_mg_precond)
{
    use_mg_precond = _use_mg_precond;
    set_mg_precond();
}

inline
bool
CGSolver::getUseMGPrecond () const
{
    return use_mg_precond;
}

inline
void
CGSolver::setVerbose (int _verbose)
{
    verbose = _verbose;
}

inline
int
CGSolver::getVerbose ()
{
    return verbose;
}

inline
void
CGSolver::setCGSolver(Solver cg_solver_)
{
  cg_solver = cg_solver_;
}

inline
CGSolver::Solver
CGSolver::getCGSolver() const
{
  return cg_solver;
}


#endif /*_CGSOLVER_H_*/
