//BL_COPYRIGHT_NOTICE

#ifndef _MCCGSOLVER_H_
#define _MCCGSOLVER_H_

//
// $Id: MCCGSolver.H,v 1.4 1999-01-04 21:31:48 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cmath>
#else
#include <math.h>
#endif

#include <Boolean.H>
#include <Array.H>
#include <Pointers.H>

#include <MultiFab.H>

#include <MCLinOp.H>
#include <MCMultiGrid.H>

//@Man:
/*@Memo:
  A CGColver solves the linear equation, L(phi)=rhs, for a MCLinOp L and
  MultiFabs rhs and phi, using the conjugate gradient algorithm, either
  vanilla, or preconditioned via a single V-cycle of multigrid.
*/	
/*@Doc:
        A MCCGSolver object solves the linear equation, L(phi)=rhs for a MCLinOp
	L, and MultiFabs phi and rhs.  A MCCGSolver is constructed with a
	fully initialized 2D or 3D MCLinOp, and responds to "solve" requests of
	various signatures.  The solve request (implicitly) includes a flag
	as to whether the system is to be solved with homogeneous boundary
	conditions or no.  If homogeneous BCs are requested, all boundary
	information within the MCLinOp is used, except that the values of
	boundary FabSets are zeroed.  The algorithm follows closely that
	described of p15 of the SIAM book, "Templates for the Solution of
	Linear Systems".  Before solving a system, a small number of control
	parameters may be modified (maximum number of allowed iterations,
	reporting verbosity, etc).  Also, either on construction or via
	member access functions, the user may choose to use a V-cycle
	1-step MCMultiGrid preconditioner within the CG loop by setting
	the appropriate bool flag (see members/ctrs below).  The MG
	preconditioner used is just an instantiation of a MCMultiGrid class
	object (the MCMultiGrid class is documented separately).

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

	\item In general, a MCLinOp has the ability to apply high-order
	interpolation methods for computing boundary values.  However, if
	the operator used in this MCCGSolver is non-symmetric, cg performance
	suffers tremendously.  As a result, BC's should be applied with
	interpolants of order less than or equal to 2.  We do not enforce
	this condition rigourosly, since it would interfere with the user's
	definition of the operator.  If cg performance is bad however,
	setting the MCLinOp's maxorder=2 is a good place to start.

        \end{itemize}
	
	Default settings:
	The user has access to the following options in the cg algorithm:
	\begin{itemize}
	\item maxiter(40) Maximum number of cg cycles allowed to solve the
	system before abandoning
	\item verbose(0) Verbosity (1-results, 2-progress, 3-detailed progress)
	\item use_mg_precond(false) Whether to use the V-cycle multigrid
	solver for the preconditioner system
	\end{itemize}
*/

class MCCGSolver
{
public:
    //
    //@ManDoc: constructor
    //
    MCCGSolver (MCLinOp& _Lp,
		bool     _use_mg_precond = false,
		int      _lev=0);
    //
    //@ManDoc: destructor
    //
    virtual ~MCCGSolver ();
    //
    //@ManDoc: solve the system, Lp(solnL)=rhsL to relative err, tolerance
    //
    virtual void solve (MultiFab&       solnL,
			const MultiFab& rhsL,
			Real            eps_rel = -1.0,
			Real            eps_abs = -1.0,
			MCBC_Mode       bc_mode = MCInhomogeneous_BC);
    //
    //@ManDoc: set maximum allowed number of CG iterations
    //
    void setMaxIter (int _maxiter);
    //
    //@ManDoc: get maximum allowed number of CG iterations
    //
    int getMaxIter () const;
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

protected:
    //
    //@ManDoc: construct work space, initialize parameters
    //
    static void initialize ();
    //
    //@ManDoc: if (use_mg_precond == 1) then define the MCMultiGrid * mg_precond
    //
    void set_mg_precond ();
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
    Real axp (MultiFab& w,
	      MultiFab& p,
	      MCBC_Mode bc_mode);
    //
    //@ManDoc: compute a suitable norm of the residual, used to check stopping criteria
    //
    Real norm (const MultiFab& res);
    //
    //@ManDoc: MCMultiGrid solver to be used as preconditioner
    //
    MCMultiGrid* mg_precond;
    //
    //@ManDoc: flag, whether initialized yet
    //
    static int initialized;
    //
    //@ManDoc: default maximum number of allowed iterations, verbosity
    //
    static int def_maxiter, def_verbose;
    //
    //@ManDoc: flag: use multigrid as a preconditioner
    //
    bool use_mg_precond;
    //
    //@ManDoc: current maximum number of allowed iterations, verbosity
    //
    int maxiter, verbose;
    //
    //@ManDoc: reference to operator for linear system to be solved
    //
    MCLinOp &Lp;
    //
    //@ManDoc: integer to choose which "level" of the linear operator to use
    //
    int lev;
    
private:
    //
    //@ManDoc: disable copy constructor, = operator
    //
    MCCGSolver (const MCCGSolver&);
    MCCGSolver& operator= (const MCCGSolver&);
};

inline
void
MCCGSolver::setMaxIter (int _maxiter)
{
    maxiter = _maxiter;
}

inline
int
MCCGSolver::getMaxIter () const
{
    return maxiter;
}

inline
void
MCCGSolver::setUseMGPrecond (bool _use_mg_precond)
{
    use_mg_precond = _use_mg_precond;
    set_mg_precond();
}

inline
bool
MCCGSolver::getUseMGPrecond () const
{
    return use_mg_precond;
}

inline
void
MCCGSolver::setVerbose (int _verbose)
{
    verbose = _verbose;
}

inline
int
MCCGSolver::getVerbose ()
{
    return verbose;
}

#endif /*_MCCGSOLVER_H_*/
