
#ifndef AMREX_CGSOLVER_H_
#define AMREX_CGSOLVER_H_

#include <cmath>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LinOp.H>
#include <AMReX_ABecLaplacian.H>

namespace amrex {

class MultiGrid;

/**
        \brief A CGSolver solves the linear equation, L(phi)=rhs, for a LinOp L and
        MultiFabs rhs and phi, using the conjugate gradient algorithm, either
        vanilla, or preconditioned via a single V-cycle of multigrid.

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

        This algorithm implements solution of equations in a way that
        requires linearity of the operator.  In "residual correction form",
        equations are solved only for the "correction" to the initial guess
        that is required to solve the desired system.  In particular,
        boundary conditions are assumed to be satisfied after a single
        application of the linear operator (therefore, the correction is
        homogeneous at the boundary).  As a result, after putting the
        problem in residual correction form, the entire system CG system
        has homigeneous boundary conditions.  This solver therefore cannot
        incorporate fully nonlinear systems.

        In general, a LinOp has the ability to apply high-order
        interpolation methods for computing boundary values.  However, if
        the operator used in this CGSolver is non-symmetric, cg performance
        suffers tremendously.  As a result, BC's should be applied with
        interpolants of order less than or equal to 2.  We do not enforce
        this condition rigourously, since it would interfere with the user's
        definition of the operator.  If cg performance is bad however,
        setting the LinOp's maxorder=2 is a good place to start.

	CGSolver can fail to converge, or more likely diverge, if you
	run into a solvability problem (likely the accuracy cannot be
	achieved).  CGSolver returns error codes 0 means good, 1--7 means
	loss-of-accuracy, 8 means iteration count exceeded.

        Default settings:
        The user has access to the following options in the cg algorithm:

        maxiter(40) Maximum number of cg cycles allowed to solve the
        system before abandoning

        verbose(0) Verbosity (1-results, 2-progress, 3-detailed progress)

        use_mg_precond(false) Whether to use the V-cycle multigrid
        solver for the preconditioner system

	unstable_criterion(10) if norm of residual grows by more than
	this factor, it is taken as signal that you've run into a solvability
	problem.

        This class does NOT provide a copy constructor or assignment operator.
*/

class CGSolver
{
public:

    enum Solver { CG, BiCGStab, CABiCGStab, CABiCGStabQuad };

    /**
    * \brief The Constructor.
    *
    * \param _lp
    * \param _use_mg_precond
    * \param _lev
    */
    CGSolver (LinOp& _lp,
              bool   _use_mg_precond = false,
              int    _lev = 0);

    /**
    * \brief The Destructor.
    */
    ~CGSolver ();

    /**
    * \brief solve the system, Lp(solnL)=rhsL to relative err, tolerance
    * RETURNS AN INT!!!! indicating success or failure.
    * 0 means success
    * 1 means failed for loss of precision
    * 2 means iterations exceeded
    *
    * \param solnL
    * \param rhsL
    * \param eps_rel
    * \param eps_abs
    * \param bc_mode
    */
    int solve (MultiFab&       solnL,
	       const MultiFab& rhsL,
	       Real            eps_rel    = -1.0,
	       Real            eps_abs    = -1.0,
	       LinOp::BC_Mode  bc_mode    = LinOp::Inhomogeneous_BC);

    /**
    * \brief Set maximum allowed number of CG iterations.
    *
    * \param _maxiter
    */
    void setMaxIter (int _maxiter) { maxiter = _maxiter; }

    /**
    * \brief Get maximum allowed number of CG iterations.
    */
    int getMaxIter () const { return maxiter; }

    /**
    * \brief Set flag determining whether MG preconditioning is used.
    *
    * \param _use_mg_precond
    */
    void setUseMGPrecond (bool _use_mg_precond)
    {
        use_mg_precond = _use_mg_precond; set_mg_precond();
    }

    /**
    * \brief Get  flag determining whether MG preconditioning is used.
    */
    bool getUseMGPrecond () const { return use_mg_precond; }

    void setCGSolver (Solver s) { cg_solver = s; }
    Solver getCGSolver () const { return cg_solver; }


    /**
    * \brief Set the verbosity value.
    *
    * \param _verbose
    */
    void setVerbose (int _verbose) { verbose = _verbose; }

    /**
    * \brief Return the verbosity value.
    */
    int getVerbose () const { return verbose; }

protected:

    int solve_cg (MultiFab&       solnL,
		  const MultiFab& rhsL,
		  Real            eps_rel,
		  Real            eps_abs,
		  LinOp::BC_Mode  bc_mode);

    int solve_bicgstab (MultiFab&       solnL,
                        const MultiFab& rhsL,
                        Real            eps_rel,
                        Real            eps_abs,
                        LinOp::BC_Mode  bc_mode);

    int solve_cabicgstab (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          Real            eps_rel,
                          Real            eps_abs,
                          LinOp::BC_Mode  bc_mode);

    int solve_cabicgstab_quad (MultiFab&       solnL,
                               const MultiFab& rhsL,
                               Real            eps_rel,
                               Real            eps_abs,
                               LinOp::BC_Mode  bc_mode);

    int jbb_precond (MultiFab&       sol,
                     const MultiFab& rhs,
                     int             lev,
                     LinOp&          Lp);

    /**
    * \brief Construct work space, initialize parameters.
    */
    static void Initialize ();

    static void Finalize ();

    /**
    * \brief if  (use_mg_precond == 1) then define the MultiGrid * mg_precond
    */
    void set_mg_precond ();

private:

    static int    def_maxiter;
    static double def_unstable_criterion;
    static int    def_verbose;
    static Solver def_cg_solver;
    static bool   use_jbb_precond;        //!< Use JBB's new method as a preconditioner.
    static bool   use_jacobi_precond;     //!< Use Jacobi smoothing as a preconditioner.


    //! The data.
    LinOp&     Lp;             //!< Operator for linear system to be solved.
    MultiGrid* mg_precond;         //!< MultiGrid solver to be used as preconditioner
    int        maxiter;            //!< Current maximum number of allowed iterations.
    int        verbose;            //!< Current verbosity level.
    int        lev;                //!< Level of the linear operator to use
    bool       use_mg_precond;     //!< Use multigrid as a preconditioner.
    Solver     cg_solver;


    /**
    * \brief Disable copy constructor and assignment operator.
    *
    */
    CGSolver (const CGSolver&);
    CGSolver& operator= (const CGSolver&);
};

}

#endif /*_CGSOLVER_H_*/
