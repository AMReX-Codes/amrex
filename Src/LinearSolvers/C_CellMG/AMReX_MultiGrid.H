
#ifndef AMREX_MULTIGRID_H_
#define AMREX_MULTIGRID_H_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BndryData.H>
#include <AMReX_LinOp.H>
#include <AMReX_CGSolver.H>

#include <algorithm>

namespace amrex {

/**
  \brief A MultiGrid solves the linear equation, ``L(phi)=rhs``, for a LinOp L and
  MultiFabs rhs and phi using a V-type cycle of the MultiGrid algorithm

  A MultiGrid object solves the linear equation, ``L(phi)=rhs`` for a LinOp
  L, and MultiFabs phi and rhs.  A MultiGrid is constructed with a
  fully initialized 2D or 3D LinOp, and responds to "solve" requests of
  various signatures, ultimately performing a recursive "relax"
  operation over a hierachy of grid levels.  The LinOp therefore must
  employ "levels" of application, as well as be able to provide an
  implementation of the Gauss-Seidel red-black iterations on all levels.
  At the coarsest level, the user has the option of applying the
  GSRB smoother a set number of iterations.  Optionally, a Conjugate
  Gradient solver, CGSolver, may be used to solve the coarsest system.

  If the user chooses to use the conjugate gradient bottom solve,
  the absolute and relative tolerances of this solve are independently
  settable distinct from the tolerances of the mg solve.  However,
  it rarely makes sense to stop cg after a fixed number of iterations,
  so that capability was omited, wrt the mg options.  As a final
  option, the user may choose to follow the cg solve with a fixed
  number of relaxation passes (the relaxer within the class LinOp).

  The solve request (implicitly) includes a flag as to whether the
  system is to be solved with homogeneous boundary conditions or no.
  If homogeneous BCs are requested, all boundary information within
  the LinOp is used, except that the values of boundary FabSets are
  zeroed.

  Implementation Note:
  This algorithm implements solution of equations in a way that
  requires linearity of the operator.  In "residual correction form",
  equations are solved only for the "correction" to the initial guess
  that is required to solve the desired system.  In particular,
  boundary conditions are assumed to be satisfied after a single
  application of the linear operator (therefore, the correction is
  homogeneous at the boundary).  As a result, after putting the
  problem in residual correction form, the entire system MG hierarchy
  has homigeneous boundary conditions (thus avoiding the need to
  interpolate BC values, and simplifying the logic of the coarse-fine
  transfer operations).  This solver therefore cannot incorporate
  fully nonlinear systems.

  Default settings:
  There are a number of options in the multigrid algorithm details.
  In addition to changing the actual smoothers employed, the user
  has access to the following parameters (defaults shown in parentheses):

   nu_1(2) Number of passes of the pre-smoother
   nu_2(2) Number of passes of the post-smoother
   nu_0(1) Number of passes of the coarse-grid smoother per cycle
   nu_f(8) Number of passes of the bottom smoother (if not using
           the conjugate-gradient bottom solve)
   maxiter(40)  Maximum number of full multigrid cycles allowed to
                solve the system value is ignored if < 0)
   verbose(0)   Verbosity (1-results, 2-progress)
   usecg(1)     Whether to use the conjugate-gradient solver for the
                coarsest (bottom) solve of the multigrid hierarchy
   atol_b(-1.0) Absolute error tolerance (<0 => ignored) for cg
   rtol_b(.01)  Relative error tolerance (<0 => ignored) for cg
   nu_b(0)      Number of passes of the bottom smoother taken
                AFTER the cg bottom solve (value ignored if <= 0)
   numLevelsMAX(1024) maximum number of mg levels

  This class does NOT provide a copy constructor or assignment operator.
*/

class MultiGrid
{
public:

    /**
    * \brief constructor
    *
    * \param _lp
    */
    MultiGrid (LinOp& _lp);

    /**
    * \brief destructor
    */
    ~MultiGrid ();

     /**
     * \brief solve the system to relative err eps_rel, absolute err eps_abs
     *
     * \param solution
     * \param _rhs
     * \param eps_rel
     * \param eps_abs
     * \param bc_mode
     */
     void solve (MultiFab&       solution,
                 const MultiFab& _rhs,
                 Real            eps_rel = -1.0,
                 Real            eps_abs = -1.0,
                 LinOp::BC_Mode  bc_mode=LinOp::Inhomogeneous_BC);

    /**
    * \brief return the linear operator
    */
    LinOp& linOp () { return Lp; }

    /**
    * \brief set the maximum permitted multigrid iterations
    *
    * \param _maxiter
    */
    void setMaxIter (int _maxiter) { maxiter = _maxiter; }

    /**
    * \brief set the maximum number of permitted multigrid iterations --
    * if this flag is set to 1 then multigrid will do exactly this number,
    * regardless of convergence
    *
    * \param _do_fixed_number_of_iters
    */
    void setFixedIter (int _do_fixed_number_of_iters) { do_fixed_number_of_iters = _do_fixed_number_of_iters;}

    /**
    * \brief return the maximum permitted multigrid iterations
    */
    int getMaxIter () const { return maxiter; }

    /**
    * \brief return the number of multigrid iterations
    */
    int getNumIter () const;

    /**
    * \brief set the flag for whether to use CGSolver at coarsest level
    *
    * \param _usecg
    */
    void setUseCG (int _usecg) { usecg = _usecg; }

    /**
    * \brief return the flag for whether to use CGSolver at coarsest level
    */
    int getUseCG () const { return usecg; }

    /**
    * \brief set/return the number of multigrid levels
    *
    * \param _numlevels
    */
    int getNumLevels (int _numlevels);

    /**
    * \brief return the number of multigrid levels
    */
    int getNumLevels () const { return numlevels; }

    /**
    * \brief set the verbosity value
    *
    * \param _verbose
    */
    void setVerbose (int _verbose) { verbose = _verbose; }

    /**
    * \brief return the verbosity value
    */
    int getVerbose () { return verbose; }

    /**
    * \brief set the number of passes of the pre-smoother
    *
    * \param pre_smooth
    */
    void set_preSmooth (int pre_smooth) { nu_1 = pre_smooth; }

    /**
    * \brief set the number of passes of the post-smoother
    *
    * \param post_smooth
    */
    void set_postSmooth (int post_smooth) { nu_2 = post_smooth; }

    /**
    * \brief set the number of passes of the coarse-grid smoother/mg iteration
    *
    * \param cnt_relax
    */
    void set_cntRelax (int cnt_relax) { nu_0 = cnt_relax; }

    /**
    * \brief set the number of passes of the bottom mg relaxer
    *
    * \param final_smooth
    */
    void set_finalSmooth (int final_smooth) { nu_f = final_smooth; }

    /**
    * \brief Return the number of pre-smooth iterations at the level
    */
    int preSmooth () const { return nu_1; }

    /**
    * \brief Return the number of post-smooth iterations at the level
    */
    int postSmooth () const { return nu_2; }

    /**
    * \brief Return the number of level relaxations (not implemented)
    */
    int cntRelax () const { return nu_0; }

    /**
    * \brief Return the number of smoothing passes at bottom of MG (if no cg)
    */
    int finalSmooth () const { return nu_f; }

    /**
    * \brief set the maximum permitted absolute tolerance (\<0 bypasses test)
    *
    * \param atol
    */
    void set_atol_b (Real atol) { atol_b = atol; }

    /**
    * \brief get the maximum permitted absolute tolerance
    */
    Real get_atol_b () const { return atol_b; }

    /**
    * \brief set the maximum permitted relative tolerance (\<0 bypasses test)
    *
    * \param rtol
    */
    void set_rtol_b (Real rtol) { rtol_b = rtol; }

    /**
    * \brief get the maximum permitted relative tolerance
    */
    Real get_rtol_b () const { return rtol_b; }

    /**
    * \brief set the number of post-cg relax passes
    *
    * \param _nu_b
    */
    void set_nu_b (int _nu_b) { nu_b = _nu_b; }

    /**
    * \brief set the number of post-cg relax passes
    */
    int get_nu_b () const { return nu_b; }

    void set_maxiter_b (int n) { maxiter_b = n; }

    /**
    * \brief get the maximum permitted relative tolerance
    */
    int  get_maxiter_b () const { return maxiter_b; }

protected:

    /**
    * \brief Solve the linear system to relative and absolute tolerance
    *
    * \param _sol
    * \param _eps_rel
    * \param _eps_abs
    * \param bc_mode
    * \param bnorm
    * \param resnorm0
    */
    int solve_ (MultiFab&      _sol,
                Real           _eps_rel,
                Real           _eps_abs,
                LinOp::BC_Mode bc_mode,
                Real           bnorm,
                Real           resnorm0);

    /**
    * \brief Make space, set switches for new solution level
    *
    * \param level
    */
    void prepareForLevel (int level);

    /**
    * \brief Compute the number of multigrid levels, assuming ratio=2
    */
    int numLevels () const;

    /**
    * \brief Return scalar estimate of error
    *
    * \param level
    * \param bc_mode
    * \param local
    */
    Real errorEstimate (int            level,
                        LinOp::BC_Mode bc_mode,
                        bool           local = false);

    /**
    * \brief Transfer MultiFab from fine to coarse level
    *
    * \param c
    * \param f
    */
    void average (MultiFab&       c,
                  const MultiFab& f);

    /**
    * \brief Transfer MultiFab from coarse to fine level
    *
    * \param f
    * \param c
    */
    void interpolate (MultiFab&       f,
                      const MultiFab& c);

    /**
    * \brief Perform a MG V-cycle
    *
    * \param solL
    * \param rhsL
    * \param level
    * \param eps_rel
    * \param eps_abs
    * \param bc_mode
    * \param cg_time
    */
    void relax (MultiFab&      solL,
                MultiFab&      rhsL,
                int            level,
                Real           eps_rel,
                Real           eps_abs,
                LinOp::BC_Mode bc_mode,
                Real&          cg_time);

    /**
    * \brief Perform relaxation at bottom of V-cycle
    *
    * \param solL
    * \param rhsL
    * \param level
    * \param eps_rel
    * \param eps_abs
    * \param bc_mode
    * \param local_usecg
    * \param cg_time
    */
    void coarsestSmooth (MultiFab&      solL,
                         MultiFab&      rhsL,
                         int            level,
                         Real           eps_rel,
                         Real           eps_abs,
                         LinOp::BC_Mode bc_mode,
                         int            local_usecg,
                         Real&          cg_time);
private:

    //! default flag, whether to use CG at bottom of MG cycle
    static int def_usecg;

    //! \brief set flags, etc
    static void Initialize ();

    static void Finalize ();

    //! default number of level, pre-, post- and bottom smooths
    static int def_nu_0, def_nu_1, def_nu_2, def_nu_f;

    //! default number of post-cg relax passes
    static int def_nu_b;

    //! default maximum number of complete MG cycle iterations
    static int def_maxiter;

    //! default convergence setting
    static int use_Anorm_for_convergence;

    //! default verbosity
    static int def_verbose;

    //! default relative, absolute tolerance for cg solve
    static Real def_rtol_b, def_atol_b;
    static int def_maxiter_b;

    //! default maximum number of mg levels
    static int def_numLevelsMAX;

    //! default flag controlling the behaviour when the cg solver returns indicating unstable
    static int def_smooth_on_cg_unstable;

    //! verbosity
    int verbose;

    //! Number of MG levels
    int numlevels;

    //! current maximum number of allowed iterations
    int maxiter;

    //! whether to use a fixed number of iterations or continue until convergence
    int do_fixed_number_of_iters;

    //! current number of level, pre-, post- and bottom smooths
    int nu_0, nu_1, nu_2, nu_f;

    //! current number of post-cg relax passes
    int nu_b;

    //! current flag, whether to use CG at bottom of MG cycle
    int usecg;

    //! relative, absolute tolerance, maxiter for cg bottom solver.
    Real rtol_b, atol_b;
    int maxiter_b;

    //! maximum number of mg levels
    int numLevelsMAX;

    //! flag controlling the behaviour when the cg solver returns indicating unstable
    int smooth_on_cg_unstable;

    //! internal temp data to store initial guess of solution
    MultiFab* initialsolution;

    //! internal temp data
    Vector< MultiFab* > res;

    //! internal temp data
    Vector< MultiFab* > rhs;

    //! internal temp data
    Vector< MultiFab* > cor;

    //! internal reference to linear operator
    LinOp &Lp;

    //! \brief Disallow copy constructor, assignment operator
    MultiGrid (const MultiGrid&);
    MultiGrid& operator= (const MultiGrid&);
};

}

#endif /*_MULTIGRID_H_*/
