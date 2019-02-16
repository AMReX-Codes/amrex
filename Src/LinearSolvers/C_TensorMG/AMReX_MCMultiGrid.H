
#ifndef AMREX_MCMULTIGRID_H_
#define AMREX_MCMULTIGRID_H_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MCLinOp.H>

namespace amrex {

/**
  /brief A MCMultiGrid solves the linear equation, L(phi)=rhs, for a MCLinOp L and
  MultiFabs rhs and phi using a V-type cycle of the MCMultiGrid algorithm

  A MCMultiGrid object solves the linear equation, L(phi)=rhs for a MCLinOp
  L, and MultiFabs phi and rhs.  A MCMultiGrid is constructed with a
  fully initialized 2D or 3D MCLinOp, and responds to "solve" requests of
  various signatures, ultimately performing a recursive "relax"
  operation over a hierachy of grid levels.  The MCLinOp therefore must
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
  number of relaxation passes (the relaxer within the class MCLinOp).

  The solve request (implicitly) includes a flag as to whether the
  system is to be solved with homogeneous boundary conditions or no.
  If homogeneous BCs are requested, all boundary information within
  the MCLinOp is used, except that the values of boundary FabSets are
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
  maxiter(40) Maximum number of full multigrid cycles allowed to
              solve the system
  numiter(-1) Number of full multigrid cycles to perform
              (should be less than maxiter for fixed number of MG cycles;
              value is ignored if < 0)
  verbose(0) Verbosity (1-results, 2-progress)
  usecg(1) Whether to use the conjugate-gradient solver for the
           coarsest (bottom) solve of the multigrid hierarchy
  atol_b(-1.0) Absolute error tolerance (<0 => ignored) for cg
  rtol_b(.01)  Relative error tolerance (<0 => ignored) for cg
  nu_b(0)      Number of passes of the bottom smoother taken
               AFTER} the cg bottom solve (value ignored if <= 0)
  numLevelsMAX(1024) maximum number of mg levels
  */

class MCMultiGrid
{
public:

    /**
    * \brief Constructor.
    *
    * \param &_lp
    */
    MCMultiGrid (MCLinOp &_lp);

    /**
    * \brief Destructor.
    */
    virtual ~MCMultiGrid ();

    /**
    * \brief Solve the system to relative err eps_rel, absolute err eps_abs.
    *
    * \param solution
    * \param _rhs
    * \param eps_rel
    * \param eps_abs
    * \param bc_mode
    */
    virtual void solve (MultiFab&       solution,
			const MultiFab& _rhs,
			Real            eps_rel = -1.0,
			Real            eps_abs = -1.0,
			MCBC_Mode       bc_mode = MCInhomogeneous_BC);

    /**
    * \brief Return the linear operator.
    */
    MCLinOp& linOp ();

    /**
    * \brief Set the maximum permitted multigrid iterations.
    *
    * \param _maxiter
    */
    void setMaxIter (int _maxiter);

    /**
    * \brief Return the maximum permitted multigrid iterations.
    */
    int getMaxIter () const;

    /**
    * \brief Set the number of multigrid iterations to perform.
    *
    * \param _numiter
    */
    void setNumIter (int _numiter);

    /**
    * \brief Return the number of multigrid iterations.
    */
    int getNumIter () const;

    /**
    * \brief Set the flag for whether to use CGSolver at coarsest level.
    *
    * \param _usecg
    */
    void setUseCG (int _usecg);

    /**
    * \brief Return the flag for whether to use CGSolver at coarsest level.
    */
    int getUseCG () const;

    /**
    * \brief Set/return the number of multigrid levels.
    *
    * \param _numlevels
    */
    int getNumLevels (int _numlevels);

    /**
    * \brief Return the number of multigrid levels.
    */
    int getNumLevels () const;

    /**
    * \brief Return the number of components in solution.
    */
    int getNumComponents () const;

    /**
    * \brief Set the verbosity value.
    *
    * \param _verbose
    */
    void setVerbose (int _verbose);

    /**
    * \brief Return the verbosity value.
    */
    int getVerbose ();

    /**
    * \brief Set the number of passes of the pre-smoother.
    *
    * \param pre_smooth
    */
    void set_preSmooth (int pre_smooth);

    /**
    * \brief Set the number of passes of the post-smoother.
    *
    * \param post_smooth
    */
    void set_postSmooth (int post_smooth);

    /**
    * \brief Set the number of passes of the coarse-grid smoother/mg iteration.
    *
    * \param cnt_relax
    */
    void set_cntRelax (int cnt_relax);

    /**
    * \brief Set the number of passes of the bottom mg relaxer.
    *
    * \param final_smooth
    */
    void set_finalSmooth (int final_smooth);

    /**
    * \brief Return the number of pre-smooth iterations at the level.
    */
    int preSmooth () const;

    /**
    * \brief Return the number of post-smooth iterations at the level.
    */
    int postSmooth () const;

    /**
    * \brief Return the number of level relaxations  (not implemented).
    */
    int cntRelax () const;

    /**
    * \brief Return the number of smoothing passes at bottom of MG (if no cg).
    */
    int finalSmooth () const;

    /**
    * \brief Set the maximum permitted absolute tolerance (\<0 bypasses test).
    *
    * \param atol
    */
    void set_atol_b (Real atol);

    /**
    * \brief Get the maximum permitted absolute tolerance.
    */
    Real get_atol_b () const;

    /**
    * \brief Set the maximum permitted relative tolerance (\<0 bypasses test).
    *
    * \param rtol
    */
    void set_rtol_b (Real rtol);

    /**
    * \brief Get the maximum permitted relative tolerance.
    */
    Real get_rtol_b () const;

    /**
    * \brief Set the number of post-cg relax passes.
    *
    * \param _nu_b
    */
    void set_nu_b (int _nu_b);

    /**
    * \brief Set the number of post-cg relax passes.
    */
    int get_nu_b () const;

protected:

    /**
    * \brief Solve the linear system to relative and absolute tolerance
    *
    * \param _sol
    * \param _eps_rel
    * \param _eps_abs
    * \param bc_mode
    * \param level
    */
    virtual int solve_ (MultiFab& _sol,
			Real      _eps_rel,
			Real      _eps_abs,
			MCBC_Mode bc_mode,
			int       level);

    /**
    * \brief Put the system in r-c form
    *
    * \param newrhs
    * \param oldrhs
    * \param initialsolution
    * \param inisol
    * \param bc_mode
    * \param level
    */
    void residualCorrectionForm (MultiFab&       newrhs,
				 const MultiFab& oldrhs,
				 MultiFab&       initialsolution,
				 const MultiFab& inisol,
				 MCBC_Mode       bc_mode,
				 int             level);

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
    virtual Real errorEstimate (int       level,
				MCBC_Mode bc_mode,
                                bool      local = false);

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
    */
    void relax (MultiFab& solL,
		MultiFab& rhsL,
		int       level,
		Real      eps_rel,
		Real      eps_abs,
		MCBC_Mode bc_mode);

    /**
    * \brief Perform relaxation at bottom of V-cycle
    *
    * \param solL
    * \param rhsL
    * \param level
    * \param eps_rel
    * \param eps_abs
    * \param bc_mode
    */
    virtual void coarsestSmooth (MultiFab& solL,
				 MultiFab& rhsL,
				 int       level,
				 Real      eps_rel,
				 Real      eps_abs,
				 MCBC_Mode bc_mode);
protected:

    //! Default flag, whether to use CG at bottom of MG cycle.
    static int def_usecg;

    /**
    * \brief Set flags, etc.
    */
    static void Initialize ();

    static void Finalize ();

    //! Default number of level, pre-, post and bottom smooths.
    static int def_nu_0, def_nu_1, def_nu_2, def_nu_f;

    //! Default number of post-cg relax passes.
    static int def_nu_b;

    //! Default maximum number of complete MG cycle iterations.
    static int def_maxiter;

    //! Default number of complete MG cycle iterations to perform.
    static int def_numiter;

    //! Default verbosity.
    static int def_verbose;

    //! Default relative, absolute tolerances for cg solve.
    static Real def_rtol_b, def_atol_b;

    //! Default maximum number of mg levels.
    static int def_numLevelsMAX;

    //! verbosity
    int verbose;

    //! Number of MG levels.  numlevels-1 is coarsest grid allowed.
    int numlevels;

    //! Number of components in solution.
    int numcomps;

    //! Current maximum number of allowed iterations.
    int maxiter;

    //! Current number of iterations to take.
    int numiter;

    //! Current number of level, pre-, post- and bottom smooths.
    int nu_0, nu_1, nu_2, nu_f;

    //! Current number of post-cg relax passes.
    int nu_b;

    //! Current flag, whether to use CG at bottom of MG cycle.
    int usecg;

    //! Relative, absolute tolerances.
    Real rtol_b, atol_b;

    //! Maximum number of mg levels.
    int numLevelsMAX;

    //! Internal temp data to store initial guess of solution.
    MultiFab* initialsolution;

    //! Internal temp data.
    Vector< MultiFab* > res;

    //! Internal temp data.
    Vector< MultiFab* > rhs;

    //! Internal temp data.
    Vector< MultiFab* > cor;

    //! Internal reference to linear operator.
    MCLinOp &Lp;

private:

    /**
    * \brief Disallow copy constructor, assignment operator
    *
    */
    MCMultiGrid (const MCMultiGrid&);
    MCMultiGrid& operator= (const MCMultiGrid&);
};

inline
MCLinOp&
MCMultiGrid::linOp ()
{
    return Lp;
}

inline
void
MCMultiGrid::setMaxIter (int _maxiter)
{
    maxiter = _maxiter;
}

inline
int
MCMultiGrid::getMaxIter () const
{
    return maxiter;
}

inline
void
MCMultiGrid::setNumIter (int _numiter)
{
    numiter = _numiter;
}

inline
int
MCMultiGrid::getNumIter () const
{
    return numiter;
}

inline
void
MCMultiGrid::setUseCG (int _usecg)
{
    usecg = _usecg;
}

inline
int
MCMultiGrid::getUseCG () const
{
    return usecg;
}

inline
int
MCMultiGrid::getNumLevels () const
{
    return numlevels;
}

inline
int
MCMultiGrid::getNumComponents () const
{
    return numcomps;
}

inline
void
MCMultiGrid::setVerbose (int _verbose)
{
    verbose = _verbose;
}

inline
int
MCMultiGrid::getVerbose ()
{
    return verbose;
}

inline
void
MCMultiGrid::set_preSmooth (int pre_smooth)
{
    nu_1 = pre_smooth;
}

inline
void
MCMultiGrid::set_postSmooth (int post_smooth)
{
    nu_2 = post_smooth;
}

inline
void
MCMultiGrid::set_cntRelax (int cnt_relax)
{
    nu_0 = cnt_relax;
}

inline
void
MCMultiGrid::set_finalSmooth (int final_smooth)
{
    nu_f = final_smooth;
}

inline
int
MCMultiGrid::preSmooth () const
{
    return nu_1;
}

inline
int
MCMultiGrid::postSmooth () const
{
    return nu_2;
}

inline
int
MCMultiGrid::cntRelax () const
{
    return nu_0;
}

inline
int
MCMultiGrid::finalSmooth () const
{
    return nu_f;
}

inline
void
MCMultiGrid::set_atol_b (Real atol)
{
    atol_b = atol;
}

inline
Real
MCMultiGrid::get_atol_b () const
{
    return atol_b;
}

inline
void
MCMultiGrid::set_rtol_b (Real rtol)
{
    rtol_b = rtol;
}

inline
Real
MCMultiGrid::get_rtol_b () const
{
    return rtol_b;
}

inline
void
MCMultiGrid::set_nu_b (int _nu_b)
{
    nu_b = _nu_b;
}

inline
int
MCMultiGrid::get_nu_b () const
{
    return nu_b;
}

}

#endif /*_MCMULTIGRID_H_*/
