//BL_COPYRIGHT_NOTICE

#ifndef _MULTIGRID_H_
#define _MULTIGRID_H_

//
// $Id: MultiGrid.H,v 1.3 1998-07-07 17:25:15 lijewski Exp $
//

#include <Boolean.H>
#include <Array.H>
#include <PArray.H>
#include <Pointers.H>
#include <MultiFab.H>
#include <BndryData.H>
#include <LinOp.H>

//@Man:
/*@Memo:
  A MultiGrid solves the linear equation, L(phi)=rhs, for a LinOp L and
  MultiFabs rhs and phi using a V-type cycle of the MultiGrid algorithm
  */        
/*@Doc:
  A MultiGrid object solves the linear equation, L(phi)=rhs for a LinOp
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
  \begin{itemize}
  \item nu\_1(2) Number of passes of the pre-smoother
  \item nu\_2(2) Number of passes of the post-smoother
  \item nu\_0(1) Number of passes of the coarse-grid smoother per
  cycle
  \item nu\_f(8) Number of passes of the bottom smoother (if not using
  the conjugate-gradient bottom solve)
  \item maxiter(40) Maximum number of full multigrid cycles allowed to
  solve the system
  \item numiter(-1) Number of full multigrid cycles to perform
  (should be less than maxiter for fixed number of MG cycles;
  value is ignored if < 0)
  \item verbose(0) Verbosity (1-results, 2-progress)
  \item usecg(1) Whether to use the conjugate-gradient solver for the
  coarsest (bottom) solve of the multigrid hierarchy
  \item atol\_b(-1.0) Absolute error tolerance (<0 => ignored) for cg
  \item rtol\_b(.01) Relative error tolerance (<0 => ignored) for cg
  \item nu\_b(0) Number of passes of the bottom smoother taken
  {\it AFTER} the cg bottom solve (value ignored if <= 0)
  \item numLevelsMAX(1024) maximum number of mg levels
  \end{itemize}
        
  This class does NOT provide a copy constructor or assignment operator.
*/

class MultiGrid
{
public:
    //
    //@ManDoc: constructor
    //
    MultiGrid (LinOp& _Lp);
    //
    //@ManDoc: destructor
    //
    virtual ~MultiGrid ();
    //
    //@ManDoc: solve the system to relative err eps\_rel, absolute err eps\_abs
    //
    virtual void solve (MultiFab&       solution,
                        const MultiFab& _rhs,
                        Real            eps_rel = -1.0,
                        Real            eps_abs = -1.0,
                        LinOp::BC_Mode  bc_mode=LinOp::Inhomogeneous_BC);
    //
    //@ManDoc: return the linear operator
    //
    LinOp& linOp ();
    //
    //@ManDoc: set the maximum permitted multigrid iterations
    //
    void setMaxIter (int _maxiter);
    //
    //@ManDoc: return the maximum permitted multigrid iterations
    //
    int getMaxIter () const;
    //
    //@ManDoc: set the number of multigrid iterations to perform
    //
    void setNumIter (int _numiter);
    //
    //@ManDoc: return the number of multigrid iterations
    //
    int getNumIter () const;
    //
    //@ManDoc: set the flag for whether to use CGSolver at coarsest level
    //
    void setUseCG (int _usecg);
    //
    //@ManDoc: return the flag for whether to use CGSolver at coarsest level
    //
    int getUseCG () const;
    //
    //@ManDoc: set/return the number of multigrid levels
    //
    int getNumLevels (int _numlevels);
    //
    //@ManDoc: return the number of multigrid levels
    //
    int getNumLevels () const;
    //
    //@ManDoc: set the verbosity value
    //
    void setVerbose (int _verbose);
    //
    //@ManDoc: return the verbosity value
    //
    int getVerbose ();
    //
    //@ManDoc: set the number of passes of the pre-smoother
    //
    void set_preSmooth (int pre_smooth);
    //
    //@ManDoc: set the number of passes of the post-smoother
    //
    void set_postSmooth (int post_smooth);
    //
    //@ManDoc: set the number of passes of the coarse-grid smoother/mg iteration
    //
    void set_cntRelax (int cnt_relax);
    //
    //@ManDoc: set the number of passes of the bottom mg relaxer
    //
    void set_finalSmooth (int final_smooth);
    //
    //@ManDoc: Return the number of pre-smooth iterations at the level
    //
    int preSmooth () const;
    //
    //@ManDoc: Return the number of post-smooth iterations at the level
    //
    int postSmooth () const;
    //
    //@ManDoc: Return the number of level relaxations (not implemented)
    //
    int cntRelax () const;
    //
    //@ManDoc: Return the number of smoothing passes at bottom of MG (if no cg)
    //
    int finalSmooth () const;
    //
    //@ManDoc: set the maximum permitted absolute tolerance (<0 bypasses test)
    //
    void set_atol_b (Real atol);
    //
    //@ManDoc: get the maximum permitted absolute tolerance
    //
    Real get_atol_b () const;
    //
    //@ManDoc: set the maximum permitted relative tolerance (<0 bypasses test)
    //
    void set_rtol_b (Real rtol);
    //
    //@ManDoc: get the maximum permitted relative tolerance
    //
    Real get_rtol_b () const;
    //
    //@ManDoc: set the number of post-cg relax passes
    //
    void set_nu_b (int _nu_b);
    //
    //@ManDoc: set the number of post-cg relax passes
    //
    int get_nu_b () const;
    
protected:
    //@ManMemo: Internal members
    
    //
    //@ManDoc: Solve the linear system to relative and absolute tolerance
    //
    virtual int solve_ (MultiFab&      _sol,
                        Real           _eps_rel,
                        Real           _eps_abs,
                        LinOp::BC_Mode bc_mode,
                        int            level);
    //
    //@ManDoc: Put the system in r-c form
    //
    void residualCorrectionForm (MultiFab&       newrhs,
                                 const MultiFab& oldrhs,
                                 MultiFab&       initialsolution,
                                 const MultiFab& inisol,
                                 LinOp::BC_Mode  bc_mode,
                                 int             level);
    //
    //@ManDoc: Make space, set switches for new solution level
    //
    void prepareForLevel (int level);
    //
    //@ManDoc: Compute the number of multigrid levels, assuming ratio=2
    //
    int numLevels () const;
    //
    //@ManDoc: Return scalar estimate of error
    //
    virtual Real errorEstimate (int            level,
                                LinOp::BC_Mode bc_mode);
    //
    //@ManDoc: Transfer MultiFab from fine to coarse level
    //
    void average (MultiFab&       c,
                  const MultiFab& f);
    //
    //@ManDoc: Transfer MultiFab from coarse to fine level
    //
    void interpolate (MultiFab&       f,
                      const MultiFab& c);
    //
    //@ManDoc: Perform a MG V-cycle
    //
    void relax (MultiFab&      solL,
                MultiFab&      rhsL,
                int            level,
                Real           eps_rel,
                Real           eps_abs,
                LinOp::BC_Mode bc_mode);
    //
    //@ManDoc: Perform relaxation at bottom of V-cycle
    //
    virtual void coarsestSmooth (MultiFab&      solL,
                                 MultiFab&      rhsL,
                                 int            level,
                                 Real           eps_rel,
                                 Real           eps_abs,
                                 LinOp::BC_Mode bc_mode);
    //
    //@ManDoc: default flag, whether to use CG at bottom of MG cycle
    //
    static int def_usecg;
    //
    //@ManDoc: set flags, etc
    //
    static void initialize ();
    //
    //@ManDoc: default number of level, pre-, post- and bottom smooths
    //
    static int def_nu_0, def_nu_1, def_nu_2, def_nu_f;
    //
    //@ManDoc: default number of post-cg relax passes
    //
    static int def_nu_b;
    //
    //@ManDoc: default maximum number of complete MG cycle iterations
    //
    static int def_maxiter;
    //
    //@ManDoc: default number of complete MG cycle iterations to perform
    //
    static int def_numiter;
    //
    //@ManDoc: default verbosity
    //
    static int def_verbose;
    //
    //@ManDoc: default relative, absolute tolerance for cg solve
    //
    static Real def_rtol_b, def_atol_b;
    //
    //@ManDoc: default maximum number of mg levels
    //
    static int def_numLevelsMAX;
    //
    //@ManDoc: verbosity
    //
    int verbose;
    //
    //@ManDoc: Number of MG levels
    //
    int numlevels;
    //
    //@ManDoc: current maximum number of allowed iterations
    //
    int maxiter;
    //
    //@ManDoc: current number of iterations to take
    //
    int numiter;
    //
    //@ManDoc: current number of level, pre-, post- and bottom smooths
    //
    int nu_0, nu_1, nu_2, nu_f;
    //
    //@ManDoc: current number of post-cg relax passes
    //
    int nu_b;
    //
    //@ManDoc: current flag, whether to use CG at bottom of MG cycle
    //
    int usecg;
    //
    //@ManDoc: relative, absolute tolerance
    //
    Real rtol_b, atol_b;
    //
    //@ManDoc: maximum number of mg levels
    //
    int numLevelsMAX;
    //
    //@ManDoc: internal temp data to store initial guess of solution
    //
    MultiFab* initialsolution;
    //
    //@ManDoc: internal temp data
    //
    Array< MultiFab* > res;
    //
    //@ManDoc: internal temp data
    //
    Array< MultiFab* > rhs;
    //
    //@ManDoc: internal temp data
    //
    Array< MultiFab* > cor;
    //
    //@ManDoc: internal reference to linear operator
    //
    LinOp &Lp;
    
private:
    //
    // Flag, whether initialized.
    //
    static bool initialized;
    //
    // Disallow copy constructor, assignment operator
    //
    MultiGrid (const MultiGrid&);
    MultiGrid& operator= (const MultiGrid&);
};

inline
LinOp&
MultiGrid::linOp ()
{
    return Lp;
}

inline
void
MultiGrid::setMaxIter (int _maxiter)
{
    maxiter = _maxiter;
}

inline
int
MultiGrid::getMaxIter () const
{
    return maxiter;
}

inline
void
MultiGrid::setNumIter (int _numiter)
{
    numiter = _numiter;
}

inline
int
MultiGrid::getNumIter () const
{
    return numiter;
}

inline
void
MultiGrid::setUseCG (int _usecg)
{
    usecg = _usecg;
}

inline
int
MultiGrid::getUseCG () const
{
    return usecg;
}

inline
int
MultiGrid::getNumLevels () const
{
    return numlevels;
}

inline
void
MultiGrid::setVerbose (int _verbose)
{
    verbose = _verbose;
}

inline
int
MultiGrid::getVerbose ()
{
    return verbose;
}

inline
void
MultiGrid::set_preSmooth (int pre_smooth)
{
    nu_1 = pre_smooth;
}

inline
void
MultiGrid::set_postSmooth (int post_smooth)
{
    nu_2 = post_smooth;
}

inline
void
MultiGrid::set_cntRelax (int cnt_relax)
{
    nu_0 = cnt_relax;
}

inline
void
MultiGrid::set_finalSmooth (int final_smooth)
{
    nu_f = final_smooth;
}

inline
int
MultiGrid::preSmooth () const
{
    return nu_1;
}

inline
int
MultiGrid::postSmooth () const
{
    return nu_2;
}

inline
int
MultiGrid::cntRelax () const
{
    return nu_0;
}

inline
int
MultiGrid::finalSmooth () const
{
    return nu_f;
}

inline
void
MultiGrid::set_atol_b (Real atol)
{
    atol_b = atol;
}

inline
Real
MultiGrid::get_atol_b () const
{
    return atol_b;
}

inline
void
MultiGrid::set_rtol_b (Real rtol)
{
    rtol_b = rtol;
}

inline
Real
MultiGrid::get_rtol_b () const
{
    return rtol_b;
}

inline
void
MultiGrid::set_nu_b (int _nu_b)
{
    nu_b = _nu_b;
}

inline
int
MultiGrid::get_nu_b () const
{
    return nu_b;
}

inline
int
MultiGrid::getNumLevels (int _numlevels)
{
    assert(_numlevels >= 0);
    int oldnumlevels = numlevels;
    numlevels = Min(_numlevels, numLevels());
    return oldnumlevels;
}

#endif /*_MULTIGRID_H_*/
