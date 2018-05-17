
#include <algorithm>
#include <cstdlib>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MCCGSolver.H>
#include <AMReX_MG_F.H>
#include <AMReX_MCMultiGrid.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

namespace
{
    bool initialized = false;
}
//
// Set default values for these in Initialize()!!!
//
int  MCMultiGrid::def_nu_0;
int  MCMultiGrid::def_nu_1;
int  MCMultiGrid::def_nu_2;
int  MCMultiGrid::def_nu_f;
int  MCMultiGrid::def_maxiter;
int  MCMultiGrid::def_numiter;
int  MCMultiGrid::def_verbose;
int  MCMultiGrid::def_usecg;
Real MCMultiGrid::def_rtol_b;
Real MCMultiGrid::def_atol_b;
int  MCMultiGrid::def_nu_b;
int  MCMultiGrid::def_numLevelsMAX;

void
MCMultiGrid::Initialize ()
{
    if (initialized) return;
    //
    // Set defaults here!!!
    //
    MCMultiGrid::def_nu_0         = 1;
    MCMultiGrid::def_nu_1         = 2;
    MCMultiGrid::def_nu_2         = 2;
    MCMultiGrid::def_nu_f         = 8;
    MCMultiGrid::def_maxiter      = 40;
    MCMultiGrid::def_numiter      = -1;
    MCMultiGrid::def_verbose      = 0;
    MCMultiGrid::def_usecg        = 1;
    MCMultiGrid::def_rtol_b       = 0.01;
    MCMultiGrid::def_atol_b       = -1.0;
    MCMultiGrid::def_nu_b         = 0;
    MCMultiGrid::def_numLevelsMAX = 1024;

    ParmParse pp("mg");

    pp.query("maxiter",      def_maxiter);
    pp.query("numiter",      def_numiter);
    pp.query("nu_0",         def_nu_0);
    pp.query("nu_1",         def_nu_1);
    pp.query("nu_2",         def_nu_2);
    pp.query("nu_f",         def_nu_f);
    pp.query("v",            def_verbose);
    pp.query("usecg",        def_usecg);
    pp.query("rtol_b",       def_rtol_b);
    pp.query("bot_atol",     def_atol_b);
    pp.query("nu_b",         def_nu_b);
    pp.query("numLevelsMAX", def_numLevelsMAX);

    if (ParallelDescriptor::IOProcessor() && def_verbose)
    {
	std::cout << "def_nu_0         = " << def_nu_0         << '\n';
	std::cout << "def_nu_1         = " << def_nu_1         << '\n';
	std::cout << "def_nu_2         = " << def_nu_2         << '\n';
	std::cout << "def_nu_f         = " << def_nu_f         << '\n';
	std::cout << "def_maxiter      = " << def_maxiter      << '\n';
	std::cout << "def_usecg        = " << def_usecg        << '\n';
	std::cout << "def_rtol_b       = " << def_rtol_b       << '\n';
	std::cout << "def_atol_b       = " << def_atol_b       << '\n';
	std::cout << "def_nu_b         = " << def_nu_b         << '\n';
	std::cout << "def_numLevelsMAX = " << def_numLevelsMAX << '\n';
	std::cout << "def_verbose      = " << def_verbose      << '\n';
    }

    amrex::ExecOnFinalize(MCMultiGrid::Finalize);

    initialized = true;
}

void
MCMultiGrid::Finalize ()
{
    initialized = false;
}

static
Real
norm_inf (const MultiFab& res, bool local = false)
{
    Real restot = 0.0;
#ifdef _OPENMP
#pragma omp parallel reduction(max:restot)
#endif
    for (MFIter mfi(res,true); mfi.isValid(); ++mfi) 
    {
      restot = std::max(restot, res[mfi].norm(mfi.tilebox(), 0, 0, res.nComp()));
    }
    if ( !local )
        ParallelDescriptor::ReduceRealMax(restot);
    return restot;
}

static
void
Spacer (std::ostream& os, int lev)
{
    for (int k = 0; k < lev; k++)
    {
        os << "   ";
    }
}

MCMultiGrid::MCMultiGrid (MCLinOp &_lp)
    :
    initialsolution(0),
    Lp(_lp)
{
    Initialize();

    maxiter = def_maxiter;
    numiter = def_numiter;
    nu_0 = def_nu_0;
    nu_1 = def_nu_1;
    nu_2 = def_nu_2;
    nu_f = def_nu_f;
    usecg = def_usecg;
    verbose = def_verbose;
    rtol_b = def_rtol_b;
    atol_b = def_atol_b;
    nu_b = def_nu_b;
    numLevelsMAX = def_numLevelsMAX;
    numlevels = numLevels();
    numcomps = _lp.numberComponents();
    if ( ParallelDescriptor::IOProcessor() && (verbose > 2) )
    {
	BoxArray tmp = Lp.boxArray();
	std::cout << "MCMultiGrid: numlevels = " << numlevels 
		  << ": ngrid = " << tmp.size() << ", npts = [";
	for ( int i = 0; i < numlevels; ++i ) 
        {
	    if ( i > 0 ) tmp.coarsen(2);
	    std::cout << tmp.d_numPts() << " ";
        }
	std::cout << "]" << '\n';

	std::cout << "MCMultiGrid: " << numlevels
	     << " multigrid levels created for this solve" << '\n';
    }

    if ( ParallelDescriptor::IOProcessor() && (verbose > 4) )
    {
	std::cout << "Grids: " << '\n';
	BoxArray tmp = Lp.boxArray();
	for (int i = 0; i < numlevels; ++i)
	{
            Orientation face(0, Orientation::low);
            const DistributionMapping& map = Lp.bndryData().bndryValues(face).DistributionMap();
	    if (i > 0)
		tmp.coarsen(2);
	    std::cout << " Level: " << i << '\n';
	    for (int k = 0; k < tmp.size(); k++)
	    {
		const Box& b = tmp[k];
		std::cout << "  [" << k << "]: " << b << "   ";
		for (int j = 0; j < BL_SPACEDIM; j++)
		    std::cout << b.length(j) << ' ';
                std::cout << ":: " << map[k] << '\n';
	    }
	}
    }
}

MCMultiGrid::~MCMultiGrid ()
{
    delete initialsolution;

    for (int i = 0; i < cor.size(); ++i)
    {
	delete res[i];
	delete rhs[i];
	delete cor[i];
    }
}

Real
MCMultiGrid::errorEstimate (int       level,
			    MCBC_Mode bc_mode,
                            bool      local)
{
    Lp.residual(*res[level], *rhs[level], *cor[level], level, bc_mode);
    return norm_inf(*res[level], local);
}

void
MCMultiGrid::prepareForLevel (int level)
{
    //
    // Build this level by allocating reqd internal MultiFabs if necessary
    //
    if (cor.size() > level)
	return;
    res.resize(level+1, (MultiFab*)0);
    rhs.resize(level+1, (MultiFab*)0);
    cor.resize(level+1, (MultiFab*)0);
    Lp.prepareForLevel(level);
    if (cor[level] == 0)
    {
	const BoxArray& ba = Lp.boxArray(level);
	const DistributionMapping& dm = Lp.DistributionMap();
	res[level] = new MultiFab(ba,dm,numcomps,1);
	rhs[level] = new MultiFab(ba,dm,numcomps,1);
	cor[level] = new MultiFab(ba,dm,numcomps,1);
	if (level == 0)
        {
	    initialsolution = new MultiFab(ba, dm, numcomps, 1);
	    
	}
    }
}

void
MCMultiGrid::residualCorrectionForm (MultiFab&       resL,
				     const MultiFab& rhsL,
				     MultiFab&       solnL,
				     const MultiFab& inisol,
				     MCBC_Mode       bc_mode,
				     int             level)
{
    //
    // Using the linearity of the operator, Lp, we can solve this system
    // instead by solving for the correction required to the initial guess.
    //
    initialsolution->copy(inisol);
    solnL.copy(inisol);
    Lp.residual(resL, rhsL, solnL, level, bc_mode);
}

void
MCMultiGrid::solve (MultiFab&       _sol,
		    const MultiFab& _rhs,
		    Real            _eps_rel,
		    Real            _eps_abs,
		    MCBC_Mode       bc_mode)
{
    BL_ASSERT(numcomps == _sol.nComp());
    //
    // Prepare memory for new level, and solve the general boundary
    // value problem to within relative error _eps_rel.  Customized
    // to solve at level=0
    //
    int level = 0;
    prepareForLevel(level);
    residualCorrectionForm(*rhs[level],_rhs,*cor[level],_sol,bc_mode,level);
    if (!solve_(_sol, _eps_rel, _eps_abs, MCHomogeneous_BC, level))
      amrex::Error("MCMultiGrid::solve(): failed to converge!");
}

int
MCMultiGrid::solve_ (MultiFab& _sol,
		     Real      eps_rel,
		     Real      eps_abs,
		     MCBC_Mode bc_mode,
		     int       level)
{
  //
  // Relax system maxiter times, stop if relative error <= _eps_rel or
  // if absolute err <= _abs_eps
  //
  const Real strt_time = ParallelDescriptor::second();
  //
  // Elide a reduction by doing these together.
  //
  Real tmp[2] = { norm_inf(*rhs[level],true), errorEstimate(level,bc_mode,true) };

  ParallelDescriptor::ReduceRealMax(tmp,2);

  const Real norm_rhs  = tmp[0];
  const Real error0    = tmp[1];
  int        returnVal = 0;
  Real       error     = error0;

  if ( ParallelDescriptor::IOProcessor() && (verbose > 0) )
  {
      Spacer(std::cout, level);
      std::cout << "MCMultiGrid: Initial rhs                = " << norm_rhs << '\n';
      std::cout << "MCMultiGrid: Initial error (error0)     = " << error0 << '\n';
  }
  
  if ( ParallelDescriptor::IOProcessor() && eps_rel < 1.0e-16 && eps_rel > 0 )
  {
      std::cout << "MCMultiGrid: Tolerance "
                << eps_rel
                << " < 1e-16 is probably set too low" << '\n';
  }
  //
  // Initialize correction to zero at this level (auto-filled at levels below)
  //
  (*cor[level]).setVal(0.0);
  //
  // Note: if eps_rel, eps_abs < 0 then that test is effectively bypassed.
  //
  int        nit         = 1;
  const Real new_error_0 = norm_rhs;
  //const Real norm_Lp     = Lp.norm(0, level);


  for ( ;
        error > eps_abs &&
          error > eps_rel*norm_rhs &&
          nit <= maxiter;
        ++nit)
  {
    relax(*cor[level], *rhs[level], level, eps_rel, eps_abs, bc_mode);

    error = errorEstimate(level,bc_mode);
	
    if ( ParallelDescriptor::IOProcessor() && verbose > 1 )
    {
      const Real rel_error = (error0 != 0) ? error/new_error_0 : 0;
      Spacer(std::cout, level);
      std::cout << "MCMultiGrid: Iteration   "
                << nit
                << " error/error0 = "
                << rel_error << '\n';
    }
  }

  Real run_time = (ParallelDescriptor::second() - strt_time);
  if ( verbose > 0 )
  {
      if ( ParallelDescriptor::IOProcessor() )
      {
          const Real rel_error = (error0 != 0) ? error/error0 : 0;
          Spacer(std::cout, level);
          std::cout << "MCMultiGrid: Final Iter. "
                    << nit-1
                    << " error/error0 = "
                    << rel_error;
      }

      if ( verbose > 1 )
      {
        
        ParallelDescriptor::ReduceRealMax(run_time);

        if ( ParallelDescriptor::IOProcessor() )
          std::cout << ", Solve time: " << run_time << '\n';
      }
  }

  if ( ParallelDescriptor::IOProcessor() && (verbose > 0) )
  {
    if ( error < eps_rel*norm_rhs )
    {
      std::cout << "   Converged res < eps_rel*bnorm\n";
    } 
    else if ( error < eps_abs )
    {
      std::cout << "   Converged res < eps_abs\n";
    }
  }
  //
  // Omit ghost update since maybe not initialized in calling routine.
  // Add to boundary values stored in initialsolution.
  //
  _sol.copy(*cor[level]);
  _sol.plus(*initialsolution,0,_sol.nComp(),0);

  if ( error <= eps_rel*(norm_rhs) ||
       error <= eps_abs )
    returnVal = 1;

  //
  // Otherwise, failed to solve satisfactorily
  //
  return returnVal;
}

int
MCMultiGrid::getNumLevels (int _numlevels)
{
    BL_ASSERT(_numlevels >= 0);
    int oldnumlevels = numlevels;
    numlevels = std::min(_numlevels, numLevels());
    return oldnumlevels;
}

int
MCMultiGrid::numLevels () const
{
    int ng = Lp.numGrids();
    int lv = numLevelsMAX - 1;
    //
    // The routine `falls through' since coarsening and refining
    // a unit box does not yield the initial box.
    //
    const BoxArray &bs = Lp.boxArray(0);

    for (int i = 0; i < ng; ++i)
    {
	int llv = 0;
	Box tmp = bs[i];

	if (tmp.shortside() <= 3 )
	    amrex::Error("MCMultiGrid::numLevels(): fine grid too small");

	for (;;)
        {
	    Box ctmp  = tmp; ctmp.coarsen(2);
	    Box rctmp = ctmp; rctmp.refine(2);
	    if (tmp != rctmp)
		break;

	    if (ctmp.numPts() == 1)
		break;
            //
	    // wyc -- change so grid does not get too small
	    // wyc -- this is necessary for one-sided derivs in tensor-solve
	    // wyc -- when BBMGoUR is implemented, this will only be
	    // wyc -- necessary on finest level
            //
	    if (ctmp.shortside() <= 3)
		break;
	    llv++;
	    tmp = ctmp;
	}
        //
	// Set number of levels so that every box can be refined to there.
        //
	lv = (lv < llv)? lv : llv;
    }
    return lv+1;		// including coarsest
}

void
MCMultiGrid::relax (MultiFab& solL,
		    MultiFab& rhsL,
		    int       level,
		    Real      eps_rel,
		    Real      eps_abs,
		    MCBC_Mode bc_mode)
{
  //
  // Recursively relax system.  Equivalent to multigrid V-cycle.
  // At coarsest grid, call coarsestSmooth.
  //
  if (level < numlevels - 1 ) {
    for (int i = preSmooth() ; i > 0 ; i--) {
      Lp.smooth(solL, rhsL, level, bc_mode);
    }
    
    Lp.residual(*res[level], rhsL, solL, level, bc_mode);
    prepareForLevel(level+1);
    average(*rhs[level+1], *res[level]);
    cor[level+1]->setVal(0.0);
    for (int i = cntRelax(); i > 0 ; i--) {
      relax(*cor[level+1],*rhs[level+1],level+1,eps_rel,eps_abs,bc_mode);
    }
    interpolate(solL, *cor[level+1]);
    for (int i = postSmooth(); i > 0 ; i--) {
      Lp.smooth(solL, rhsL, level, bc_mode);
    }
  } else {
    coarsestSmooth(solL, rhsL, level, eps_rel, eps_abs, bc_mode);
  }
}

void
MCMultiGrid::coarsestSmooth (MultiFab& solL,
			     MultiFab& rhsL,
			     int       level,
			     Real      eps_rel,
			     Real      eps_abs,
			     MCBC_Mode bc_mode)
{
    prepareForLevel(level);
        
    if (usecg == 0)
    {
        for (int i = finalSmooth(); i > 0; i--)
            Lp.smooth(solL, rhsL, level, bc_mode);
    }
    else
    {
        bool use_mg_precond = false;
        MCCGSolver cg(Lp, use_mg_precond, level);
        cg.solve(solL, rhsL, rtol_b, atol_b, bc_mode);
        for(int i=0; i<nu_b; i++)
            Lp.smooth(solL, rhsL, level, bc_mode);
    }
}

void
MCMultiGrid::average (MultiFab&       c,
		      const MultiFab& f)
{
    //
    // Use Fortran function to average down (restrict) f to c.
    //
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter cmfi(c,true); cmfi.isValid(); ++cmfi)
    {
        const Box&       bx   = cmfi.tilebox();
	int              nc   = c.nComp();
        FArrayBox&       cfab = c[cmfi];
        const FArrayBox& ffab = f[cmfi];
	FORT_AVERAGE(
	    cfab.dataPtr(),ARLIM(cfab.loVect()),ARLIM(cfab.hiVect()),
	    ffab.dataPtr(),ARLIM(ffab.loVect()),ARLIM(ffab.hiVect()),
	    bx.loVect(), bx.hiVect(), &nc);
    }
}

void
MCMultiGrid::interpolate (MultiFab&       f,
			  const MultiFab& c)
{
    //
    // Use fortran function to interpolate up (prolong) c to f
    // Note: returns f=f+P(c) , i.e. ADDS interp'd c to f
    //
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter fmfi(c,true); fmfi.isValid(); ++fmfi)
    {
	const Box&       bx   = fmfi.tilebox();
	int              nc   = f.nComp();
        const FArrayBox& cfab = c[fmfi];
        FArrayBox&       ffab = f[fmfi];
	FORT_INTERP(
	    ffab.dataPtr(),ARLIM(ffab.loVect()),ARLIM(ffab.hiVect()),
	    cfab.dataPtr(),ARLIM(cfab.loVect()),ARLIM(cfab.hiVect()),
	    bx.loVect(), bx.hiVect(), &nc);
    }
}

}
