//
// $Id: MultiGrid.cpp,v 1.34 2004-12-15 18:39:56 car Exp $
// 
#include <winstd.H>

#include <algorithm>
#include <cstdlib>

#include <ParmParse.H>
#include <Utility.H>
#include <ParallelDescriptor.H>
#include <CGSolver.H>
#include <MG_F.H>
#include <MultiGrid.H>
#include <Profiler.H>

bool MultiGrid::initialized     = false;
int MultiGrid::def_nu_0         = 1;
int MultiGrid::def_nu_1         = 2;
int MultiGrid::def_nu_2         = 2;
int MultiGrid::def_nu_f         = 8;
int MultiGrid::def_maxiter      = 40;
int MultiGrid::def_maxiter_b    = 80;
int MultiGrid::def_numiter      = -1;
int MultiGrid::def_verbose      = 0;
int MultiGrid::def_usecg        = 1;
#ifndef CG_USE_OLD_CONVERGENCE_CRITERIA
Real MultiGrid::def_rtol_b      = 0.0001;
#else
Real MultiGrid::def_rtol_b      = 0.01;
#endif
Real MultiGrid::def_atol_b      = -1.0;
int MultiGrid::def_nu_b         = 0;
int MultiGrid::def_numLevelsMAX = 1024;
int MultiGrid::def_smooth_on_cg_unstable = 0;
CGSolver::Solver MultiGrid::def_cg_solver    = CGSolver::BiCGStab;

static
Real
norm_inf (const MultiFab& res)
{
    Real restot = 0.0;
    for (MFIter mfi(res); mfi.isValid(); ++mfi) 
    {
        restot = std::max(restot, res[mfi].norm(mfi.validbox(), 0));
    }
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

void
MultiGrid::initialize ()
{
    ParmParse pp("mg");

    initialized = true;

    pp.query("maxiter", def_maxiter);
    pp.query("maxiter_b", def_maxiter_b);
    pp.query("numiter", def_numiter);
    pp.query("nu_0", def_nu_0);
    pp.query("nu_1", def_nu_1);
    pp.query("nu_2", def_nu_2);
    pp.query("nu_f", def_nu_f);
    pp.query("v", def_verbose);
    pp.query("verbose", def_verbose);
    pp.query("usecg", def_usecg);
    pp.query("rtol_b", def_rtol_b);
    pp.query("bot_atol", def_atol_b);
    pp.query("nu_b", def_nu_b);
    pp.query("numLevelsMAX", def_numLevelsMAX);
    pp.query("smooth_on_cg_unstable", def_smooth_on_cg_unstable);
    int ii;
    if (pp.query("cg_solver", ii ))
    {
        switch (ii)
        {
        case 0: def_cg_solver = CGSolver::CG; break;
        case 1: def_cg_solver = CGSolver::BiCGStab; break;
        case 2: def_cg_solver = CGSolver::CG_Alt; break;
        default:
            BoxLib::Error("MultiGrid::initialize(): bad cg_solver value");
        }
    }

    if (ParallelDescriptor::IOProcessor() && def_verbose)
    {
        std::cout << "MultiGrid settings...\n";
        std::cout << "   def_nu_0 =         " << def_nu_0         << '\n';
        std::cout << "   def_nu_1 =         " << def_nu_1         << '\n';
        std::cout << "   def_nu_2 =         " << def_nu_2         << '\n';
        std::cout << "   def_nu_f =         " << def_nu_f         << '\n';
        std::cout << "   def_maxiter =      " << def_maxiter      << '\n';
        std::cout << "   def_usecg =        " << def_usecg        << '\n';
        std::cout << "   def_maxiter_b =      " << def_maxiter_b  << '\n';
        std::cout << "   def_rtol_b =       " << def_rtol_b       << '\n';
        std::cout << "   def_atol_b =       " << def_atol_b       << '\n';
        std::cout << "   def_nu_b =         " << def_nu_b         << '\n';
        std::cout << "   def_numLevelsMAX = " << def_numLevelsMAX << '\n';
        std::cout << "   def_smooth_on_cg_unstable = " << def_smooth_on_cg_unstable << '\n';
        std::cout << "   def_cg_solver = " << def_cg_solver << '\n';
    }
}

MultiGrid::MultiGrid (LinOp &_Lp)
    :
    initialsolution(0),
    Lp(_Lp)
{
    if (!initialized)
        initialize();

    maxiter      = def_maxiter;
    numiter      = def_numiter;
    nu_0         = def_nu_0;
    nu_1         = def_nu_1;
    nu_2         = def_nu_2;
    nu_f         = def_nu_f;
    usecg        = def_usecg;
    verbose      = def_verbose;
    maxiter_b   = def_maxiter_b;
    rtol_b       = def_rtol_b;
    atol_b       = def_atol_b;
    nu_b         = def_nu_b;
    numLevelsMAX = def_numLevelsMAX;
    smooth_on_cg_unstable = def_smooth_on_cg_unstable;
    cg_solver    = def_cg_solver;
    numlevels    = numLevels();
    if ( ParallelDescriptor::IOProcessor() && verbose == 1 )
    {
	BoxArray tmp = Lp.boxArray();
	std::cout << "MultiGrid: numlevels = " << numlevels 
		  << ": ngrid = " << tmp.size() << ", npts = [";
	for ( int i = 0; i < numlevels; ++i ) 
        {
	    if ( i > 0 ) tmp.coarsen(2);
	    std::cout << tmp.d_numPts() << " ";
        }
	std::cout << "]" << std::endl;
    }
    if ( ParallelDescriptor::IOProcessor() && verbose > 2 )
    {
	std::cout << "MultiGrid: " << numlevels
	     << " multigrid levels created for this solve" << '\n';
	std::cout << "Grids: " << '\n';
	BoxArray tmp = Lp.boxArray();
	for (int i = 0; i < numlevels; ++i)
	{
	    if (i > 0)
		tmp.coarsen(2);
	    std::cout << " Level: " << i << '\n';
	    for (int k = 0; k < tmp.size(); k++)
	    {
		const Box& b = tmp[k];
		std::cout << "  [" << k << "]: " << b << "   ";
		for (int j = 0; j < BL_SPACEDIM; j++)
		    std::cout << b.length(j) << ' ';
		std::cout << '\n';
	    }
	}
    }
}

MultiGrid::~MultiGrid ()
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
MultiGrid::errorEstimate (int            level,
                          LinOp::BC_Mode bc_mode)
{
    Lp.residual(*res[level], *rhs[level], *cor[level], level, bc_mode);
    return norm_inf(*res[level]);
}

void
MultiGrid::prepareForLevel (int level)
{
    //
    // Build this level by allocating reqd internal MultiFabs if necessary.
    //
    if (cor.size() > level) return;

    res.resize(level+1, (MultiFab*)0);
    rhs.resize(level+1, (MultiFab*)0);
    cor.resize(level+1, (MultiFab*)0);

    Lp.prepareForLevel(level);

    if (cor[level] == 0)
    {
        res[level] = new MultiFab(Lp.boxArray(level), 1, 1, Fab_allocate);
        rhs[level] = new MultiFab(Lp.boxArray(level), 1, 1, Fab_allocate);
        cor[level] = new MultiFab(Lp.boxArray(level), 1, 1, Fab_allocate);
        if (level == 0)
        {
            initialsolution = new MultiFab(Lp.boxArray(0), 1, 1, Fab_allocate);
        }
    }
}

void
MultiGrid::residualCorrectionForm (MultiFab&       resL,
                                   const MultiFab& rhsL,
                                   MultiFab&       solnL,
                                   const MultiFab& inisol,
                                   LinOp::BC_Mode  bc_mode,
                                   int             level)
{
    //
    // Using the linearity of the operator, Lp, we can solve this system
    // instead by solving for the correction required to the initial guess.
    //

    initialsolution->copy(inisol);
    solnL.copy(inisol);
    Lp.residual(resL, rhsL, solnL, level, bc_mode);
    solnL.setVal(0.0);
}

void
MultiGrid::solve (MultiFab&       _sol,
                  const MultiFab& _rhs,
                  Real            _eps_rel,
                  Real            _eps_abs,
                  LinOp::BC_Mode  bc_mode)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::solve()");
    //
    // Prepare memory for new level, and solve the general boundary
    // value problem to within relative error _eps_rel.  Customized
    // to solve at level=0.
    //
    const int level = 0;
    prepareForLevel(level);
    residualCorrectionForm(*rhs[level],_rhs,*cor[level],_sol,bc_mode,level);
    if (!solve_(_sol, _eps_rel, _eps_abs, LinOp::Homogeneous_BC, level))
    {
        BoxLib::Error("MultiGrid:: failed to converge!");
    }
}

int
MultiGrid::solve_ (MultiFab&      _sol,
                   Real           eps_rel,
                   Real           eps_abs,
                   LinOp::BC_Mode bc_mode,
                   int            level)
{
  //
  // Relax system maxiter times, stop if relative error <= _eps_rel or
  // if absolute err <= _abs_eps
  //
  int        returnVal = 0;
  const Real error0    = errorEstimate(level, bc_mode);
  Real       error     = error0;
  if (ParallelDescriptor::IOProcessor() && verbose)
  {
      Spacer(std::cout, level);
      std::cout << "MultiGrid: Initial error (error0) = " << error0 << '\n';
  }

  if (ParallelDescriptor::IOProcessor() && eps_rel < 1.0e-16 && eps_rel > 0)
  {
      std::cout << "MultiGrid: Tolerance "
                << eps_rel
                << " < 1e-16 is probably set too low" << '\n';
  }
  //
  // Initialize correction to zero at this level (auto-filled at levels below)
  //
  (*cor[level]).setVal(0.0);
  //
  // Note: if eps_rel, eps_abs < 0 then that test is effectively bypassed
  //
  const Real norm_rhs    = norm_inf(*rhs[level]);
  const Real norm_Lp     = Lp.norm(0, level);
  const Real new_error_0 = norm_rhs;
  Real       norm_cor    = 0.0;
  int        nit         = 1;
  for ( ;
        error > eps_abs &&
            error > eps_rel*(norm_Lp*norm_cor+norm_rhs) &&
            nit <= maxiter;
        ++nit)
  {
      relax(*cor[level], *rhs[level], level, eps_rel, eps_abs, bc_mode);
      norm_cor = norm_inf(*cor[level]);
      error = errorEstimate(level, bc_mode);
	
      if (ParallelDescriptor::IOProcessor() && verbose > 1 )
      {
          const Real rel_error = (error0 != 0) ? error/new_error_0 : 0;
          Spacer(std::cout, level);
          std::cout << "MultiGrid: Iteration "
                    << nit
                    << " error/error0 "
                    << rel_error << '\n';
      }
  }

  if ( ParallelDescriptor::IOProcessor() && verbose )
  {
      const Real rel_error = (error0 != 0) ? error/error0 : 0;
      Spacer(std::cout, level);
      std::cout << "MultiGrid: iterations(" << nit 
		<< ") rel_error( " << rel_error << ")" << std::endl;
  }

  if ( nit == numiter                               ||
       error <= eps_rel*(norm_Lp*norm_cor+norm_rhs) ||
       error <= eps_abs )
  {
      //
      // Omit ghost update since maybe not initialized in calling routine.
      // Add to boundary values stored in initialsolution.
      //
      _sol.copy(*cor[level]);
      _sol.plus(*initialsolution,0,_sol.nComp(),0);
      returnVal = 1;
  }
  //
  // Otherwise, failed to solve satisfactorily
  //
  return returnVal;
}

int
MultiGrid::numLevels () const
{
    int ng = Lp.numGrids();
    int lv = numLevelsMAX;
    //
    // The routine `falls through' since coarsening and refining
    // a unit box does not yield the initial box.
    //
    const BoxArray& bs = Lp.boxArray(0);

    for (int i = 0; i < ng; ++i)
    {
        int llv = 0;
        Box tmp = bs[i];
        for (;;)
        {
            Box ctmp  = tmp;   ctmp.coarsen(2);
            Box rctmp = ctmp; rctmp.refine(2);
            if (tmp != rctmp || ctmp.numPts() == 1)
                break;
            llv++;
            tmp = ctmp;
        }
        //
        // Set number of levels so that every box can be refined to there.
        //
        if (lv >= llv)
            lv = llv;
    }

    return lv+1; // Including coarsest.
}

void
MultiGrid::relax (MultiFab&      solL,
                  MultiFab&      rhsL,
                  int            level,
                  Real           eps_rel,
                  Real           eps_abs,
                  LinOp::BC_Mode bc_mode)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::relax()");
    //
    // Recursively relax system.  Equivalent to multigrid V-cycle.
    // At coarsest grid, call coarsestSmooth.
    //
    if (level < numlevels - 1 )
    {
        for (int i = preSmooth() ; i > 0 ; i--)
        {
            Lp.smooth(solL, rhsL, level, bc_mode);
        }
        Lp.residual(*res[level], rhsL, solL, level, bc_mode);
        prepareForLevel(level+1);
        average(*rhs[level+1], *res[level]);
        cor[level+1]->setVal(0.0);
        for (int i = cntRelax(); i > 0 ; i--)
        {
            relax(*cor[level+1],*rhs[level+1],level+1,eps_rel,eps_abs,bc_mode);
        }
        interpolate(solL, *cor[level+1]);
        for (int i = postSmooth(); i > 0 ; i--)
        {
            Lp.smooth(solL, rhsL, level, bc_mode);
        }
    }
    else
    {
        coarsestSmooth(solL, rhsL, level, eps_rel, eps_abs, bc_mode, usecg);
    }
}

void
MultiGrid::coarsestSmooth (MultiFab&      solL,
                           MultiFab&      rhsL,
                           int            level,
                           Real           eps_rel,
                           Real           eps_abs,
                           LinOp::BC_Mode bc_mode,
                           int            local_usecg)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "coarsestSmooth()");

    prepareForLevel(level);
    if (local_usecg == 0)
    {
        Real error0;
        if (verbose)
        {
            error0 = errorEstimate(level, bc_mode);
            if (ParallelDescriptor::IOProcessor())
                std::cout << "   Bottom Smoother: Initial error (error0) = " 
                          << error0 << '\n';
        }

        for (int i = finalSmooth(); i > 0; i--)
        {
            Lp.smooth(solL, rhsL, level, bc_mode);

            if (verbose > 1 || (i == 1 && verbose))
            {
                Real error = errorEstimate(level, bc_mode);
                const Real rel_error = (error0 != 0) ? error/error0 : 0;
                if (ParallelDescriptor::IOProcessor())
                    std::cout << "   Bottom Smoother: Iteration "
                              << i
                              << " error/error0 "
                              << rel_error << '\n';
            }
        }
    }
    else
    {
        bool use_mg_precond = false;
	CGSolver cg(Lp, use_mg_precond, level);
	// cg.setExpert(true);
	cg.setMaxIter(maxiter_b);
	int ret = cg.solve(solL, rhsL, rtol_b, atol_b, bc_mode, cg_solver);
	if (ret != 0)
        {
            if (smooth_on_cg_unstable)
            {
                //
                // If the CG solver returns a nonzero value indicating 
                // the problem is unstable.  Assume this is not an accuracy 
                // issue and pound on it with the smoother.
	      // if ret == 8, then you have failure to converge
                //
                if (ParallelDescriptor::IOProcessor() && def_verbose)
                {
                    std::cout << "MultiGrid::coarsestSmooth(): CGSolver returns nonzero. Smoothing...." << std::endl;
                }

                coarsestSmooth(solL, rhsL, level, eps_rel, eps_abs, bc_mode, 0);
            }
            else
            {
                //
                // cg failure probably indicates loss of precision accident.
	      // if ret == 8, then you have failure to converge
                // setting solL to 0 should be ok.
                //
                solL.setVal(0);
                if (ParallelDescriptor::IOProcessor() && def_verbose)
                {
                    std::cout << "MultiGrid::coarsestSmooth(): setting coarse corr to zero" << std::endl;
                }
            }
	}
        for (int i = 0; i < nu_b; i++)
        {
            Lp.smooth(solL, rhsL, level, bc_mode);
        }
    }
}

void
MultiGrid::average (MultiFab&       c,
                    const MultiFab& f)
{
    //
    // Use Fortran function to average down (restrict) f to c.
    //
    for (MFIter cmfi(c); cmfi.isValid(); ++cmfi)
    {
        BL_ASSERT(c.boxArray().get(cmfi.index()) == cmfi.validbox());

        const Box& bx = cmfi.validbox();

        int nc = c.nComp();
        FORT_AVERAGE(c[cmfi].dataPtr(),
                     ARLIM(c[cmfi].loVect()), ARLIM(c[cmfi].hiVect()),
                     f[cmfi].dataPtr(),
                     ARLIM(f[cmfi].loVect()), ARLIM(f[cmfi].hiVect()),
                     bx.loVect(), bx.hiVect(), &nc);
    }
}

void
MultiGrid::interpolate (MultiFab&       f,
                        const MultiFab& c)
{
    //
    // Use fortran function to interpolate up (prolong) c to f
    // Note: returns f=f+P(c) , i.e. ADDS interp'd c to f.
    //
    for (MFIter fmfi(f); fmfi.isValid(); ++fmfi)
    {
        const Box& bx = c.boxArray()[fmfi.index()];
        int nc        = f.nComp();

        FORT_INTERP(f[fmfi].dataPtr(),
                    ARLIM(f[fmfi].loVect()), ARLIM(f[fmfi].hiVect()),
                    c[fmfi].dataPtr(),
                    ARLIM(c[fmfi].loVect()), ARLIM(c[fmfi].hiVect()),
                    bx.loVect(), bx.hiVect(), &nc);
    }
}
