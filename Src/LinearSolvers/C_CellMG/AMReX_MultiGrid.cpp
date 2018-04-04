
#include <algorithm>
#include <cstdlib>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_CGSolver.H>
#include <AMReX_MG_F.H>
#include <AMReX_MultiGrid.H>

namespace amrex {

namespace
{
    bool initialized = false;
}
//
// Set default values for these in Initialize()!!!
//
int              MultiGrid::def_nu_0;
int              MultiGrid::def_nu_1;
int              MultiGrid::def_nu_2;
int              MultiGrid::def_nu_f;
int              MultiGrid::def_nu_b;
int              MultiGrid::def_usecg;
Real             MultiGrid::def_rtol_b;
Real             MultiGrid::def_atol_b;
int              MultiGrid::def_verbose;
int              MultiGrid::def_maxiter;
int              MultiGrid::def_maxiter_b;
int              MultiGrid::def_numLevelsMAX;
int              MultiGrid::def_smooth_on_cg_unstable;
int              MultiGrid::use_Anorm_for_convergence;

void
MultiGrid::Initialize ()
{
    if ( initialized ) return;
    //
    // Set defaults here!!!
    //
    MultiGrid::def_nu_0                  = 1;
    MultiGrid::def_nu_1                  = 2;
    MultiGrid::def_nu_2                  = 2;
    MultiGrid::def_nu_f                  = 8;
    MultiGrid::def_nu_b                  = 0;
    MultiGrid::def_usecg                 = 1;
#ifdef CG_USE_OLD_CONVERGENCE_CRITERIA
    MultiGrid::def_rtol_b                = 0.01;
#else
    MultiGrid::def_rtol_b                = 0.0001;
#endif
    MultiGrid::def_atol_b                = -1.0;
    MultiGrid::def_verbose               = 0;
    MultiGrid::def_maxiter               = 40;
    MultiGrid::def_maxiter_b             = 120;
    MultiGrid::def_numLevelsMAX          = 1024;
    MultiGrid::def_smooth_on_cg_unstable = 1;

    // This has traditionally been part of the stopping criteria, but for testing against
    //  other solvers it is convenient to be able to turn it off
    MultiGrid::use_Anorm_for_convergence = 1;

    ParmParse pp("mg");

    pp.query("v",                     def_verbose);
    pp.query("nu_0",                  def_nu_0);
    pp.query("nu_1",                  def_nu_1);
    pp.query("nu_2",                  def_nu_2);
    pp.query("nu_f",                  def_nu_f);
    pp.query("nu_b",                  def_nu_b);
    pp.query("usecg",                 def_usecg);
    pp.query("rtol_b",                def_rtol_b);
    pp.query("verbose",               def_verbose);
    pp.query("maxiter",               def_maxiter);
    pp.query("bot_atol",              def_atol_b);
    pp.query("maxiter_b",             def_maxiter_b);
    pp.query("numLevelsMAX",          def_numLevelsMAX);
    pp.query("smooth_on_cg_unstable", def_smooth_on_cg_unstable);

    pp.query("use_Anorm_for_convergence", use_Anorm_for_convergence);
#ifndef CG_USE_OLD_CONVERGENCE_CRITERIA
    if ( ParallelDescriptor::IOProcessor() && def_verbose > 2 )
    {
        if ( use_Anorm_for_convergence == 0 )
            std::cout << "It might be a good idea to define CG_USE_OLD_CONVERGENCE_CRITERIA\n";
    }
#endif

    if ( ParallelDescriptor::IOProcessor() && (def_verbose > 2) )
    {
        std::cout << "MultiGrid settings...\n";
        std::cout << "   def_nu_0                  = " << def_nu_0                  << '\n';
        std::cout << "   def_nu_1                  = " << def_nu_1                  << '\n';
        std::cout << "   def_nu_2                  = " << def_nu_2                  << '\n';
        std::cout << "   def_nu_f                  = " << def_nu_f                  << '\n';
        std::cout << "   def_nu_b                  = " << def_nu_b                  << '\n';
        std::cout << "   def_usecg                 = " << def_usecg                 << '\n';
        std::cout << "   def_rtol_b                = " << def_rtol_b                << '\n';
        std::cout << "   def_atol_b                = " << def_atol_b                << '\n';
        std::cout << "   def_maxiter               = " << def_maxiter               << '\n';
        std::cout << "   def_maxiter_b             = " << def_maxiter_b             << '\n';
        std::cout << "   def_numLevelsMAX          = " << def_numLevelsMAX          << '\n';
        std::cout << "   def_smooth_on_cg_unstable = " << def_smooth_on_cg_unstable << '\n';
        std::cout << "   use_Anorm_for_convergence = " << use_Anorm_for_convergence << '\n';
    }

    amrex::ExecOnFinalize(MultiGrid::Finalize);

    initialized = true;
}

void
MultiGrid::Finalize ()
{
    ;
}

static
Real
norm_inf (const MultiFab& res, bool local = false)
{
    return res.norm0(0, 0, local);
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

MultiGrid::MultiGrid (LinOp &_lp)
    :
    initialsolution(0),
    Lp(_lp)
{
    Initialize();

    maxiter      = def_maxiter;
    nu_0         = def_nu_0;
    nu_1         = def_nu_1;
    nu_2         = def_nu_2;
    nu_f         = def_nu_f;
    usecg        = def_usecg;
    verbose      = def_verbose;
    maxiter_b    = def_maxiter_b;
    rtol_b       = def_rtol_b;
    atol_b       = def_atol_b;
    nu_b         = def_nu_b;
    numLevelsMAX = def_numLevelsMAX;
    smooth_on_cg_unstable = def_smooth_on_cg_unstable;
    numlevels    = numLevels();

    do_fixed_number_of_iters = 0;

    if ( ParallelDescriptor::IOProcessor() && (verbose > 2) )
    {
	BoxArray tmp = Lp.boxArray();
	std::cout << "MultiGrid: numlevels = " << numlevels 
		  << ": ngrid = " << tmp.size() << ", npts = [";
	for ( int i = 0; i < numlevels; ++i ) 
        {
	    if ( i > 0 ) tmp.coarsen(2);
	    std::cout << tmp.d_numPts() << " ";
        }
	std::cout << "]" << '\n';

	std::cout << "MultiGrid: " << numlevels
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
                          LinOp::BC_Mode bc_mode,
                          bool           local)
{
    Lp.residual(*res[level], *rhs[level], *cor[level], level, bc_mode);
    return norm_inf(*res[level], local);
}

void
MultiGrid::prepareForLevel (int level)
{
    //
    // Build this level by allocating reqd internal MultiFabs if necessary.
    //
    if ( cor.size() > level ) return;

    res.resize(level+1, (MultiFab*)0);
    rhs.resize(level+1, (MultiFab*)0);
    cor.resize(level+1, (MultiFab*)0);

    Lp.prepareForLevel(level);

    if ( cor[level] == 0 )
    {
	const DistributionMapping& dm = Lp.DistributionMap();
	res[level] = new MultiFab(Lp.boxArray(level), dm, 1, Lp.NumGrow(), MFInfo(), FArrayBoxFactory());
	rhs[level] = new MultiFab(Lp.boxArray(level), dm, 1, Lp.NumGrow(), MFInfo(), FArrayBoxFactory());
	cor[level] = new MultiFab(Lp.boxArray(level), dm, 1, Lp.NumGrow(), MFInfo(), FArrayBoxFactory());
	if ( level == 0 )
	{
	    initialsolution = new MultiFab(Lp.boxArray(0), dm, 1, Lp.NumGrow(), MFInfo(), FArrayBoxFactory());
	}
    }
}

void
MultiGrid::solve (MultiFab&       _sol,
                  const MultiFab& _rhs,
                  Real            _eps_rel,
                  Real            _eps_abs,
                  LinOp::BC_Mode  bc_mode)
{
    //
    // Prepare memory for new level, and solve the general boundary
    // value problem to within relative error _eps_rel.  Customized
    // to solve at level=0.
    //
    const int level = 0;
    prepareForLevel(level);

    //
    // Copy the initial guess, which may contain inhomogeneous boundray conditions,
    // into both "initialsolution" (to be added back later) and into "cor[0]" which
    // we will only use here to compute the residual, then will set back to 0 below
    //
    initialsolution->copy(_sol);
    cor[level]->copy(_sol);

    //
    // Put the problem in residual-correction form: we will now use "rhs[level
    // the initial residual (rhs[0]) rather than the initial RHS (_rhs)
    // to begin the solve.
    //
    Lp.residual(*rhs[level],_rhs,*cor[level],level,bc_mode);

    //
    // Now initialize correction to zero at this level (auto-filled at levels below)
    //
    (*cor[level]).setVal(0.0); //

    //
    // Elide a reduction by doing these together.
    //
    Real tmp[2] = { norm_inf(_rhs,true), norm_inf(*rhs[level],true) };
    ParallelDescriptor::ReduceRealMax(tmp,2);
    if ( ParallelDescriptor::IOProcessor() && verbose > 0)
    {
        Spacer(std::cout, level);
        std::cout << "MultiGrid: Initial rhs                = " << tmp[0] << '\n';
        std::cout << "MultiGrid: Initial residual           = " << tmp[1] << '\n';
    }

    if (tmp[1] == 0.0)
	return;

    //
    // We can now use homogeneous bc's because we have put the problem into residual-correction form.
    //
    if ( !solve_(_sol, _eps_rel, _eps_abs, LinOp::Homogeneous_BC, tmp[0], tmp[1]) )
        amrex::Error("MultiGrid:: failed to converge!");
}

int
MultiGrid::solve_ (MultiFab&      _sol,
                   Real           eps_rel,
                   Real           eps_abs,
                   LinOp::BC_Mode bc_mode,
                   Real           bnorm,
                   Real           resnorm0)
{
    BL_PROFILE("MultiGrid::solve_()");

  //
  // If do_fixed_number_of_iters = 1, then do maxiter iterations without checking for convergence 
  // 
  // If do_fixed_number_of_iters = 0, then relax system maxiter times, 
  //    and stop if relative error <= _eps_rel or if absolute err <= _abs_eps
  //
  const Real strt_time = ParallelDescriptor::second();

  const int level = 0;

  //
  // We take the max of the norms of the initial RHS and the initial residual in order to capture both cases
  //
  Real norm_to_test_against;
  bool using_bnorm;
  if (bnorm >= resnorm0) 
  {
      norm_to_test_against = bnorm;
      using_bnorm          = true;
  } else {
      norm_to_test_against = resnorm0;
      using_bnorm          = false;
  } 

  int        returnVal = 0;
  Real       error     = resnorm0;

  //
  // Note: if eps_rel, eps_abs < 0 then that test is effectively bypassed
  //
  if ( ParallelDescriptor::IOProcessor() && eps_rel < 1.0e-16 && eps_rel > 0 )
  {
      std::cout << "MultiGrid: Tolerance "
                << eps_rel
                << " < 1e-16 is probably set too low" << '\n';
  }

  //
  // We initially define norm_cor based on the initial solution only so we can use it in the very first iteration
  //    to decide whether the problem is already solved (this is relevant if the previous solve used was only solved
  //    according to the Anorm test and not the bnorm test).
  //
  Real       norm_cor    = norm_inf(*initialsolution,true);
  ParallelDescriptor::ReduceRealMax(norm_cor);

  int        nit         = 1;
  const Real norm_Lp     = Lp.norm(0, level);
  Real       cg_time     = 0;

  if ( use_Anorm_for_convergence == 1 ) 
  {
     //
     // Don't need to go any further -- no iterations are required
     //
     if (error <= eps_abs || error < eps_rel*(norm_Lp*norm_cor+norm_to_test_against)) 
     {
         if ( ParallelDescriptor::IOProcessor() && (verbose > 0) )
         {
             std::cout << "   Problem is already converged -- no iterations required\n";
         }
         return 1;
     }

     for ( ;
           ( (error > eps_abs &&
              error > eps_rel*(norm_Lp*norm_cor+norm_to_test_against)) ||
             (do_fixed_number_of_iters == 1) )
             && nit <= maxiter;
           ++nit)
     {
         relax(*cor[level], *rhs[level], level, eps_rel, eps_abs, bc_mode, cg_time);

         Real tmp[2] = { norm_inf(*cor[level],true), errorEstimate(level,bc_mode,true) };

         ParallelDescriptor::ReduceRealMax(tmp,2);

         norm_cor = tmp[0];
         error    = tmp[1];

         if ( ParallelDescriptor::IOProcessor() && verbose > 1 )
         {
             const Real rel_error = error / norm_to_test_against;
             Spacer(std::cout, level);
             if (using_bnorm)
             {
                 std::cout << "MultiGrid: Iteration   "
                           << nit
                           << " resid/bnorm = "
                           << rel_error << '\n';
             } else {
                 std::cout << "MultiGrid: Iteration   "
                           << nit
                           << " resid/resid0 = "
                           << rel_error << '\n';
             }
         }
     }
  }
  else
  {
     //
     // Don't need to go any further -- no iterations are required
     //
     if (error <= eps_abs || error < eps_rel*norm_to_test_against) 
     {
         if ( ParallelDescriptor::IOProcessor() && (verbose > 0) )
         {
             std::cout << "   Problem is already converged -- no iterations required\n";
         }
         return 1;
     }

     for ( ;
           ( (error > eps_abs &&
              error > eps_rel*norm_to_test_against) ||
             (do_fixed_number_of_iters == 1) )
             && nit <= maxiter;
           ++nit)
     {
         relax(*cor[level], *rhs[level], level, eps_rel, eps_abs, bc_mode, cg_time);

         error = errorEstimate(level, bc_mode);
	
         if ( ParallelDescriptor::IOProcessor() && verbose > 1 )
         {
             const Real rel_error = error / norm_to_test_against;
             Spacer(std::cout, level);
             if (using_bnorm)
             {
                 std::cout << "MultiGrid: Iteration   "
                           << nit
                           << " resid/bnorm = "
                           << rel_error << '\n';
             } else {
                 std::cout << "MultiGrid: Iteration   "
                           << nit
                           << " resid/resid0 = "
                           << rel_error << '\n';
             }
         }
     }
  }

  Real run_time = (ParallelDescriptor::second() - strt_time);

  if ( verbose > 0 )
  {
      if ( ParallelDescriptor::IOProcessor() )
      {
          const Real rel_error = error / norm_to_test_against;
          Spacer(std::cout, level);
          if (using_bnorm)
          {
              std::cout << "MultiGrid: Iteration   "
                        << nit-1
                        << " resid/bnorm = "
                        << rel_error << '\n';
          } else {
              std::cout << "MultiGrid: Iteration   "
                        << nit-1
                        << " resid/resid0 = "
                        << rel_error << '\n';
             }
      }

      if ( verbose > 1 )
      {
          Real tmp[2] = { run_time, cg_time };

          ParallelDescriptor::ReduceRealMax(tmp,2);

          if ( ParallelDescriptor::IOProcessor() )
              std::cout << ", Solve time: " << tmp[0] << ", CG time: " << tmp[1];
      }

      if ( ParallelDescriptor::IOProcessor() ) std::cout << '\n';
  }

  if ( ParallelDescriptor::IOProcessor() && (verbose > 0) )
  {
      if ( do_fixed_number_of_iters == 1)
      {
          std::cout << "   Did fixed number of iterations: " << maxiter << std::endl;
      } 
      else if ( error < eps_rel*norm_to_test_against )
      {
          std::cout << "   Converged res < eps_rel*max(bnorm,res_norm)\n";
      } 
      else if ( (use_Anorm_for_convergence == 1) && (error < eps_rel*norm_Lp*norm_cor) )
      {
          std::cout << "   Converged res < eps_rel*Anorm*sol\n";
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

  if ( use_Anorm_for_convergence == 1 ) 
  {
     if ( do_fixed_number_of_iters == 1                ||
          error <= eps_rel*(norm_Lp*norm_cor+norm_to_test_against) ||
          error <= eps_abs )
       returnVal = 1;
  } 
  else 
  {
     if ( do_fixed_number_of_iters == 1 ||
          error <= eps_rel*(norm_to_test_against)   ||
          error <= eps_abs )
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
    int lv = numLevelsMAX-1;
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
            if ( tmp != rctmp || ctmp.numPts() == 1 )
                break;
            llv++;
            tmp = ctmp;
        }
        //
        // Set number of levels so that every box can be refined to there.
        //
        if ( lv >= llv )
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
                  LinOp::BC_Mode bc_mode,
                  Real&          cg_time)
{
    BL_PROFILE("MultiGrid::relax()");
    //
    // Recursively relax system.  Equivalent to multigrid V-cycle.
    // At coarsest grid, call coarsestSmooth.
    //
    if ( level < numlevels - 1 )
    {
        if ( verbose > 2 )
        {
           Real rnorm = errorEstimate(level, bc_mode);
           if (ParallelDescriptor::IOProcessor())
           {
              std::cout << "  AT LEVEL " << level << '\n';
              std::cout << "    DN:Norm before smooth " << rnorm << '\n';;
           }
        }
        for (int i = preSmooth() ; i > 0 ; i--)
        {
            Lp.smooth(solL, rhsL, level, bc_mode);
        }
        Lp.residual(*res[level], rhsL, solL, level, bc_mode);

        if ( verbose > 2 )
        {
           Real rnorm = norm_inf(*res[level]);
           if (ParallelDescriptor::IOProcessor())
              std::cout << "    DN:Norm after  smooth " << rnorm << '\n';
        }

        prepareForLevel(level+1);
        average(*rhs[level+1], *res[level]);
        cor[level+1]->setVal(0.0);
        for (int i = cntRelax(); i > 0 ; i--)
        {
            relax(*cor[level+1],*rhs[level+1],level+1,eps_rel,eps_abs,bc_mode,cg_time);
        }
        interpolate(solL, *cor[level+1]);

        if ( verbose > 2 )
        {
           Lp.residual(*res[level], rhsL, solL, level, bc_mode);
           Real rnorm = norm_inf(*res[level]);
           if ( ParallelDescriptor::IOProcessor() )
           {
              std::cout << "  AT LEVEL " << level << '\n';
              std::cout << "    UP:Norm before  smooth " << rnorm << '\n';
           }
        }

        for (int i = postSmooth(); i > 0 ; i--)
        {
            Lp.smooth(solL, rhsL, level, bc_mode);
        }
        if ( verbose > 2 )
        {
           Lp.residual(*res[level], rhsL, solL, level, bc_mode);
           Real rnorm = norm_inf(*res[level]);
           if ( ParallelDescriptor::IOProcessor() ) 
             std::cout << "    UP:Norm after  smooth " << rnorm << '\n';
        }
    }
    else
    {
        if ( verbose > 2 )
        {
           Real rnorm = norm_inf(rhsL);
           if ( ParallelDescriptor::IOProcessor() )
           {
              std::cout << "  AT LEVEL " << level << '\n';
              std::cout << "    DN:Norm before bottom " << rnorm << '\n';
           }
        }

        coarsestSmooth(solL, rhsL, level, eps_rel, eps_abs, bc_mode, usecg, cg_time);

        if ( verbose > 2 )
        {
           Lp.residual(*res[level], rhsL, solL, level, bc_mode);
           Real rnorm = norm_inf(*res[level]);
           if ( ParallelDescriptor::IOProcessor() ) 
              std::cout << "    UP:Norm after  bottom " << rnorm << '\n';
        }
    }
}

void
MultiGrid::coarsestSmooth (MultiFab&      solL,
                           MultiFab&      rhsL,
                           int            level,
                           Real           eps_rel,
                           Real           eps_abs,
                           LinOp::BC_Mode bc_mode,
                           int            local_usecg,
                           Real&          cg_time)
{
    BL_PROFILE("MultiGrid::coarsestSmooth()");
    prepareForLevel(level);

    if ( local_usecg == 0 )
    {
        Real error0 = 0;
        if ( verbose > 0 )
        {
            error0 = errorEstimate(level, bc_mode);
            if ( ParallelDescriptor::IOProcessor() )
                std::cout << "   Bottom Smoother: Initial error (error0) = " 
                          << error0 << '\n';
        }

        for (int i = finalSmooth(); i > 0; i--)
        {
            Lp.smooth(solL, rhsL, level, bc_mode);

            if ( verbose > 1 || (i == 1 && verbose) )
            {
                Real error = errorEstimate(level, bc_mode);
                const Real rel_error = (error0 != 0) ? error/error0 : 0;
                if ( ParallelDescriptor::IOProcessor() )
                    std::cout << "   Bottom Smoother: Iteration "
                              << i
                              << " error/error0 = "
                              << rel_error << '\n';
            }
        }
    }
    else
    {
        bool use_mg_precond = false;
	CGSolver cg(Lp, use_mg_precond, level);
	cg.setMaxIter(maxiter_b);

        const Real stime = ParallelDescriptor::second();

	int ret = cg.solve(solL, rhsL, rtol_b, atol_b, bc_mode);
        //
        // The whole purpose of cg_time is to accumulate time spent in CGSolver.
        //
        cg_time += (ParallelDescriptor::second() - stime);

	if ( ret != 0 )
        {
            if ( smooth_on_cg_unstable )
            {
                //
                // If the CG solver returns a nonzero value indicating 
                // the problem is unstable.  Assume this is not an accuracy 
                // issue and pound on it with the smoother.
                // if ret == 8, then you have failure to converge
                //
                if ( ParallelDescriptor::IOProcessor() && (verbose > 0) )
                    std::cout << "MultiGrid::coarsestSmooth(): CGSolver returns nonzero. Smoothing ...\n";

                coarsestSmooth(solL, rhsL, level, eps_rel, eps_abs, bc_mode, 0, cg_time);
            }
            else
            {
                //
                // cg failure probably indicates loss of precision accident.
                // if ret == 8, then you have failure to converge
                // setting solL to 0 should be ok.
                //
                solL.setVal(0);
                if ( ParallelDescriptor::IOProcessor() && (verbose > 0) )
                {
                    std::cout << "MultiGrid::coarsestSmooth(): setting coarse corr to zero" << '\n';
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
    BL_PROFILE("MultiGrid::average()");
    //
    // Use Fortran function to average down (restrict) f to c.
    //
    const bool tiling = true;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter cmfi(c,tiling); cmfi.isValid(); ++cmfi)
    {
        BL_ASSERT(c.boxArray().get(cmfi.index()) == cmfi.validbox());

        const int        nc   = c.nComp();
        const Box&       bx   = cmfi.tilebox();
        FArrayBox&       cfab = c[cmfi];
        const FArrayBox& ffab = f[cmfi];

        FORT_AVERAGE(cfab.dataPtr(),
                     ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                     ffab.dataPtr(),
                     ARLIM(ffab.loVect()), ARLIM(ffab.hiVect()),
                     bx.loVect(), bx.hiVect(), &nc);
    }
}

void
MultiGrid::interpolate (MultiFab&       f,
                        const MultiFab& c)
{
    BL_PROFILE("MultiGrid::interpolate()");
    //
    // Use fortran function to interpolate up (prolong) c to f
    // Note: returns f=f+P(c) , i.e. ADDS interp'd c to f.
    //
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(c,true); mfi.isValid(); ++mfi)
    {
        const Box&         bx = mfi.tilebox();
        const int          nc = f.nComp();
        const FArrayBox& cfab = c[mfi];
        FArrayBox&       ffab = f[mfi];

        FORT_INTERP(ffab.dataPtr(),
                    ARLIM(ffab.loVect()), ARLIM(ffab.hiVect()),
                    cfab.dataPtr(),
                    ARLIM(cfab.loVect()), ARLIM(cfab.hiVect()),
                    bx.loVect(), bx.hiVect(), &nc);
    }
}

int
MultiGrid::getNumLevels (int _numlevels)
{
    BL_ASSERT(_numlevels >= 0);
    int oldnumlevels = numlevels;
    numlevels = std::min(_numlevels, numLevels());
    return oldnumlevels;
}

}
