//BL_COPYRIGHT_NOTICE

//
// $Id: CGSolver.cpp,v 1.15 2000-08-24 22:39:02 car Exp $
//

// Conjugate gradient support

#ifdef BL3_PROFILING
#include <BoxLib3/Profiler.H>
#endif

#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <CG_F.H>
#include <CGSolver.H>
#include <MultiGrid.H>

int CGSolver::initialized = 0;
int CGSolver::def_maxiter = 40;
int CGSolver::def_verbose = 0;
CGSolver::Solver CGSolver::def_cg_solver = CG;
double CGSolver::def_unstable_criterion = 10.;

namespace
{
void
spacer(ostream& os, int lev)
{
  for (int k = 0; k < lev; k++)
    {
      os << "   ";
    }
}
}

void
CGSolver::initialize ()
{
    ParmParse pp("cg");

    pp.query("maxiter", def_maxiter);
    pp.query("v", def_verbose);
    pp.query("unstable_criterion",def_unstable_criterion);
    int ii = 0;
    pp.query("cg_solver", ii);
    switch (ii)
      {
      case 0: def_cg_solver = CG; break;
      case 1: def_cg_solver = BiCGStab; break;
      case 2: def_cg_solver = CG_Alt; break;
      default:
	BoxLib::Error("CGSolver::initialize(): bad cg_solver");
      }

    if (ParallelDescriptor::IOProcessor() && def_verbose)
    {
        cout << "CGSolver settings...\n";
	cout << "   def_maxiter            = " << def_maxiter << '\n';
	cout << "   def_unstable_criterion = " << def_unstable_criterion << '\n';
	cout << "   def_cg_solver = " << def_cg_solver << '\n';
    }
    
    initialized = 1;
}

CGSolver::CGSolver (LinOp& _Lp,
            bool   _use_mg_precond,
            int    _lev)
    :
    Lp(_Lp),
    mg_precond(0),
    lev(_lev),
    use_mg_precond(_use_mg_precond),
    isExpert(false)
{
    if (!initialized)
        initialize();
    maxiter = def_maxiter;
    verbose = def_verbose;
    cg_solver = def_cg_solver;
    set_mg_precond();
}

void
CGSolver::set_mg_precond ()
{
    delete mg_precond;
    if (use_mg_precond)
    {
        mg_precond = new MultiGrid(Lp);
        mg_precond->setNumIter(1);
    }
}

CGSolver::~CGSolver ()
{
    delete mg_precond;
}

Real
CGSolver::norm (const MultiFab& res)
{
    //
    // Compute max-norm.
    //
    const int p       = 0;
    Real restot = 0.0;
    const BoxArray& gbox = res.boxArray();

    for (ConstMultiFabIterator mfi(res); mfi.isValid(); ++mfi) 
    {
        BL_ASSERT(mfi.validbox() == gbox[mfi.index()]);

        Real resk = mfi().norm(mfi.validbox(), p);

        if (p == 0)
        {
            restot = Max(restot, resk);
        }
        else if (p == 2)
        {
            restot += resk*resk;
        }
        else
        {
            BoxLib::Error("BOGUS P IN NORM" );
        } 
    }
    if (p == 0)
    {
        ParallelDescriptor::ReduceRealMax(restot);
    }
    else if (p == 2)
    {
        ParallelDescriptor::ReduceRealSum(restot);
        restot = sqrt(restot);
    }

    return restot;
}

int
CGSolver::solve (MultiFab&       sol,
                 const MultiFab& rhs,
                 Real            eps_rel,
                 Real            eps_abs,
                 LinOp::BC_Mode  bc_mode,
		 Solver solver)
{
  if ( solver == Automatic ) solver = cg_solver;
  if ( solver == CG )       return solve_cg(sol, rhs, eps_rel, eps_abs, bc_mode);
  if ( solver == BiCGStab ) return solve_bicgstab(sol, rhs, eps_rel, eps_abs, bc_mode);
  if ( solver == CG_Alt )   return solve_00(sol, rhs, eps_rel, eps_abs, bc_mode);
  return solve_00(sol, rhs, eps_rel, eps_abs, bc_mode);
}

int
CGSolver::solve_00 (MultiFab&       sol,
                 const MultiFab& rhs,
                 Real            eps_rel,
                 Real            eps_abs,
                 LinOp::BC_Mode  bc_mode)
{
#ifdef BL3_PROFILING
  BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::solve_00()");
#endif
    //
    // algorithm:
    //
    //   k=0;r=rhs-A*soln_0;
    //   while (||r_k||^2_2 > eps^2*||r_o||^2_2 && k < maxiter {
    //      k++
    //      solve Mz_k-1 = r_k-1 (if preconditioning, else z_k-1 = r_k-1)
    //      rho_k-1 = r_k-1^T z_k-1
    //      if (k=1) { p_1 = z_0 }
    //      else { beta = rho_k-1/rho_k-2; p = z + beta*p }
    //      w = Ap
    //      alpha = rho_k-1/p^tw
    //      x += alpha p
    //      r -= alpha w
    //   }
    //
    BL_ASSERT(sol.boxArray() == Lp.boxArray(lev));
    BL_ASSERT(rhs.boxArray() == Lp.boxArray(lev));

    int nghost = 1; int ncomp = sol.nComp();
    MultiFab* s = new MultiFab(sol.boxArray(), ncomp, nghost, Fab_allocate);
    MultiFab* r = new MultiFab(sol.boxArray(), ncomp, nghost, Fab_allocate);
    MultiFab* z = new MultiFab(sol.boxArray(), ncomp, nghost, Fab_allocate);
    MultiFab* w = new MultiFab(sol.boxArray(), ncomp, nghost, Fab_allocate);
    MultiFab* p = new MultiFab(sol.boxArray(), ncomp, nghost, Fab_allocate);
    //
    // Copy initial guess into a temp multifab guaranteed to have ghost cells.
    //
    int srccomp=0;  int destcomp=0;  ncomp=1;  nghost=0;
    s->copy(sol,srccomp,destcomp,ncomp);
    /*
      This routine assumes the LinOp is linear, and that when bc_mode =
      LinOp::Homogeneous_BC, LinOp::apply() on a zero vector will return a zero
      vector.  Given that, we define the problem we solve here from the
      original equation:

      Lp(sol) = rhs --> Lp(s) + Lp(sol,bc_mode=LinOp::Homogeneous_BC) = rhs

      where s is set to the incoming solution guess.  Rewriting,

      Lp(sol,bc_mode=LinOp::Homogeneous_BC) = r     [ = rhs - Lp(s) ].

      CG needs the residual of this equation on our initial guess.  But
      because we made the above assumption,

      r - Lp(sol,bc_mode=LinOp::Homogeneous_BC) = r = rhs - Lp(s)

      Which is simply the residual of the original equation evaluated at
      the initial guess.  Thus we get by with only one call to Lp.residual.
      Without this assumption, we'd need two.
    */
    Lp.residual((*r),  rhs, (*s), lev, bc_mode);
    //
    // Set initial guess for correction to 0.
    //
    sol.setVal(0.0);
    //
    // Set bc_mode=homogeneous
    //
    LinOp::BC_Mode temp_bc_mode=LinOp::Homogeneous_BC;

    Real rnorm  = norm(*r);
    Real rnorm0 = rnorm;
    Real minrnorm = rnorm;
    int ret = 0; // will return this value if all goes well

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
      spacer(cout, lev);
      cout << "CGsolver: Initial error (error0) =  " << rnorm0 << '\n';
    }

    Real beta, rho, rhoold = 0.0;
    /*
      The MultiFab copies used below to update z and p require nghost=0
      to avoid the possibility of filling valid regions with uninitialized
      data in the invalid regions of neighboring grids.  The default
      behavior in MultiFab copies will likely be changed in the future.
    */
    //
    // If eps_rel or eps_abs < 0: that test is effectively bypassed.
    //
    int nit = 1;
    for (;
         nit <= maxiter && rnorm > eps_rel*rnorm0 && rnorm > eps_abs;
         ++nit)
    {
        if (use_mg_precond)
        {
            //
            // Solve Mz_k-1 = r_k-1  and  rho_k-1 = r_k-1^T z_k-1
            //
            z->setVal(0.);
            mg_precond->solve(*z, *r, eps_rel, eps_abs, temp_bc_mode);
        }
        else
        {
            //
            // No preconditioner, z_k-1 = r_k-1  and  rho_k-1 = r_k-1^T r_k-1.
            //
            srccomp=0;  destcomp=0;  ncomp=1;
            z->copy((*r), srccomp, destcomp, ncomp);
        }

        rho = 0.0;
        int ncomp = z->nComp();
        const BoxArray& gbox = r->boxArray();

        for (MultiFabIterator rmfi(*r); rmfi.isValid(); ++rmfi)
        {
            DependentMultiFabIterator zmfi(rmfi, (*z));
            Real trho;
            BL_ASSERT(rmfi.validbox() == gbox[rmfi.index()]);
            FORT_CGXDOTY(&trho,zmfi().dataPtr(), 
                         ARLIM(zmfi().loVect()),ARLIM(zmfi().hiVect()),
                         rmfi().dataPtr(), 
                         ARLIM(rmfi().loVect()),ARLIM(rmfi().hiVect()),
                         rmfi.validbox().loVect(),rmfi.validbox().hiVect(),
                         &ncomp);
            rho += trho;
        }
        ParallelDescriptor::ReduceRealSum(rho);

        if (nit == 1)
        {
            //
            // k=1, p_1 = z_0
            //
            srccomp=0;  destcomp=0;  ncomp=1;  nghost=0;
            p->copy(*z, srccomp, destcomp, ncomp);
        }
        else
        {
            //
            // k>1, beta = rho_k-1/rho_k-2 and  p = z + beta*p
            //
            beta = rho/rhoold;
            advance(*p, beta, *z);
        }
        //
        //  w = Ap, and compute Transpose(p).w
        //
        Real pw = axp(*w, *p, temp_bc_mode);
        //
        // alpha = rho_k-1/p^tw
        //
	Real alpha;
	if( pw != 0. ){
	  alpha = rho/pw;
	}
	else {
	  ret = 1;
	  break;
	}
        
        if (ParallelDescriptor::IOProcessor() && verbose > 2)
        {
	  spacer(cout, lev);
            cout << "CGSolver:"
                 << " nit " << nit
                 << " pw "  << pw 
                 << " rho " << rho
                 << " alpha " << alpha;
            if (nit == 1)
            {
                cout << " beta undefined ...";
            }
            else
            {
                cout << " beta " << beta << " ...";
            }
        }
        //
        // x += alpha p  and  r -= alpha w
        //
        rhoold = rho;
        update(sol, alpha, *r, *p, *w);
        rnorm = norm(*r);
	if( rnorm > def_unstable_criterion*minrnorm ){
	  ret = 2;
	  break;
	}
	else if( rnorm < minrnorm ){
	  minrnorm = rnorm;
	}

        if (ParallelDescriptor::IOProcessor())
        {
            if (verbose > 1 ||
                (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
                  (eps_abs > 0. && rnorm < eps_abs)) && verbose))
            {
	      spacer(cout, lev);
                cout << "CGSolver: Iteration "
                     << nit
                     << " error/error0 "
                     << rnorm/rnorm0 << '\n';
            }
        }
    }
    
    if (ParallelDescriptor::IOProcessor())
      {
	if (verbose > 0 ||
	    (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
	      (eps_abs > 0. && rnorm < eps_abs)) && verbose))
	  {
	    spacer(cout, lev);
	    cout << "CGSolver: Final: Iteration "
		 << nit
		 << " error/error0 "
		 << rnorm/rnorm0 << '\n';
	  }
      }
    if( ret != 0 && isExpert == false ){
      BoxLib::Error("CGSolver:: apparent accuracy problem; try expert setting or change unstable_criterion");
    }
    if ( ret==0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
    {
        BoxLib::Error("CGSolver:: failed to converge!");
    }
    //
    // Omit ghost update since maybe not initialized in calling routine.
    // BoxLib_1.99 has no MultiFab::plus(MultiFab&) member, which would
    // operate only in valid regions; do explicitly.  Add to boundary
    // values stored in initialsolution.
    //
    if( ret == 0 ){
      srccomp=0; ncomp=1; nghost=0;
      sol.plus(*s,srccomp,ncomp,nghost);
    }

    delete s;
    delete r;
    delete w;
    delete p;
    delete z;
    return ret;
}

void
CGSolver::advance (MultiFab&       p,
                   Real            beta,
                   const MultiFab& z)
{
    //
    // Compute p = z  +  beta p
    //
    const BoxArray& gbox = Lp.boxArray(lev);
    int ncomp = p.nComp();

    for (MultiFabIterator pmfi(p); pmfi.isValid(); ++pmfi)
    {
        DependentMultiFabIterator zmfi(pmfi, z);

        BL_ASSERT(zmfi.validbox() == gbox[zmfi.index()]);

        FORT_CGADVCP(pmfi().dataPtr(),
                     ARLIM(pmfi().loVect()), ARLIM(pmfi().hiVect()),
                     zmfi().dataPtr(),
                     ARLIM(zmfi().loVect()), ARLIM(zmfi().hiVect()),
                     &beta, zmfi.validbox().loVect(), zmfi.validbox().hiVect(),
                     &ncomp);
    }
}

void
CGSolver::update (MultiFab&       sol,
                  Real            alpha,
                  MultiFab&       r,
                  const MultiFab& p,
                  const MultiFab& w)
{
    //
    // compute x =+ alpha p  and  r -= alpha w
    //
    const BoxArray& gbox = Lp.boxArray(lev);
    int ncomp = r.nComp();

    for (MultiFabIterator solmfi(sol); solmfi.isValid(); ++solmfi)
    {
        DependentMultiFabIterator rmfi(solmfi, r);
        DependentMultiFabIterator pmfi(solmfi, p);
        DependentMultiFabIterator wmfi(solmfi, w);

        BL_ASSERT(solmfi.validbox() == gbox[solmfi.index()]);

        FORT_CGUPDATE(solmfi().dataPtr(),
                      ARLIM(solmfi().loVect()), ARLIM(solmfi().hiVect()),
                      rmfi().dataPtr(),
                      ARLIM(rmfi().loVect()),   ARLIM(rmfi().hiVect()),
                      &alpha,
                      wmfi().dataPtr(),
                      ARLIM(wmfi().loVect()), ARLIM(wmfi().hiVect()),
                      pmfi().dataPtr(),
                      ARLIM(pmfi().loVect()), ARLIM(pmfi().hiVect()),
                      solmfi.validbox().loVect(), solmfi.validbox().hiVect(),
                      &ncomp);
    }
}

Real
CGSolver::axp (MultiFab&      w,
               MultiFab&      p,
               LinOp::BC_Mode bc_mode)
{
    //
    // Compute w = A.p, and return Transpose(p).w
    //
    Real pw = 0.0;
    const BoxArray& gbox = Lp.boxArray(lev);
    Lp.apply(w, p, lev, bc_mode);
    int ncomp = p.nComp();

    for (MultiFabIterator pmfi(p); pmfi.isValid(); ++pmfi)
    {
        DependentMultiFabIterator wmfi(pmfi, w);
        Real tpw;
        BL_ASSERT(pmfi.validbox() == gbox[pmfi.index()]);
        FORT_CGXDOTY(&tpw,
                     pmfi().dataPtr(),
                     ARLIM(pmfi().loVect()), ARLIM(pmfi().hiVect()),
                     wmfi().dataPtr(),
                     ARLIM(wmfi().loVect()), ARLIM(wmfi().hiVect()),
                     pmfi.validbox().loVect(), pmfi.validbox().hiVect(),
                     &ncomp);
        pw += tpw;
    }

    ParallelDescriptor::ReduceRealSum(pw);

    return pw;
}

namespace
{
void
sxay (MultiFab& ss, const MultiFab& xx, Real a, const MultiFab& yy)
{
#ifdef BL3_PROFILING
  BL3_PROFILE("CGSolver::sxay");
#endif
  const int ncomp = ss.nComp();

  for (MultiFabIterator smfi(ss); smfi.isValid(); ++smfi)
    {
      DependentMultiFabIterator xmfi(smfi, xx);
      DependentMultiFabIterator ymfi(smfi, yy);

      FORT_CGSXAY(smfi().dataPtr(),   ARLIM(smfi().loVect()), ARLIM(smfi().hiVect()),
		    xmfi().dataPtr(), ARLIM(xmfi().loVect()), ARLIM(xmfi().hiVect()),
		    &a,
		    ymfi().dataPtr(), ARLIM(ymfi().loVect()), ARLIM(ymfi().hiVect()),
		    smfi.validbox().loVect(), smfi.validbox().hiVect(),
		    &ncomp);
    }
}

Real
dotxy(const MultiFab& r, const MultiFab& z)
{
#ifdef BL3_PROFILING
  BL3_PROFILE("CGSolver::dotxy");
#endif
  int ncomp = z.nComp();
  Real rho = 0.0;
  for ( ConstMultiFabIterator rmfi(r); rmfi.isValid(); ++rmfi )
    {
      ConstDependentMultiFabIterator zmfi(rmfi, z);
      Real trho;
      FORT_CGXDOTY(&trho,
		   zmfi().dataPtr(), ARLIM(zmfi().loVect()),ARLIM(zmfi().hiVect()),
		   rmfi().dataPtr(), ARLIM(rmfi().loVect()),ARLIM(rmfi().hiVect()),
		   rmfi.validbox().loVect(),rmfi.validbox().hiVect(),
		   &ncomp);
      rho += trho;
    }
  ParallelDescriptor::ReduceRealSum(rho);
  return rho;
}
}

int
CGSolver::solve_bicgstab (MultiFab&       sol,
		       const MultiFab& rhs,
		       Real            eps_rel,
		       Real            eps_abs,
		       LinOp::BC_Mode  bc_mode)
{
#ifdef BL3_PROFILING
    BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::solve_bicgstab()");
#endif
  const int nghost = 1;
  const int ncomp  = 1;
  BL_ASSERT(sol.boxArray() == Lp.boxArray(lev));
  BL_ASSERT(rhs.boxArray() == Lp.boxArray(lev));
  BL_ASSERT(sol.nComp() == 1);

  MultiFab sorig(sol.boxArray(), ncomp, nghost);
  MultiFab s(sol.boxArray(), ncomp, nghost);
  MultiFab sh(sol.boxArray(), ncomp, nghost);
  MultiFab r(sol.boxArray(), ncomp, nghost);
  MultiFab rh(sol.boxArray(), ncomp, nghost);
  MultiFab p(sol.boxArray(), ncomp, nghost);
  MultiFab ph(sol.boxArray(), ncomp, nghost);
  MultiFab v(sol.boxArray(), ncomp, nghost);
  MultiFab t(sol.boxArray(), ncomp, nghost);

  //  cout << "eps_rel = " << eps_rel << endl;
  //  cout << "eps_abs = " << eps_abs << endl;
  //  cout << "lp.norm = " << Lp.norm(0, lev) << endl;
  //  cout << "sol.norm  = " << norm(sol) << endl;
  //  cout << "rhs.norm  = " << norm(rhs) << endl;

  sorig.copy(sol);
  Lp.residual(r, rhs, sorig, lev, bc_mode);
  rh.copy(r);
  sol.setVal(0.0);
  const LinOp::BC_Mode temp_bc_mode=LinOp::Homogeneous_BC;

  Real rnorm = norm(r);
  const Real rnorm0 = rnorm;

  if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
      spacer(cout, lev);
      cout << "CGSolver_bicgstab: Initial error (error0) =  " << rnorm0 << '\n';
    }
  int ret = 0;			// will return this value if all goes well
  Real rho_1;
  Real alpha, omega;
  int nit = 1;
  for (; nit <= maxiter; ++nit)
    {
      Real rho = dotxy(rh, r);
      if ( rho == 0 ) 
	{
	  BoxLib::Error("CGSolver_bicgstab fails");
	}
      if ( nit == 1 )
        {
	  p.copy(r);
        }
      else
        {
	  Real beta = (rho/rho_1)*(alpha/omega);
	  sxay(p, p, -omega, v);
	  sxay(p, r, beta, p);
        }

      if ( use_mg_precond )
        {
	  ph.setVal(0.0);
	  mg_precond->solve(ph, p, eps_rel, eps_abs, temp_bc_mode);
        }
      else
        {
	  ph.copy(p);
        }
      Lp.apply(v, ph, lev, temp_bc_mode);

      if ( Real rhTv = dotxy(rh, v) )
	{
	  alpha = rho/rhTv;
	}
      else
	{
	  ret = 1;
	  break;
	}

      sxay(s, r, -alpha, v);
      rnorm = norm(s);
      if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs )
	{
	  sxay(sol, sol, alpha, ph);
	  break;
	}
      if (ParallelDescriptor::IOProcessor())
        {
	  if (verbose > 1 ||
	      (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
		(eps_abs > 0. && rnorm < eps_abs)) && verbose))
            {
	      spacer(cout, lev);
	      cout << "CGSolver_bicgstab: Half Iter "
		   << nit
		   << " error/error0 "
		   << rnorm/rnorm0 << '\n';
            }
        }

      if ( use_mg_precond )
        {
	  sh.setVal(0.0);
	  mg_precond->solve(sh, s, eps_rel, eps_abs, temp_bc_mode);
        }
      else
        {
	  sh.copy(s);
        }
      Lp.apply(t, sh, lev, temp_bc_mode);
      if ( Real tTt = dotxy(t,t) )
	{
	  omega = dotxy(t,s)/tTt;
	}
      else
	{
	  ret = 2;
	  break;
	}
      sxay(sol, sol, alpha, ph);
      sxay(sol, sol, omega, sh);
      sxay(r, s, -omega, t);
      rnorm = norm(r);
      if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs )
	{
	  break;
	}
      if ( omega == 0 )
	{
	  ret = 3;
	  break;
	}
      if (ParallelDescriptor::IOProcessor())
        {
	  if (verbose > 1 ||
	      (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
		(eps_abs > 0. && rnorm < eps_abs)) && verbose))
            {
	      spacer(cout, lev);
	      cout << "CGSolver_bicgstab: Iteration "
		   << nit
		   << " error/error0 "
		   << rnorm/rnorm0 << '\n';
            }
        }
      rho_1 = rho;
    }
    
  //  cout << "norm(R) = " << rnorm << endl;
  //  cout << "AX+B = " << Lp.norm(0,lev)*norm(sol) + norm(rhs) << endl;
  
  if (ParallelDescriptor::IOProcessor())
    {
      if (verbose > 0 ||
	  (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
	    (eps_abs > 0. && rnorm < eps_abs)) && verbose))
	{
	  spacer(cout, lev);
	  cout << "CGSolver_bicgstab: Final: Iteration "
	       << nit
	       << " error/error0 "
	       << rnorm/rnorm0 << '\n';
	}
    }
  if( ret != 0 && isExpert == false )
    {
      BoxLib::Error("CGSolver_bicgstab:: apparent accuracy problem; try expert setting or change unstable_criterion");
    }
  if ( ret == 0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
    {
      BoxLib::Error("CGSolver_bicgstab:: failed to converge!");
    }

  if ( ret == 0 )
    {
      sol.plus(sorig, 0, 1, 0);
    }
  return ret;
}

int
CGSolver::solve_cg (MultiFab&       sol,
		    const MultiFab& rhs,
		    Real            eps_rel,
		    Real            eps_abs,
		    LinOp::BC_Mode  bc_mode)
{
#ifdef BL3_PROFILING
    BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::solve_cg()");
#endif
  const int nghost = 1;
  const int ncomp = sol.nComp();
  BL_ASSERT(sol.boxArray() == Lp.boxArray(lev));
  BL_ASSERT(rhs.boxArray() == Lp.boxArray(lev));
  BL_ASSERT(ncomp == 1 );

  MultiFab sorig(sol.boxArray(), ncomp, nghost);
  MultiFab r(sol.boxArray(), ncomp, nghost);
  MultiFab z(sol.boxArray(), ncomp, nghost);
  MultiFab q(sol.boxArray(), ncomp, nghost);
  MultiFab p(sol.boxArray(), ncomp, nghost);

  sorig.copy(sol);
  Lp.residual(r, rhs, sorig, lev, bc_mode);
  sol.setVal(0.0);
  const LinOp::BC_Mode temp_bc_mode=LinOp::Homogeneous_BC;

  Real rnorm  = norm(r);
  const Real rnorm0 = rnorm;
  Real minrnorm = rnorm;

  if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
      spacer(cout, lev);
      cout << "aCGsolver: Initial error (error0) =  " << rnorm0 << '\n';
    }

  int ret = 0;			// will return this value if all goes well
  Real rho_1;
  int nit = 1;
  for (; nit <= maxiter && rnorm > eps_rel*rnorm0 && rnorm > eps_abs; ++nit)
    {
      if (use_mg_precond)
        {
	  z.setVal(0.);
	  mg_precond->solve(z, r, eps_rel, eps_abs, temp_bc_mode);
        }
      else
        {
	  z.copy(r);
        }

      Real rho = dotxy(z,r);
      if (nit == 1)
        {
	  p.copy(z);
        }
      else
        {
	  Real beta = rho/rho_1;
	  sxay(p, z, beta, p);
        }

      Lp.apply(q, p, lev, temp_bc_mode);
      Real alpha;
      if ( Real pw = dotxy(p, q) )
	{
	  alpha = rho/pw;
	}
      else
	{
	  ret = 1;
	  break;
	}
        
      if (ParallelDescriptor::IOProcessor() && verbose > 2)
        {
	  spacer(cout, lev);
	  cout << "CGSolver_cg:"
	       << " nit " << nit
	       << " rho " << rho
	       << " alpha " << alpha;
        }
      sxay(sol, sol, alpha, p);
      sxay(  r,   r,-alpha, q);
      rnorm = norm(r);
      if( rnorm > def_unstable_criterion*minrnorm )
	{
	  ret = 2;
	  break;
	}
      else if( rnorm < minrnorm )
	{
	  minrnorm = rnorm;
	}

      if (ParallelDescriptor::IOProcessor())
        {
	  if (verbose > 1 ||
	      (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
		(eps_abs > 0. && rnorm < eps_abs)) && verbose))
            {
	      spacer(cout, lev);
	      cout << "CGSolver_cg: Iteration "
		   << nit
		   << " error/error0 "
		   << rnorm/rnorm0 << '\n';
            }
        }
      rho_1 = rho;
    }
    
  if (ParallelDescriptor::IOProcessor())
    {
      if (verbose > 0 ||
	  (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
	    (eps_abs > 0. && rnorm < eps_abs)) && verbose))
	{
	  spacer(cout, lev);
	  cout << "CGSolver_cg: Final: Iteration "
	       << nit
	       << " error/error0 "
	       << rnorm/rnorm0 << '\n';
	}
    }
  if ( ret != 0 && !isExpert )
    {
      BoxLib::Error("CGSolver_cg:: apparent accuracy problem; try expert setting or change unstable_criterion");
    }
  if ( ret == 0 &&  rnorm > eps_rel*rnorm0 && rnorm > eps_abs )
    {
      BoxLib::Error("CGSolver_cg:: failed to converge!");
    }
  if ( ret == 0 )
    {
      sol.plus(sorig, 0, 1, 0);
    }

  return ret;
}

