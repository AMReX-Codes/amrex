//
// $Id: CGSolver.cpp,v 1.40 2007-02-23 22:38:31 lijewski Exp $
//
#include <winstd.H>

#include <algorithm>
#include <iomanip>

#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <Profiler.H>
#include <Utility.H>
#include <LO_BCTYPES.H>
#include <CG_F.H>
#include <CGSolver.H>
#include <MultiGrid.H>

int              CGSolver::initialized            = 0;
int              CGSolver::def_maxiter            = 40;
int              CGSolver::def_verbose            = 1;
CGSolver::Solver CGSolver::def_cg_solver          = BiCGStab;
//CGSolver::Solver CGSolver::def_cg_solver          = CG;
double           CGSolver::def_unstable_criterion = 10.;
bool             CGSolver::use_jbb_precond        = 0;

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
CGSolver::initialize ()
{
    ParmParse pp("cg");

    pp.query("v", def_verbose);
    pp.query("maxiter", def_maxiter);
    pp.query("verbose", def_verbose);
    pp.query("use_jbb_precond", use_jbb_precond);
    pp.query("unstable_criterion",def_unstable_criterion);

    int ii;
    if (pp.query("cg_solver", ii))
    {
        switch (ii)
        {
        case 0: def_cg_solver = CG;       break;
        case 1: def_cg_solver = BiCGStab; break;
        default:
            BoxLib::Error("CGSolver::initialize(): bad cg_solver");
        }
    }

    if (ParallelDescriptor::IOProcessor() && def_verbose)
    {
        std::cout << "CGSolver settings ...\n";
	std::cout << "   def_maxiter            = " << def_maxiter            << '\n';
	std::cout << "   def_unstable_criterion = " << def_unstable_criterion << '\n';
	std::cout << "   def_cg_solver          = " << def_cg_solver          << '\n';
	std::cout << "   use_jbb_precond        = " << use_jbb_precond        << '\n';
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
    use_mg_precond(_use_mg_precond)
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
    }
}

CGSolver::~CGSolver ()
{
    delete mg_precond;
}

static
Real
norm_inf (const MultiFab& res, bool local = false)
{
    Real restot = 0.0;
    for (MFIter mfi(res); mfi.isValid(); ++mfi) 
    {
        restot = std::max(restot, res[mfi].norm(mfi.validbox(), 0));
    }
    if (!local)
        ParallelDescriptor::ReduceRealMax(restot);
    return restot;
}

int
CGSolver::solve (MultiFab&       sol,
                 const MultiFab& rhs,
                 Real            eps_rel,
                 Real            eps_abs,
                 LinOp::BC_Mode  bc_mode,
		 Solver          solver)
{
    int ret = -1;

    switch (solver)
    {
    case 0:
        ret = solve_cg(sol, rhs, eps_rel, eps_abs, bc_mode);
        break;
    case 1:
        ret = solve_bicgstab(sol, rhs, eps_rel, eps_abs, bc_mode);
        break;
    default:
        BoxLib::Error("CGSolver::solve(): unknown solver");
    }

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
    const BoxArray& zbox = z.boxArray();

    for (MFIter pmfi(p); pmfi.isValid(); ++pmfi)
    {
        BL_ASSERT(zbox[pmfi.index()] == gbox[pmfi.index()]);

        FORT_CGADVCP(p[pmfi].dataPtr(),
                     ARLIM(p[pmfi].loVect()), ARLIM(p[pmfi].hiVect()),
                     z[pmfi].dataPtr(),
                     ARLIM(z[pmfi].loVect()), ARLIM(z[pmfi].hiVect()),
                     &beta,
                     zbox[pmfi.index()].loVect(), zbox[pmfi.index()].hiVect(),
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

    for (MFIter solmfi(sol); solmfi.isValid(); ++solmfi)
    {
        BL_ASSERT(solmfi.validbox() == gbox[solmfi.index()]);

        FORT_CGUPDATE(sol[solmfi].dataPtr(),
                      ARLIM(sol[solmfi].loVect()), ARLIM(sol[solmfi].hiVect()),
                      r[solmfi].dataPtr(),
                      ARLIM(r[solmfi].loVect()),   ARLIM(r[solmfi].hiVect()),
                      &alpha,
                      w[solmfi].dataPtr(),
                      ARLIM(w[solmfi].loVect()), ARLIM(w[solmfi].hiVect()),
                      p[solmfi].dataPtr(),
                      ARLIM(p[solmfi].loVect()), ARLIM(p[solmfi].hiVect()),
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

    for (MFIter pmfi(p); pmfi.isValid(); ++pmfi)
    {
        Real tpw;
        BL_ASSERT(pmfi.validbox() == gbox[pmfi.index()]);
        FORT_CGXDOTY(&tpw,
                     p[pmfi].dataPtr(),
                     ARLIM(p[pmfi].loVect()), ARLIM(p[pmfi].hiVect()),
                     w[pmfi].dataPtr(),
                     ARLIM(w[pmfi].loVect()), ARLIM(w[pmfi].hiVect()),
                     pmfi.validbox().loVect(), pmfi.validbox().hiVect(),
                     &ncomp);
        pw += tpw;
    }

    ParallelDescriptor::ReduceRealSum(pw);

    return pw;
}

static
void
sxay (MultiFab& ss, const MultiFab& xx, Real a, const MultiFab& yy)
{
    const int ncomp = ss.nComp();

    for (MFIter smfi(ss); smfi.isValid(); ++smfi)
    {
        FORT_CGSXAY(ss[smfi].dataPtr(),
                    ARLIM(ss[smfi].loVect()), ARLIM(ss[smfi].hiVect()),
		    xx[smfi].dataPtr(),
                    ARLIM(xx[smfi].loVect()), ARLIM(xx[smfi].hiVect()),
		    &a,
		    yy[smfi].dataPtr(),
                    ARLIM(yy[smfi].loVect()), ARLIM(yy[smfi].hiVect()),
		    smfi.validbox().loVect(), smfi.validbox().hiVect(),
		    &ncomp);
    }
}

static
Real
dotxy (const MultiFab& r, const MultiFab& z, bool local = false)
{
    BL_PROFILE("CGSolver::dotxy");

    int ncomp = z.nComp();
    Real rho = 0.0;
    for (MFIter rmfi(r); rmfi.isValid(); ++rmfi)
    {
        Real trho;
        FORT_CGXDOTY(&trho,
                     z[rmfi].dataPtr(),
                     ARLIM(z[rmfi].loVect()),ARLIM(z[rmfi].hiVect()),
                     r[rmfi].dataPtr(),
                     ARLIM(r[rmfi].loVect()),ARLIM(r[rmfi].hiVect()),
                     rmfi.validbox().loVect(),rmfi.validbox().hiVect(),
                     &ncomp);
        rho += trho;
    }
    if (!local)
        ParallelDescriptor::ReduceRealSum(rho);
    return rho;
}

int
CGSolver::solve_bicgstab (MultiFab&       sol,
		       const MultiFab& rhs,
		       Real            eps_rel,
		       Real            eps_abs,
		       LinOp::BC_Mode  bc_mode)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::solve_bicgstab()");

    const int nghost = 1;
    const int ncomp  = 1;

    BL_ASSERT(sol.nComp() == 1);
    BL_ASSERT(sol.boxArray() == Lp.boxArray(lev));
    BL_ASSERT(rhs.boxArray() == Lp.boxArray(lev));

    MultiFab sorig(sol.boxArray(), ncomp, nghost);
    MultiFab s(sol.boxArray(), ncomp, nghost);
    MultiFab sh(sol.boxArray(), ncomp, nghost);
    MultiFab r(sol.boxArray(), ncomp, nghost);
    MultiFab rh(sol.boxArray(), ncomp, nghost);
    MultiFab p(sol.boxArray(), ncomp, nghost);
    MultiFab ph(sol.boxArray(), ncomp, nghost);
    MultiFab v(sol.boxArray(), ncomp, nghost);
    MultiFab t(sol.boxArray(), ncomp, nghost);

    if (verbose && false)
    {
        std::cout << "eps_rel = "       << eps_rel         << std::endl;
        std::cout << "eps_abs = "       << eps_abs         << std::endl;
        std::cout << "lp.norm = "       << Lp.norm(0, lev) << std::endl;
        std::cout << "sol.norm_inf = " << norm_inf(sol)   << std::endl;
        std::cout << "rhs.norm_inf = " << norm_inf(rhs)   << std::endl;
    }

    sorig.copy(sol);
    Lp.residual(r, rhs, sorig, lev, bc_mode);
    rh.copy(r);
    sol.setVal(0.0);
    const LinOp::BC_Mode temp_bc_mode=LinOp::Homogeneous_BC;

    Real       rnorm    = norm_inf(r);
    const Real rnorm0   = rnorm;
    const Real Lp_norm  = Lp.norm(0, lev);
    const Real rh_norm  = rnorm0;
    Real       sol_norm = 0.0;
  
    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        Spacer(std::cout, lev);
        std::cout << "CGSolver_bicgstab: Initial error (error0) =        " << rnorm0 << '\n';
    }
    int ret = 0, nit = 1;
    Real rho_1 = 0, alpha = 0, omega = 0;

    if ( rnorm == 0.0 || rnorm < eps_rel*(Lp_norm*sol_norm + rh_norm ) || rnorm < eps_abs )
    {
        if (verbose > 0 && ParallelDescriptor::IOProcessor())
	{
            Spacer(std::cout, lev);
            std::cout << "CGSolver_bicgstab: niter = 0,"
                      << ", rnorm = " << rnorm 
                      << ", eps_rel*(Lp_norm*sol_norm + rh_norm )" <<  eps_rel*(Lp_norm*sol_norm + rh_norm ) 
                      << ", eps_abs = " << eps_abs << std::endl;
	}
        return 0;
    }

    for (; nit <= maxiter; ++nit)
    {
        Real rho = dotxy(rh, r);
        if ( rho == 0 ) 
	{
            ret = 1;
            break;
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
            ret = 2;
            break;
	}
        sxay(sol, sol, alpha, ph);
        sxay(s, r, -alpha, v);
        rnorm = norm_inf(s);

        if (ParallelDescriptor::IOProcessor())
        {
            if (verbose > 1 ||
                (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
                  (eps_abs > 0. && rnorm < eps_abs)) && verbose))
            {
                Spacer(std::cout, lev);
                std::cout << "CGSolver_bicgstab: Half Iter "
                          << std::setw(11) << nit
                          << " rel. err. "
                          << rnorm/(rh_norm) << '\n';
            }
        }

#ifndef CG_USE_OLD_CONVERGENCE_CRITERIA
        sol_norm = norm_inf(sol);
        if ( rnorm < eps_rel*(Lp_norm*sol_norm + rh_norm ) || rnorm < eps_abs ) break;
#else
        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;
#endif

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
            ret = 3;
            break;
	}
        sxay(sol, sol, omega, sh);
        sxay(r, s, -omega, t);
        rnorm = norm_inf(r);

        if (ParallelDescriptor::IOProcessor())
        {
            if (verbose > 1 ||
                (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
                  (eps_abs > 0. && rnorm < eps_abs)) && verbose))
            {
                Spacer(std::cout, lev);
                std::cout << "CGSolver_bicgstab: Iteration "
                          << std::setw(11) << nit
                          << " rel. err. "
                          << rnorm/(rh_norm) << '\n';
            }
        }

#ifndef CG_USE_OLD_CONVERGENCE_CRITERIA
        sol_norm = norm_inf(sol);
        if ( rnorm < eps_rel*(Lp_norm*sol_norm + rh_norm ) || rnorm < eps_abs ) break;
#else
        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;
#endif
        if ( omega == 0 )
	{
            ret = 4;
            break;
	}
        rho_1 = rho;
    }

    if (ParallelDescriptor::IOProcessor())
    {
        if (verbose > 0 ||
            (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
              (eps_abs > 0. && rnorm < eps_abs)) && verbose))
	{
            Spacer(std::cout, lev);
            std::cout << "CGSolver_bicgstab: Final: Iteration "
                      << std::setw(4) << nit
                      << " rel. err. "
                      << rnorm/(rh_norm) << '\n';
	}
    }
#ifndef CG_USE_OLD_CONVERGENCE_CRITERIA
    if ( ret == 0 && rnorm > eps_rel*(Lp_norm*sol_norm + rh_norm ) && rnorm > eps_abs )
    {
        if ( ParallelDescriptor::IOProcessor() )
            BoxLib::Warning("CGSolver_bicgstab:: failed to converge!");
        ret = 8;
    }
#else
    if ( ret == 0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
    {
        if ( ParallelDescriptor::IOProcessor() )
            BoxLib::Warning("CGSolver_bicgstab:: failed to converge!");
        ret = 8;
    }
#endif

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rh_norm) )
    {
        sol.plus(sorig, 0, 1, 0);
    } 
    else 
    {
        sol.setVal(0.0);
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::solve_cg()");

    const int nghost = 1;
    const int ncomp  = sol.nComp();

    BL_ASSERT(ncomp == 1 );
    BL_ASSERT(sol.boxArray() == Lp.boxArray(lev));
    BL_ASSERT(rhs.boxArray() == Lp.boxArray(lev));

    MultiFab sorig(sol.boxArray(), ncomp, nghost);
    MultiFab r(sol.boxArray(), ncomp, nghost);
    MultiFab z(sol.boxArray(), ncomp, nghost);
    MultiFab q(sol.boxArray(), ncomp, nghost);
    MultiFab p(sol.boxArray(), ncomp, nghost);

    MultiFab r1(sol.boxArray(), ncomp, nghost);
    MultiFab z1(sol.boxArray(), ncomp, nghost);
    MultiFab r2(sol.boxArray(), ncomp, nghost);
    MultiFab z2(sol.boxArray(), ncomp, nghost);

    sorig.copy(sol);

    Lp.residual(r, rhs, sorig, lev, bc_mode);

    sol.setVal(0);

    const LinOp::BC_Mode temp_bc_mode=LinOp::Homogeneous_BC;

    Real       rnorm    = norm_inf(r);
    const Real rnorm0   = rnorm;
    Real       minrnorm = rnorm;

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        Spacer(std::cout, lev);
        std::cout << "              CG: Initial error :        " << rnorm0 << '\n';
    }

    const Real Lp_norm = Lp.norm(0, lev);
    const Real rh_norm = rnorm0;
    Real sol_norm      = 0;
    Real rho_1         = 0;
    int  ret           = 0;
    int  nit           = 1;

    if ( rnorm == 0.0 || rnorm < eps_rel*(Lp_norm*sol_norm + rh_norm ) || rnorm < eps_abs )
    {
        if (verbose > 0 && ParallelDescriptor::IOProcessor())
	{
            Spacer(std::cout, lev);
            std::cout << "       CG: niter = 0,"
                      << ", rnorm = " << rnorm 
                      << ", eps_rel*(Lp_norm*sol_norm + rh_norm )" <<  eps_rel*(Lp_norm*sol_norm + rh_norm ) 
                      << ", eps_abs = " << eps_abs << std::endl;
	}
        return 0;
    }

    for (; nit <= maxiter; ++nit)
    {
        if (use_jbb_precond && ParallelDescriptor::NProcs() > 1)
        {
            z.setVal(0);

            jbb_precond(z,r,lev,Lp);
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
            Spacer(std::cout, lev);
            std::cout << "CGSolver_cg:"
                      << " nit " << nit
                      << " rho " << rho
                      << " alpha " << alpha << '\n';
        }
        sxay(sol, sol, alpha, p);
        sxay(  r,   r,-alpha, q);
        rnorm = norm_inf(r);
        sol_norm = norm_inf(sol);

        if (ParallelDescriptor::IOProcessor())
        {
            if (verbose > 1 ||
                (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
                  (eps_abs > 0. && rnorm < eps_abs)) && verbose))
            {
                Spacer(std::cout, lev);
                std::cout << "       CG:       Iteration"
                          << std::setw(4) << nit
                          << " rel. err. "
                          << rnorm/(rh_norm) << '\n';
            }
        }

#ifndef CG_USE_OLD_CONVERGENCE_CRITERIA
        if ( rnorm < eps_rel*(Lp_norm*sol_norm + rh_norm) || rnorm < eps_abs ) break;
#else
        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;
#endif
      
        if ( rnorm > def_unstable_criterion*minrnorm )
	{
            ret = 2;
            break;
	}
        else if ( rnorm < minrnorm )
	{
            minrnorm = rnorm;
	}

        rho_1 = rho;
    }
    
    if (ParallelDescriptor::IOProcessor())
    {
        if (verbose > 0 ||
            (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
              (eps_abs > 0. && rnorm < eps_abs)) && verbose))
	{
            Spacer(std::cout, lev);
            std::cout << "       CG: Final Iteration"
                      << std::setw(4) << nit
                      << " rel. err. "
                      << rnorm/(rh_norm) << '\n';
	}
    }
#ifndef CG_USE_OLD_CONVERGENCE_CRITERIA
    if ( ret == 0 && rnorm > eps_rel*(Lp_norm*sol_norm + rh_norm) && rnorm > eps_abs )
    {
        if ( ParallelDescriptor::IOProcessor() )
            BoxLib::Warning("CGSolver_cg:: failed to converge!");
        ret = 8;
    }
#else
    if ( ret == 0 &&  rnorm > eps_rel*rnorm0 && rnorm > eps_abs )
    {
        if ( ParallelDescriptor::IOProcessor() )
            BoxLib::Warning("CGSolver_cg:: failed to converge!");
        ret = 8;
    }
#endif

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rh_norm) )
    {
        sol.plus(sorig, 0, 1, 0);
    } 
    else 
    {
        sol.setVal(0.0);
        sol.plus(sorig, 0, 1, 0);
    }

    return ret;
}

int
CGSolver::jbb_precond (MultiFab&       sol,
		       const MultiFab& rhs,
                       int             lev,
		       LinOp&          Lp)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::jbb_precond(MF)");
    //
    // This is a local routine.  No parallel is allowed to happen here.
    //
    int                  lev_loc = lev;
    const Real           eps_rel = 1.e-2;
    const Real           eps_abs = 1.e-16;
    const int            nghost  = 1;
    const int            ncomp   = sol.nComp();
    const bool           local   = true;
    const LinOp::BC_Mode bc_mode = LinOp::Homogeneous_BC;

    BL_ASSERT(ncomp == 1 );
    BL_ASSERT(sol.boxArray() == Lp.boxArray(lev_loc));
    BL_ASSERT(rhs.boxArray() == Lp.boxArray(lev_loc));

    MultiFab sorig(sol.boxArray(), ncomp, nghost);

    MultiFab r(sol.boxArray(), ncomp, nghost);
    MultiFab z(sol.boxArray(), ncomp, nghost);
    MultiFab q(sol.boxArray(), ncomp, nghost);
    MultiFab p(sol.boxArray(), ncomp, nghost);

    sorig.copy(sol);

    Lp.residual(r, rhs, sorig, lev_loc, LinOp::Homogeneous_BC, local);

    sol.setVal(0);

    Real       rnorm    = norm_inf(r,local);
    const Real rnorm0   = rnorm;
    Real       minrnorm = rnorm;

    if (verbose > 2 && ParallelDescriptor::IOProcessor())
    {
        Spacer(std::cout, lev_loc);
        std::cout << "     jbb_precond: Initial error :        " << rnorm0 << '\n';
    }

    const Real Lp_norm = Lp.norm(0, lev_loc, local);
    const Real rh_norm = rnorm0;
    Real sol_norm = 0;
    int  ret      = 0;			// will return this value if all goes well
    Real rho_1    = 0;
    int  nit      = 1;
    if ( rnorm == 0.0 || rnorm < eps_rel*(Lp_norm*sol_norm + rh_norm ) || rnorm < eps_abs )
    {
        if (verbose > 2 && ParallelDescriptor::IOProcessor())
	{
            Spacer(std::cout, lev_loc);
            std::cout << "jbb_precond: niter = 0,"
                      << ", rnorm = " << rnorm 
                      << ", eps_rel*(Lp_norm*sol_norm + rh_norm )" <<  eps_rel*(Lp_norm*sol_norm + rh_norm ) 
                      << ", eps_abs = " << eps_abs << std::endl;
	}
        return 0;
    }
    for (; nit <= maxiter; ++nit)
    {
        z.copy(r);

        Real rho = dotxy(z,r,local);
        if (nit == 1)
        {
            p.copy(z);
        }
        else
        {
            Real beta = rho/rho_1;
            sxay(p, z, beta, p);
        }

        Lp.apply(q, p, lev_loc, bc_mode, local);

        Real alpha;
        if ( Real pw = dotxy(p,q,local) )
	{
            alpha = rho/pw;
	}
        else
	{
            ret = 1;
            break;
	}
        
        if ( ParallelDescriptor::IOProcessor() && verbose > 3 )
        {
            Spacer(std::cout, lev_loc);
            std::cout << "jbb_precond:" << " nit " << nit
                      << " rho " << rho << " alpha " << alpha << '\n';
        }
        sxay(sol, sol, alpha, p);
        sxay(  r,   r,-alpha, q);
        rnorm    = norm_inf(r,   local);
        sol_norm = norm_inf(sol, local);

        if ( ParallelDescriptor::IOProcessor() )
        {
            if ( verbose > 3 ||
                 (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
                   (eps_abs > 0. && rnorm < eps_abs)) && verbose > 3) )
            {
                Spacer(std::cout, lev_loc);
                std::cout << "jbb_precond:       Iteration"
                          << std::setw(4) << nit
                          << " rel. err. "
                          << rnorm/(rh_norm) << '\n';
            }
        }

        if ( rnorm < eps_rel*(Lp_norm*sol_norm + rh_norm) || rnorm < eps_abs )
	{
            break;
	}
      
        if ( rnorm > def_unstable_criterion*minrnorm )
	{
            ret = 2;
            break;
	}
        else if ( rnorm < minrnorm )
	{
            minrnorm = rnorm;
	}

        rho_1 = rho;
    }
    
    if ( ParallelDescriptor::IOProcessor() )
    {
        if ( verbose > 2 ||
             (((eps_rel > 0. && rnorm < eps_rel*rnorm0) ||
               (eps_abs > 0. && rnorm < eps_abs)) && (verbose > 2)) )
	{
            Spacer(std::cout, lev_loc);
            std::cout << "jbb_precond: Final Iteration"
                      << std::setw(4) << nit
                      << " rel. err. "
                      << rnorm/(rh_norm) << '\n';
	}
    }
    if ( ret == 0 && rnorm > eps_rel*(Lp_norm*sol_norm + rh_norm) && rnorm > eps_abs )
    {
        if ( ParallelDescriptor::IOProcessor() )
	{
            BoxLib::Warning("jbb_precond:: failed to converge!");
	}
        ret = 8;
    }

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rh_norm) )
    {
        sol.plus(sorig, 0, 1, 0);
    } 
    else 
    {
        sol.setVal(0);
        sol.plus(sorig, 0, 1, 0);
    }

    return ret;
}
