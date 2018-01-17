
#include <algorithm>
#include <iomanip>
#include <cmath>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_VisMF.H>

#ifdef _OPENMP
#include <omp.h>
#endif

#define CG_USE_OLD_CONVERGENCE_CRITERIA 1

namespace amrex {

namespace {

static
void
sxay (MultiFab&       ss,
      const MultiFab& xx,
      Real            a,
      const MultiFab& yy,
      int             yycomp)
{
    BL_PROFILE("CGSolver::sxay()");

    const int ncomp  = 1;
    const int sscomp = 0;
    const int xxcomp = 0;
    MultiFab::LinComb(ss, 1.0, xx, xxcomp, a, yy, yycomp, sscomp, ncomp, 0);
}

inline
void
sxay (MultiFab&       ss,
      const MultiFab& xx,
      Real            a,
      const MultiFab& yy)
{
    sxay(ss,xx,a,yy,0);
}

}

MLCGSolver::MLCGSolver (MLLinOp& _lp)
    : Lp(_lp),
      amrlev(0),
      mglev(_lp.NMGLevels(0)-1)
{
}

MLCGSolver::~MLCGSolver ()
{
}

int
MLCGSolver::solve (MultiFab&       sol,
                   const MultiFab& rhs,
                   Real            eps_rel,
                   Real            eps_abs)
{
    BL_PROFILE_REGION("MLCGSolver::solve()");

    const int nghost = sol.nGrow(), ncomp = 1;

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();

    BL_ASSERT(sol.nComp() == ncomp);

    MultiFab ph(ba, dm, ncomp, nghost, MFInfo(), FArrayBoxFactory());
    MultiFab sh(ba, dm, ncomp, nghost, MFInfo(), FArrayBoxFactory());
    ph.setVal(0.0);
    sh.setVal(0.0);

    MultiFab sorig(ba, dm, ncomp, 0, MFInfo(), FArrayBoxFactory());
    MultiFab p    (ba, dm, ncomp, 0, MFInfo(), FArrayBoxFactory());
    MultiFab r    (ba, dm, ncomp, 0, MFInfo(), FArrayBoxFactory());
    MultiFab s    (ba, dm, ncomp, 0, MFInfo(), FArrayBoxFactory());
    MultiFab rh   (ba, dm, ncomp, 0, MFInfo(), FArrayBoxFactory());
    MultiFab v    (ba, dm, ncomp, 0, MFInfo(), FArrayBoxFactory());
    MultiFab t    (ba, dm, ncomp, 0, MFInfo(), FArrayBoxFactory());

    Lp.correctionResidual(amrlev, mglev, r, sol, rhs, MLLinOp::BCMode::Homogeneous);

    MultiFab::Copy(sorig,sol,0,0,1,0);
    MultiFab::Copy(rh,   r,  0,0,1,0);

    sol.setVal(0);

    Real rnorm = norm_inf(r);
    const Real rnorm0   = rnorm;

    if ( verbose > 0 && ParallelDescriptor::IOProcessor(p.color()) )
    {
        std::cout << "MLCGSolver_BiCGStab: Initial error (error0) =        " << rnorm0 << '\n';
    }
    int ret = 0, nit = 1;
    Real rho_1 = 0, alpha = 0, omega = 0;

    if ( rnorm0 == 0 || rnorm0 < eps_abs )
    {
        if ( verbose > 0 && ParallelDescriptor::IOProcessor(p.color()) )
	{
            std::cout << "MLCGSolver_BiCGStab: niter = 0,"
                      << ", rnorm = " << rnorm 
                      << ", eps_abs = " << eps_abs << std::endl;
	}
        return ret;
    }

    for (; nit <= maxiter; ++nit)
    {
        const Real rho = dotxy(rh,r);
        if ( rho == 0 ) 
	{
            ret = 1; break;
	}
        if ( nit == 1 )
        {
            MultiFab::Copy(p,r,0,0,1,0);
        }
        else
        {
            const Real beta = (rho/rho_1)*(alpha/omega);
            sxay(p, p, -omega, v);
            sxay(p, r,   beta, p);
        }
        MultiFab::Copy(ph,p,0,0,1,0);
        Lp.apply(amrlev, mglev, v, ph, MLLinOp::BCMode::Homogeneous);

        if ( Real rhTv = dotxy(rh,v) )
	{
            alpha = rho/rhTv;
	}
        else
	{
            ret = 2; break;
	}
        sxay(sol, sol,  alpha, ph);
        sxay(s,     r, -alpha,  v);

        rnorm = norm_inf(s);

        if ( verbose > 2 && ParallelDescriptor::IOProcessor(p.color()) )
        {
            std::cout << "MLCGSolver_BiCGStab: Half Iter "
                      << std::setw(11) << nit
                      << " rel. err. "
                      << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;

        MultiFab::Copy(sh,s,0,0,1,0);
        Lp.apply(amrlev, mglev, t, sh, MLLinOp::BCMode::Homogeneous);
        //
        // This is a little funky.  I want to elide one of the reductions
        // in the following two dotxy()s.  We do that by calculating the "local"
        // values and then reducing the two local values at the same time.
        //
        Real tvals[2] = { dotxy(t,t,true), dotxy(t,s,true) };

        ParallelAllReduce::Sum(tvals,2,Lp.BottomCommunicator());

        if ( tvals[0] )
	{
            omega = tvals[1]/tvals[0];
	}
        else
	{
            ret = 3; break;
	}
        sxay(sol, sol,  omega, sh);
        sxay(r,     s, -omega,  t);

        rnorm = norm_inf(r);

        if ( verbose > 2 && ParallelDescriptor::IOProcessor(p.color()) )
        {
            std::cout << "MLCGSolver_BiCGStab: Iteration "
                      << std::setw(11) << nit
                      << " rel. err. "
                      << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;

        if ( omega == 0 )
	{
            ret = 4; break;
	}
        rho_1 = rho;
    }

    if ( verbose > 0 && ParallelDescriptor::IOProcessor(p.color()) )
    {
        std::cout << "MLCGSolver_BiCGStab: Final: Iteration "
                  << std::setw(4) << nit
                  << " rel. err. "
                  << rnorm/(rnorm0) << '\n';
    }

    if ( ret == 0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
    {
        if ( ParallelDescriptor::IOProcessor(p.color()) )
            amrex::Warning("MLCGSolver_BiCGStab:: failed to converge!");
        ret = 8;
    }

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rnorm0) )
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

Real
MLCGSolver::dotxy (const MultiFab& r, const MultiFab& z, bool local)
{
    return Lp.xdoty(amrlev, mglev, r, z, local);
}

Real
MLCGSolver::norm_inf (const MultiFab& res, bool local)
{
    Real result = res.norm0(0,0,true);
    if (!local) {
        ParallelAllReduce::Max(result, Lp.BottomCommunicator());
    }
    return result;
}


}
