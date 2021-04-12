
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_MLMG.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <limits>
#include <algorithm>
#include <iomanip>
#include <cmath>


namespace amrex {

namespace {

template <class MF>
static
void
sxay (FabArray<MF>&       dst,
      const FabArray<MF>& x,
      Real            b,
      const FabArray<MF>& y,
      int             ycomp,
      int             nghost)
{
    BL_PROFILE("CGSolver::sxay()");

    Real a = 1.0;
    const int numcomp  = dst.nComp();
    const int dstcomp = 0;
    const int xcomp = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);

        if (bx.ok()) {
            auto const xfab =   x.array(mfi);
            auto const yfab =   y.array(mfi);
            auto       dfab = dst.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FUSIBLE ( bx, numcomp, i, j, k, n,
            {
                dfab(i,j,k,dstcomp+n) = a*xfab(i,j,k,xcomp+n) + b*yfab(i,j,k,ycomp+n);
            });
        }
    }

}

template <class MF>
inline
void
sxay (FabArray<MF>&       ss,
      const FabArray<MF>& xx,
      Real            a,
      const FabArray<MF>& yy,
      const int       nghost)
{
    sxay(ss,xx,a,yy,0,nghost);
}

}

template<class MF>
MLCGSolver<MF>::MLCGSolver (MLMG* a_mlmg, MLLinOp& _lp, Type _typ)
    : mlmg(a_mlmg),
      Lp(_lp),
      solver_type(_typ),
      amrlev(0),
      mglev(_lp.NMGLevels(0)-1)
{
    amrex::ignore_unused(mlmg);
}

template<class MF>
MLCGSolver<MF>::~MLCGSolver ()
{
}

template<class MF>
int
MLCGSolver<MF>::solve (MultiFab&       sol,
                   const MultiFab& rhs,
                   Real            eps_rel,
                   Real            eps_abs)
{
    if (solver_type == Type::BiCGStab) {
        return solve_bicgstab(sol,rhs,eps_rel,eps_abs);
    } else {
        return solve_cg(sol,rhs,eps_rel,eps_abs);
    }
}

template<class MF>
int
MLCGSolver<MF>::solve_bicgstab (FabArray<MF>&       sol,
                                const FabArray<MF>& rhs,
                            Real            eps_rel,
                            Real            eps_abs)
{
    BL_PROFILE("MLCGSolver::bicgstab");

    const int ncomp = sol.nComp();

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();
    const auto& factory = sol.Factory();

    FabArray<MF> ph(ba, dm, ncomp, sol.nGrowVect(), MFInfo(), factory);
    FabArray<MF> sh(ba, dm, ncomp, sol.nGrowVect(), MFInfo(), factory);
    ph.setVal(0.0);
    sh.setVal(0.0);

    FabArray<MF> sorig(ba, dm, ncomp, nghost, MFInfo(), factory);
    FabArray<MF> p    (ba, dm, ncomp, nghost, MFInfo(), factory);
    FabArray<MF> r    (ba, dm, ncomp, nghost, MFInfo(), factory);
    FabArray<MF> s    (ba, dm, ncomp, nghost, MFInfo(), factory);
    FabArray<MF> rh   (ba, dm, ncomp, nghost, MFInfo(), factory);
    FabArray<MF> v    (ba, dm, ncomp, nghost, MFInfo(), factory);
    FabArray<MF> t    (ba, dm, ncomp, nghost, MFInfo(), factory);

    Lp.correctionResidual(amrlev, mglev, r, sol, rhs, MLLinOp::BCMode::Homogeneous);

    // If singular remove mean from residual
//    if (Lp.isBottomSingular()) mlmg->makeSolvable(amrlev, mglev, r);

    // Then normalize
    Lp.normalize(amrlev, mglev, r);

    amrex::Copy(sorig,sol,0,0,ncomp,IntVect(nghost));
    amrex::Copy(rh,   r,  0,0,ncomp,IntVect(nghost));

    sol.setVal(0);

    Real rnorm = norm_inf(r);
    const Real rnorm0   = rnorm;

    if ( verbose > 0 )
    {
        amrex::Print() << "MLCGSolver_BiCGStab: Initial error (error0) =        " << rnorm0 << '\n';
    }
    int ret = 0;
    iter = 1;
    Real rho_1 = 0, alpha = 0, omega = 0;

    if ( rnorm0 == 0 || rnorm0 < eps_abs )
    {
        if ( verbose > 0 )
        {
            amrex::Print() << "MLCGSolver_BiCGStab: niter = 0,"
                           << ", rnorm = " << rnorm
                           << ", eps_abs = " << eps_abs << std::endl;
        }
        return ret;
    }

    for (; iter <= maxiter; ++iter)
    {
        const Real rho = dotxy(rh,r);
        if ( rho == 0 )
        {
            ret = 1; break;
        }
        if ( iter == 1 )
        {
            amrex::Copy(p,r,0,0,ncomp,IntVect(nghost));
        }
        else
        {
            const Real beta = (rho/rho_1)*(alpha/omega);
            sxay(p, p, -omega, v, nghost);
            sxay(p, r,   beta, p, nghost);
        }
        amrex::Copy(ph,p,0,0,ncomp,IntVect(nghost));
        Lp.apply(amrlev, mglev, v, ph, MLLinOp::BCMode::Homogeneous, MLLinOp::StateMode::Correction);
        Lp.normalize(amrlev, mglev, v);

        Real rhTv = dotxy(rh,v);
        if ( rhTv != Real(0.0) )
        {
            alpha = rho/rhTv;
        }
        else
        {
            ret = 2; break;
        }
        sxay(sol, sol,  alpha, ph, nghost);
        sxay(s,     r, -alpha,  v, nghost);

        //Subtract mean from s
//        if (Lp.isBottomSingular()) mlmg->makeSolvable(amrlev, mglev, s);

        rnorm = norm_inf(s);

        if ( verbose > 2 && ParallelDescriptor::IOProcessor() )
        {
            amrex::Print() << "MLCGSolver_BiCGStab: Half Iter "
                           << std::setw(11) << iter
                           << " rel. err. "
                           << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;

        amrex::Copy(sh,s,0,0,ncomp,nghost);
        Lp.apply(amrlev, mglev, t, sh, MLLinOp::BCMode::Homogeneous, MLLinOp::StateMode::Correction);
        Lp.normalize(amrlev, mglev, t);
        //
        // This is a little funky.  I want to elide one of the reductions
        // in the following two dotxy()s.  We do that by calculating the "local"
        // values and then reducing the two local values at the same time.
        //
        Real tvals[2] = { dotxy(t,t,true), dotxy(t,s,true) };

        BL_PROFILE_VAR("MLCGSolver::ParallelAllReduce", blp_par);
        ParallelAllReduce::Sum(tvals,2,Lp.BottomCommunicator());
        BL_PROFILE_VAR_STOP(blp_par);

        if ( tvals[0] != Real(0.0) )
        {
            omega = tvals[1]/tvals[0];
        }
        else
        {
            ret = 3; break;
        }
        sxay(sol, sol,  omega, sh, nghost);
        sxay(r,     s, -omega,  t, nghost);

//        if (Lp.isBottomSingular()) mlmg->makeSolvable(amrlev, mglev, r);

        rnorm = norm_inf(r);

        if ( verbose > 2 )
        {
            amrex::Print() << "MLCGSolver_BiCGStab: Iteration "
                           << std::setw(11) << iter
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

    if ( verbose > 0 )
    {
        amrex::Print() << "MLCGSolver_BiCGStab: Final: Iteration "
                       << std::setw(4) << iter
                       << " rel. err. "
                       << rnorm/(rnorm0) << '\n';
    }

    if ( ret == 0 && rnorm > eps_rel*rnorm0 && rnorm > eps_abs)
    {
        if ( verbose > 0 && ParallelDescriptor::IOProcessor() )
            amrex::Warning("MLCGSolver_BiCGStab:: failed to converge!");
        ret = 8;
    }

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rnorm0) )
    {
        amrex::Add(sol,sorig,0,0,ncomp,nghost);
    }
    else
    {
        sol.setVal(0);
        amrex::Add(sol,sorig,0,0,ncomp,nghost);
    }

    return ret;
}

template<class MF>
int
MLCGSolver<MF>::solve_cg (MultiFab&       sol,
                      const MultiFab& rhs,
                      Real            eps_rel,
                      Real            eps_abs)
{
    BL_PROFILE("MLCGSolver::cg");

    const int ncomp = sol.nComp();

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();
    const auto& factory = sol.Factory();

    MultiFab p(ba, dm, ncomp, sol.nGrowVect(), MFInfo(), factory);
    p.setVal(0.0);

    MultiFab sorig(ba, dm, ncomp, nghost, MFInfo(), factory);
    MultiFab r    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MultiFab z    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MultiFab q    (ba, dm, ncomp, nghost, MFInfo(), factory);

    MultiFab::Copy(sorig,sol,0,0,ncomp,nghost);

    Lp.correctionResidual(amrlev, mglev, r, sol, rhs, MLLinOp::BCMode::Homogeneous);

    sol.setVal(0);

    Real       rnorm    = norm_inf(r);
    const Real rnorm0   = rnorm;

    if ( verbose > 0 )
    {
        amrex::Print() << "MLCGSolver_CG: Initial error (error0) :        " << rnorm0 << '\n';
    }

    Real rho_1 = 0;
    int  ret = 0;
    iter = 1;

    if ( rnorm0 == 0 || rnorm0 < eps_abs )
    {
        if ( verbose > 0 ) {
            amrex::Print() << "MLCGSolver_CG: niter = 0,"
                           << ", rnorm = " << rnorm
                           << ", eps_abs = " << eps_abs << std::endl;
        }
        return ret;
    }

    for (; iter <= maxiter; ++iter)
    {
        MultiFab::Copy(z,r,0,0,ncomp,nghost);

        Real rho = dotxy(z,r);

        if ( rho == 0 )
        {
            ret = 1; break;
        }
        if (iter == 1)
        {
            MultiFab::Copy(p,z,0,0,ncomp,nghost);
        }
        else
        {
            Real beta = rho/rho_1;
            sxay(p, z, beta, p, nghost);
        }
        Lp.apply(amrlev, mglev, q, p, MLLinOp::BCMode::Homogeneous, MLLinOp::StateMode::Correction);

        Real alpha;
        Real pw = dotxy(p,q);
        if ( pw != Real(0.0))
        {
            alpha = rho/pw;
        }
        else
        {
            ret = 1; break;
        }

        if ( verbose > 2 )
        {
            amrex::Print() << "MLCGSolver_cg:"
                           << " iter " << iter
                           << " rho " << rho
                           << " alpha " << alpha << '\n';
        }
        sxay(sol, sol, alpha, p, nghost);
        sxay(  r,   r,-alpha, q, nghost);
        rnorm = norm_inf(r);

        if ( verbose > 2 )
        {
            amrex::Print() << "MLCGSolver_cg:       Iteration"
                           << std::setw(4) << iter
                           << " rel. err. "
                           << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) break;

        rho_1 = rho;
    }

    if ( verbose > 0 )
    {
        amrex::Print() << "MLCGSolver_cg: Final Iteration"
                       << std::setw(4) << iter
                       << " rel. err. "
                       << rnorm/(rnorm0) << '\n';
    }

    if ( ret == 0 &&  rnorm > eps_rel*rnorm0 && rnorm > eps_abs )
    {
        if ( verbose > 0 && ParallelDescriptor::IOProcessor() )
            amrex::Warning("MLCGSolver_cg: failed to converge!");
        ret = 8;
    }

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rnorm0) )
    {
        sol.plus(sorig, 0, ncomp, nghost);
    }
    else
    {
        sol.setVal(0);
        sol.plus(sorig, 0, ncomp, nghost);
    }

    return ret;
}

template<class MF>
Real
MLCGSolver<MF>::dotxy (const FabArray<MF>& r, const FabArray<MF>& z, bool local)
{
    BL_PROFILE_VAR_NS("MLCGSolver::ParallelAllReduce", blp_par);
    if (!local) { BL_PROFILE_VAR_START(blp_par); }
    Real result = Lp.xdoty(amrlev, mglev, r, z, local);
    if (!local) { BL_PROFILE_VAR_STOP(blp_par); }
    return result;
}

template<class MF>
Real
MLCGSolver<MF>::norm_inf (const FabArray<MF>& res, bool local)
{
    int ncomp = res.nComp();
    Real result = std::numeric_limits<Real>::lowest();
    for (int n=0; n<ncomp; n++)
      result = std::max(result,Lp.norm0(res,n,0,true));

    if (!local) {
        BL_PROFILE("MLCGSolver::ParallelAllReduce");
        ParallelAllReduce::Max(result, Lp.BottomCommunicator());
    }
    return result;
}

template class MLCGSolver<FArrayBox>;

}
