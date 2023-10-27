#include <AMReX_Config.H>

#include <AMReX_MLLinOp.H>
#include <TestMFS_MLLinOp.H>


namespace amrex {

// template <typename MFSet>
class MLCGSolver_Test
{
public:

    using MFSet = typename amrex::MultiFabSet;
    using FAB = typename MFSet::fab_type;
    using RT  = typename MFSet::value_type;

    enum struct Type { BiCGStab, CG };
    MLCGSolver_Test (MLLinOp_Test& _lp, Type _typ = Type::BiCGStab);
    ~MLCGSolver_Test ();

    MLCGSolver_Test (const MLCGSolver_Test& rhs) = delete;
    MLCGSolver_Test (MLCGSolver_Test&& rhs) = delete;
    MLCGSolver_Test& operator= (const MLCGSolver_Test& rhs) = delete;
    MLCGSolver_Test& operator= (MLCGSolver_Test&& rhs) = delete;

    void setSolver (Type _typ) noexcept { solver_type = _typ; }

    /**
    * solve the system, Lp(solnL)=rhsL to relative err, tolerance
    * Returns an int indicating success or failure.
    * 0 means success
    * 1 means failed for loss of precision
    * 2 means iterations exceeded
    */
    int solve (MFSet& solnL, const MFSet& rhsL, RT eps_rel, RT eps_abs);

    void setVerbose (int _verbose) { verbose = _verbose; }
    [[nodiscard]] int getVerbose () const { return verbose; }

    void setMaxIter (int _maxiter) { maxiter = _maxiter; }
    [[nodiscard]] int getMaxIter () const { return maxiter; }

    void setNGhost(int _nghost) {nghost = IntVect(_nghost);}
    [[nodiscard]] int getNGhost() {return nghost[0];}

    [[nodiscard]] RT dotxy (const MFSet& r, const MFSet& z, bool local = false);
    [[nodiscard]] RT norm_inf (const MFSet& res, bool local = false);
    int solve_bicgstab (MFSet& solnL, const MFSet& rhsL, RT eps_rel, RT eps_abs);
    int solve_cg (MFSet& solnL, const MFSet& rhsL, RT eps_rel, RT eps_abs);

    [[nodiscard]] int getNumIters () const noexcept { return iter; }

private:

    MLLinOp_Test& Lp;
    Type solver_type;
    const int amrlev = 0;
    const int mglev;
    int verbose   = 0;
    int maxiter   = 100;
    IntVect nghost = IntVect(0);
    int iter = -1;
};

// template <typename MFSet>
MLCGSolver_Test::MLCGSolver_Test (MLLinOp_Test& _lp, Type _typ)
    : Lp(_lp), solver_type(_typ), mglev(_lp.NMGLevels(0)-1)
{}

// template <typename MFSet> MLCGSolver_Test::~MLCGSolver_Test () = default;

// template <typename MFSet>
int
MLCGSolver_Test::solve (MFSet& sol, const MFSet& rhs, RT eps_rel, RT eps_abs)
{
    if (solver_type == Type::BiCGStab) {
        return solve_bicgstab(sol,rhs,eps_rel,eps_abs);
    } else {
        return solve_cg(sol,rhs,eps_rel,eps_abs);
    }
}

// template <typename MFSet>
int
MLCGSolver_Test::solve_bicgstab (MFSet& sol, const MFSet& rhs, RT eps_rel, RT eps_abs)
{
    BL_PROFILE("MLCGSolver::bicgstab");

    const int ncomp = sol.nComp();

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();
    const auto& factory = sol.Factory();

    MFSet ph(ba, dm, ncomp, sol.nGrowVect(), MFInfo(), factory);
    MFSet sh(ba, dm, ncomp, sol.nGrowVect(), MFInfo(), factory);
    ph.setVal(RT(0.0));
    sh.setVal(RT(0.0));

    MFSet sorig(ba, dm, ncomp, nghost, MFInfo(), factory);
    MFSet p    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MFSet r    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MFSet s    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MFSet rh   (ba, dm, ncomp, nghost, MFInfo(), factory);
    MFSet v    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MFSet t    (ba, dm, ncomp, nghost, MFInfo(), factory);

    Lp.correctionResidual(amrlev, mglev, r, sol, rhs, MLLinOp_Test::BCMode::Homogeneous);

    // Then normalize
    Lp.normalize(amrlev, mglev, r);

    sorig.LocalCopy(sol,0,0,ncomp,nghost);
    rh.LocalCopy   (r  ,0,0,ncomp,nghost);

    sol.setVal(RT(0.0));

    RT rnorm = norm_inf(r);
    const RT rnorm0 = rnorm;

    if ( verbose > 0 )
    {
        amrex::Print() << "MLCGSolver_BiCGStab: Initial error (error0) =        " << rnorm0 << '\n';
    }
    int ret = 0;
    iter = 1;
    RT rho_1 = 0, alpha = 0, omega = 0;

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
        const RT rho = dotxy(rh,r);
        if ( rho == 0 )
        {
            ret = 1; break;
        }
        if ( iter == 1 )
        {
            p.LocalCopy(r,0,0,ncomp,nghost);
        }
        else
        {
            const RT beta = (rho/rho_1)*(alpha/omega);
            MFSet::Saxpy(p, -omega, v, 0, 0, ncomp, nghost); // p += -omega*v
            MFSet::Xpay(p, beta, r, 0, 0, ncomp, nghost); // p = r + beta*p
        }
        ph.LocalCopy(p,0,0,ncomp,nghost);
        Lp.apply(amrlev, mglev, v, ph, MLLinOp_Test::BCMode::Homogeneous, MLLinOp_Test::StateMode::Correction);
        Lp.normalize(amrlev, mglev, v);

        RT rhTv = dotxy(rh,v);
        if ( rhTv != RT(0.0) )
        {
            alpha = rho/rhTv;
        }
        else
        {
            ret = 2; break;
        }
        MFSet::Saxpy(sol, alpha, ph, 0, 0, ncomp, nghost); // sol += alpha * ph
        MFSet::LinComb(s, RT(1.0), r, 0, -alpha, v, 0, 0, ncomp, nghost); // s = r - alpha * v

        rnorm = norm_inf(s);

        if ( verbose > 2 && ParallelDescriptor::IOProcessor() )
        {
            amrex::Print() << "MLCGSolver_BiCGStab: Half Iter "
                           << std::setw(11) << iter
                           << " rel. err. "
                           << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) { break; }

        sh.LocalCopy(s,0,0,ncomp,nghost);
        Lp.apply(amrlev, mglev, t, sh, MLLinOp_Test::BCMode::Homogeneous, MLLinOp_Test::StateMode::Correction);
        Lp.normalize(amrlev, mglev, t);
        //
        // This is a little funky.  I want to elide one of the reductions
        // in the following two dotxy()s.  We do that by calculating the "local"
        // values and then reducing the two local values at the same time.
        //
        RT tvals[2] = { dotxy(t,t,true), dotxy(t,s,true) };

        BL_PROFILE_VAR("MLCGSolver::ParallelAllReduce", blp_par);
        ParallelAllReduce::Sum(tvals,2,Lp.BottomCommunicator());
        BL_PROFILE_VAR_STOP(blp_par);

        if ( tvals[0] != RT(0.0) )
        {
            omega = tvals[1]/tvals[0];
        }
        else
        {
            ret = 3; break;
        }
        MFSet::Saxpy(sol, omega, sh, 0, 0, ncomp, nghost); // sol += omega * sh
        MFSet::LinComb(r, RT(1.0), s, 0, -omega, t, 0, 0, ncomp, nghost); // r = s - omega * t

        rnorm = norm_inf(r);

        if ( verbose > 2 )
        {
            amrex::Print() << "MLCGSolver_BiCGStab: Iteration "
                           << std::setw(11) << iter
                           << " rel. err. "
                           << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) { break; }

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
        if ( verbose > 0 && ParallelDescriptor::IOProcessor() ) {
            amrex::Warning("MLCGSolver_BiCGStab:: failed to converge!");
        }
        ret = 8;
    }

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rnorm0) )
    {
        sol.LocalAdd(sorig, 0, 0, ncomp, nghost);
    }
    else
    {
        sol.setVal(RT(0.0));
        sol.LocalAdd(sorig, 0, 0, ncomp, nghost);
    }

    return ret;
}

// template <typename MFSet>
int
MLCGSolver_Test::solve_cg (MFSet& sol, const MFSet& rhs, RT eps_rel, RT eps_abs)
{
    BL_PROFILE("MLCGSolver::cg");

    const int ncomp = sol.nComp();

    const BoxArray& ba = sol.boxArray();
    const DistributionMapping& dm = sol.DistributionMap();
    const auto& factory = sol.Factory();

    MFSet p(ba, dm, ncomp, sol.nGrowVect(), MFInfo(), factory);
    p.setVal(RT(0.0));

    MFSet sorig(ba, dm, ncomp, nghost, MFInfo(), factory);
    MFSet r    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MFSet z    (ba, dm, ncomp, nghost, MFInfo(), factory);
    MFSet q    (ba, dm, ncomp, nghost, MFInfo(), factory);

    sorig.LocalCopy(sol,0,0,ncomp,nghost);

    Lp.correctionResidual(amrlev, mglev, r, sol, rhs, MLLinOp_Test::BCMode::Homogeneous);

    sol.setVal(RT(0.0));

    RT       rnorm    = norm_inf(r);
    const RT rnorm0   = rnorm;

    if ( verbose > 0 )
    {
        amrex::Print() << "MLCGSolver_CG: Initial error (error0) :        " << rnorm0 << '\n';
    }

    RT rho_1 = 0;
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
        z.LocalCopy(r,0,0,ncomp,nghost);

        RT rho = dotxy(z,r);

        if ( rho == 0 )
        {
            ret = 1; break;
        }
        if (iter == 1)
        {
            p.LocalCopy(z,0,0,ncomp,nghost);
        }
        else
        {
            RT beta = rho/rho_1;
            MFSet::Xpay(p, beta, z, 0, 0, ncomp, nghost); // p = z + beta * p
        }
        Lp.apply(amrlev, mglev, q, p, MLLinOp_Test::BCMode::Homogeneous, MLLinOp_Test::StateMode::Correction);

        RT alpha;
        RT pw = dotxy(p,q);
        if ( pw != RT(0.0))
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
        MFSet::Saxpy(sol, alpha, p, 0, 0, ncomp, nghost); // sol += alpha * p
        MFSet::Saxpy(r, -alpha, q, 0, 0, ncomp, nghost); // r += -alpha * q
        rnorm = norm_inf(r);

        if ( verbose > 2 )
        {
            amrex::Print() << "MLCGSolver_cg:       Iteration"
                           << std::setw(4) << iter
                           << " rel. err. "
                           << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) { break; }

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
        if ( verbose > 0 && ParallelDescriptor::IOProcessor() ) {
            amrex::Warning("MLCGSolver_cg: failed to converge!");
        }
        ret = 8;
    }

    if ( ( ret == 0 || ret == 8 ) && (rnorm < rnorm0) )
    {
        sol.LocalAdd(sorig, 0, 0, ncomp, nghost);
    }
    else
    {
        sol.setVal(RT(0.0));
        sol.LocalAdd(sorig, 0, 0, ncomp, nghost);
    }

    return ret;
}

// template <typename MFSet>
auto
MLCGSolver_Test::dotxy (const MFSet& r, const MFSet& z, bool local) -> RT
{
    BL_PROFILE_VAR_NS("MLCGSolver::ParallelAllReduce", blp_par);
    if (!local) { BL_PROFILE_VAR_START(blp_par); }
    RT result = Lp.xdoty(amrlev, mglev, r, z, local);
    if (!local) { BL_PROFILE_VAR_STOP(blp_par); }
    return result;
}

// template <typename MFSet>
auto
MLCGSolver_Test::norm_inf (const MFSet& res, bool local) -> RT
{
    int ncomp = res.nComp();
    RT result = res.norminf(0,ncomp,IntVect(0),true);
    if (!local) {
        BL_PROFILE("MLCGSolver::ParallelAllReduce");
        ParallelAllReduce::Max(result, Lp.BottomCommunicator());
    }
    return result;
}

// using MLCGSolver = MLCGSolver_Test<MultiFab>;

}
