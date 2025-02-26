
#ifndef AMREX_MLCGSOLVER_H_
#define AMREX_MLCGSOLVER_H_
#include <AMReX_Config.H>

#include <AMReX_MLLinOp.H>

namespace amrex {

template <typename MF>
class MLCGSolverT
{
public:

    using FAB = typename MLLinOpT<MF>::FAB;
    using RT  = typename MLLinOpT<MF>::RT;

    enum struct Type { BiCGStab, CG };

    MLCGSolverT (MLLinOpT<MF>& _lp, Type _typ = Type::BiCGStab);
    ~MLCGSolverT ();

    MLCGSolverT (const MLCGSolverT<MF>& rhs) = delete;
    MLCGSolverT (MLCGSolverT<MF>&& rhs) = delete;
    MLCGSolverT<MF>& operator= (const MLCGSolverT<MF>& rhs) = delete;
    MLCGSolverT<MF>& operator= (MLCGSolverT<MF>&& rhs) = delete;

    void setSolver (Type _typ) noexcept { solver_type = _typ; }

    /**
    * solve the system, Lp(solnL)=rhsL to relative err, tolerance
    * Returns an int indicating success or failure.
    * 0 means success
    * 1 means failed for loss of precision
    * 2 means iterations exceeded
    */
    int solve (MF& solnL, const MF& rhsL, RT eps_rel, RT eps_abs);

    void setVerbose (int _verbose) { verbose = _verbose; }
    [[nodiscard]] int getVerbose () const { return verbose; }

    void setMaxIter (int _maxiter) { maxiter = _maxiter; }
    [[nodiscard]] int getMaxIter () const { return maxiter; }

    void setPrintIdentation (std::string s) { print_ident = std::move(s); }

    /**
    * Is the initial guess provided to the solver zero ?
    * If so, set this to true.
    * The solver will avoid a few operations if this is true.
    * Default is false.
    */
    void setInitSolnZeroed (bool _sol_zeroed) { initial_vec_zeroed = _sol_zeroed; }
    [[nodiscard]] bool getInitSolnZeroed () const { return initial_vec_zeroed; }

    void setNGhost(int _nghost) {nghost = IntVect(_nghost);}
    [[nodiscard]] int getNGhost() {return nghost[0];}

    [[nodiscard]] RT dotxy (const MF& r, const MF& z, bool local = false);
    [[nodiscard]] RT norm_inf (const MF& res, bool local = false);
    int solve_bicgstab (MF& solnL, const MF& rhsL, RT eps_rel, RT eps_abs);
    int solve_cg (MF& solnL, const MF& rhsL, RT eps_rel, RT eps_abs);

    [[nodiscard]] int getNumIters () const noexcept { return iter; }

private:

    MLLinOpT<MF>& Lp;
    Type solver_type;
    const int amrlev = 0;
    const int mglev;
    int verbose   = 0;
    int maxiter   = 100;
    IntVect nghost = IntVect(0);
    int iter = -1;
    bool initial_vec_zeroed = false;
    std::string print_ident;
};

template <typename MF>
MLCGSolverT<MF>::MLCGSolverT (MLLinOpT<MF>& _lp, Type _typ)
    : Lp(_lp), solver_type(_typ), mglev(_lp.NMGLevels(0)-1)
{}

template <typename MF> MLCGSolverT<MF>::~MLCGSolverT () = default;

template <typename MF>
int
MLCGSolverT<MF>::solve (MF& sol, const MF& rhs, RT eps_rel, RT eps_abs)
{
    if (solver_type == Type::BiCGStab) {
        return solve_bicgstab(sol,rhs,eps_rel,eps_abs);
    } else {
        return solve_cg(sol,rhs,eps_rel,eps_abs);
    }
}

template <typename MF>
int
MLCGSolverT<MF>::solve_bicgstab (MF& sol, const MF& rhs, RT eps_rel, RT eps_abs)
{
    BL_PROFILE("MLCGSolver::bicgstab");

    const int ncomp = nComp(sol);

    MF p = Lp.make(amrlev, mglev, nGrowVect(sol));
    MF r = Lp.make(amrlev, mglev, nGrowVect(sol));
    setVal(p, RT(0.0)); // Make sure all entries are initialized to avoid errors
    setVal(r, RT(0.0));

    MF rh    = Lp.make(amrlev, mglev, nghost);
    MF v     = Lp.make(amrlev, mglev, nghost);
    MF t     = Lp.make(amrlev, mglev, nghost);


    MF sorig;

    if ( initial_vec_zeroed ) {
        LocalCopy(r,rhs,0,0,ncomp,nghost);
    } else {
        sorig = Lp.make(amrlev, mglev, nghost);

        Lp.correctionResidual(amrlev, mglev, r, sol, rhs, MLLinOpT<MF>::BCMode::Homogeneous);

        LocalCopy(sorig,sol,0,0,ncomp,nghost);
        setVal(sol, RT(0.0));
    }

    // Then normalize
    Lp.normalize(amrlev, mglev, r);
    LocalCopy(rh, r, 0,0,ncomp,nghost);

    RT rnorm = norm_inf(r);
    const RT rnorm0 = rnorm;

    if ( verbose > 0 )
    {
        amrex::Print() << print_ident << "MLCGSolver_BiCGStab: Initial error (error0) =        " << rnorm0 << '\n';
    }
    int ret = 0;
    iter = 1;
    RT rho_1 = 0, alpha = 0, omega = 0;

    if ( rnorm0 == 0 || rnorm0 < eps_abs )
    {
        if ( verbose > 0 )
        {
            amrex::Print() << print_ident << "MLCGSolver_BiCGStab: niter = 0,"
                           << ", rnorm = " << rnorm
                           << ", eps_abs = " << eps_abs << '\n';
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
            LocalCopy(p,r,0,0,ncomp,nghost);
        }
        else
        {
            const RT beta = (rho/rho_1)*(alpha/omega);
            Saxpy(p, -omega, v, 0, 0, ncomp, nghost); // p += -omega*v
            Xpay(p, beta, r, 0, 0, ncomp, nghost); // p = r + beta*p
        }
        Lp.apply(amrlev, mglev, v, p, MLLinOpT<MF>::BCMode::Homogeneous, MLLinOpT<MF>::StateMode::Correction);
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
        Saxpy(sol, alpha, p, 0, 0, ncomp, nghost); // sol += alpha * p
        Saxpy(r,  -alpha, v, 0, 0, ncomp, nghost); // r += -alpha * v

        rnorm = norm_inf(r);

        if ( verbose > 2 && ParallelDescriptor::IOProcessor() )
        {
            amrex::Print() << print_ident << "MLCGSolver_BiCGStab: Half Iter "
                           << std::setw(11) << iter
                           << " rel. err. "
                           << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) { break; }

        Lp.apply(amrlev, mglev, t, r, MLLinOpT<MF>::BCMode::Homogeneous, MLLinOpT<MF>::StateMode::Correction);
        Lp.normalize(amrlev, mglev, t);
        //
        // This is a little funky.  I want to elide one of the reductions
        // in the following two dotxy()s.  We do that by calculating the "local"
        // values and then reducing the two local values at the same time.
        //
        RT tvals[2] = { dotxy(t,t,true), dotxy(t,r,true) };

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
        Saxpy(sol, omega, r, 0, 0, ncomp, nghost); // sol += omega * r
        Saxpy(r,  -omega, t, 0, 0, ncomp, nghost); // r += -omega * t

        rnorm = norm_inf(r);

        if ( verbose > 2 )
        {
            amrex::Print() << print_ident << "MLCGSolver_BiCGStab: Iteration "
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
        amrex::Print() << print_ident << "MLCGSolver_BiCGStab: Final: Iteration "
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
        if ( !initial_vec_zeroed ) {
            LocalAdd(sol, sorig, 0, 0, ncomp, nghost);
        }
        if (ret == 8) { ret = 9; }
    }
    else
    {
        setVal(sol, RT(0.0));
        if ( !initial_vec_zeroed ) {
            LocalAdd(sol, sorig, 0, 0, ncomp, nghost);
        }
    }

    return ret;
}

template <typename MF>
int
MLCGSolverT<MF>::solve_cg (MF& sol, const MF& rhs, RT eps_rel, RT eps_abs)
{
    BL_PROFILE("MLCGSolver::cg");

    const int ncomp = nComp(sol);

    MF p = Lp.make(amrlev, mglev, nGrowVect(sol));
    setVal(p, RT(0.0));

    MF r     = Lp.make(amrlev, mglev, nghost);
    MF q     = Lp.make(amrlev, mglev, nghost);

    MF sorig;

    if ( initial_vec_zeroed ) {
        LocalCopy(r,rhs,0,0,ncomp,nghost);
    } else {
        sorig = Lp.make(amrlev, mglev, nghost);

        Lp.correctionResidual(amrlev, mglev, r, sol, rhs, MLLinOpT<MF>::BCMode::Homogeneous);

        LocalCopy(sorig,sol,0,0,ncomp,nghost);
        setVal(sol, RT(0.0));
    }

    RT       rnorm    = norm_inf(r);
    const RT rnorm0   = rnorm;

    if ( verbose > 0 )
    {
        amrex::Print() << print_ident << "MLCGSolver_CG: Initial error (error0) :        " << rnorm0 << '\n';
    }

    RT rho_1 = 0;
    int  ret = 0;
    iter = 1;

    if ( rnorm0 == 0 || rnorm0 < eps_abs )
    {
        if ( verbose > 0 ) {
            amrex::Print() << print_ident << "MLCGSolver_CG: niter = 0,"
                           << ", rnorm = " << rnorm
                           << ", eps_abs = " << eps_abs << '\n';
        }
        return ret;
    }

    for (; iter <= maxiter; ++iter)
    {
        RT rho = dotxy(r,r);

        if ( rho == 0 )
        {
            ret = 1; break;
        }
        if (iter == 1)
        {
            LocalCopy(p,r,0,0,ncomp,nghost);
        }
        else
        {
            RT beta = rho/rho_1;
            Xpay(p, beta, r, 0, 0, ncomp, nghost); // p = r + beta * p
        }
        Lp.apply(amrlev, mglev, q, p, MLLinOpT<MF>::BCMode::Homogeneous, MLLinOpT<MF>::StateMode::Correction);

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
            amrex::Print() << print_ident << "MLCGSolver_cg:"
                           << " iter " << iter
                           << " rho " << rho
                           << " alpha " << alpha << '\n';
        }
        Saxpy(sol, alpha, p, 0, 0, ncomp, nghost); // sol += alpha * p
        Saxpy(r, -alpha, q, 0, 0, ncomp, nghost); // r += -alpha * q
        rnorm = norm_inf(r);

        if ( verbose > 2 )
        {
            amrex::Print() << print_ident << "MLCGSolver_cg:       Iteration"
                           << std::setw(4) << iter
                           << " rel. err. "
                           << rnorm/(rnorm0) << '\n';
        }

        if ( rnorm < eps_rel*rnorm0 || rnorm < eps_abs ) { break; }

        rho_1 = rho;
    }

    if ( verbose > 0 )
    {
        amrex::Print() << print_ident << "MLCGSolver_cg: Final Iteration"
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
        if ( !initial_vec_zeroed ) {
            LocalAdd(sol, sorig, 0, 0, ncomp, nghost);
        }
        if (ret == 8) { ret = 9; }
    }
    else
    {
        setVal(sol, RT(0.0));
        if ( !initial_vec_zeroed ) {
            LocalAdd(sol, sorig, 0, 0, ncomp, nghost);
        }
    }

    return ret;
}

template <typename MF>
auto
MLCGSolverT<MF>::dotxy (const MF& r, const MF& z, bool local) -> RT
{
    BL_PROFILE_VAR_NS("MLCGSolver::ParallelAllReduce", blp_par);
    if (!local) { BL_PROFILE_VAR_START(blp_par); }
    RT result = Lp.xdoty(amrlev, mglev, r, z, local);
    if (!local) { BL_PROFILE_VAR_STOP(blp_par); }
    return result;
}

template <typename MF>
auto
MLCGSolverT<MF>::norm_inf (const MF& res, bool local) -> RT
{
    int ncomp = nComp(res);
    RT result = norminf(res,0,ncomp,IntVect(0),true);
    if (!local) {
        BL_PROFILE("MLCGSolver::ParallelAllReduce");
        ParallelAllReduce::Max(result, Lp.BottomCommunicator());
    }
    return result;
}

using MLCGSolver = MLCGSolverT<MultiFab>;

}

#endif /*_CGSOLVER_H_*/
