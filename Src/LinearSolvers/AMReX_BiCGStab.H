#ifndef AMREX_BICGSTAB_H_
#define AMREX_BICGSTAB_H_
#include <AMReX_Config.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>

#include <cmath>
#include <iomanip>
#include <limits>

namespace amrex {

/**
 * \brief BiCGStab solver with right preconditioning.
 *
 * The operator type `M` must provide the same interface as for
 * amrex::GMRES: `makeVecRHS`, `makeVecLHS`, `norm2`, `dotProduct`, `assign`,
 * `increment`, `linComb`, `scale`, `setToZero`, `apply` and `precond`.
 */
template <typename V, typename M>
class BiCGStab
{
public:

    using RT = typename M::RT; // double or float

    BiCGStab () = default;

    //! Defines with a reference to M. It's the user's responsibility to
    //! keep the M object alive for BiCGStab to be functional. This function
    //! must be called before solve() can be called.
    void define (M& linop) { m_linop = &linop; }

    /**
     * \brief Solve the linear system
     *
     * The initial guess is zero unless setInitialGuessNonzero(true) was
     * called, in which case it is the value of \p a_sol on entry.
     *
     * \param a_sol     unknowns, i.e., x in A x = b.
     * \param a_rhs     RHS, i.e., b in A x = b.
     * \param a_tol_rel relative tolerance, relative to the initial residual
     *                  norm, which is ||b|| for a zero initial guess.
     * \param a_tol_abs absolute tolerance.
     * \param a_its     optional argument specifying the maximum number of iterations.
     */
    void solve (V& a_sol, V const& a_rhs, RT a_tol_rel, RT a_tol_abs, int a_its=-1);

    //! Use \p a_sol passed to solve() as the initial guess (default: false,
    //! i.e., start from zero).
    void setInitialGuessNonzero (bool flag) { m_nonzero_guess = flag; }

    //! Sets verbosity.
    void setVerbose (int v) { m_verbose = v; }

    //! Sets the max number of iterations
    void setMaxIters (int niters) { m_maxiter = niters; }

    //! Gets the number of iterations.
    [[nodiscard]] int getNumIters () const { return m_its; }

    /**
     * \brief Gets the solver status.
     *
     * 0 converged, 1 max iterations reached, 2 rho = 0, 3 rh.v = 0,
     * 4 t.t = 0, 5 omega = 0. A NaN or infinity counts as the matching
     * breakdown.
     */
    [[nodiscard]] int getStatus () const { return m_status; }

    //! Gets the 2-norm of the residual.
    [[nodiscard]] RT getResidualNorm () const { return m_res; }

private:

    int m_verbose = 0;
    bool m_nonzero_guess = false;
    int m_maxiter = 2000;
    int m_its = 0;
    int m_status = -1;
    RT m_res = std::numeric_limits<RT>::max();
    RT m_rtol = RT(0);
    RT m_atol = RT(0);

    M* m_linop = nullptr;
};

template <typename V, typename M>
void BiCGStab<V,M>::solve (V& a_sol, V const& a_rhs, RT a_tol_rel, RT a_tol_abs, int a_its)
{
    BL_PROFILE("BiCGStab::solve()");

    AMREX_ALWAYS_ASSERT(m_linop != nullptr);

    m_rtol = a_tol_rel;
    m_atol = a_tol_abs;
    if (a_its < 0) { a_its = m_maxiter; }

    m_its = 0;
    m_status = -1;

    V r  = m_linop->makeVecRHS();
    V rh = m_linop->makeVecRHS();
    V p  = m_linop->makeVecRHS();
    V ph = m_linop->makeVecLHS();
    V v  = m_linop->makeVecRHS();
    V sh = m_linop->makeVecLHS();
    V t  = m_linop->makeVecRHS();

    if (m_nonzero_guess) {
        // r = b - A x. apply() needs an x made by makeVecLHS.
        m_linop->assign(ph, a_sol);
        m_linop->apply(v, ph);
        m_linop->linComb(r, RT(1), a_rhs, RT(-1), v);
    } else {
        m_linop->setToZero(a_sol);
        m_linop->assign(r, a_rhs);
    }
    m_linop->assign(rh, r);

    RT rnorm = m_linop->norm2(r);
    RT const rnorm0 = rnorm;
    m_res = rnorm;

    if (m_verbose > 0) {
        amrex::Print() << "BiCGStab: Initial residual (error0) = " << rnorm0 << '\n';
    }

    if (rnorm0 == RT(0) || rnorm0 < m_atol) {
        m_status = 0;
        return;
    }

    RT rho_1 = RT(0), alpha = RT(0), omega = RT(0);

    auto converged = [&] (RT rn) {
        return rn == RT(0) || rn < m_rtol*rnorm0 || rn < m_atol;
    };

    for (int iter = 1; iter <= a_its; ++iter)
    {
        RT const rho = m_linop->dotProduct(rh, r);
        if (rho == RT(0) || !std::isfinite(rho)) { m_status = 2; break; }

        if (iter == 1) {
            m_linop->assign(p, r);
        } else {
            RT const beta = (rho/rho_1)*(alpha/omega);
            m_linop->increment(p, v, -omega);    // p -= omega*v
            m_linop->scale(p, beta);
            m_linop->increment(p, r, RT(1));     // p = r + beta*(p - omega*v)
        }

        m_linop->precond(ph, p);
        m_linop->apply(v, ph);

        RT const rhv = m_linop->dotProduct(rh, v);
        if (rhv == RT(0) || !std::isfinite(rhv)) { m_status = 3; break; }
        alpha = rho/rhv;

        m_linop->increment(a_sol, ph, alpha);
        m_linop->increment(r, v, -alpha);        // r is now s

        rnorm = m_linop->norm2(r);
        m_its = iter;
        m_res = rnorm;

        if (m_verbose > 1) {
            amrex::Print() << "BiCGStab: Half Iter " << std::setw(11) << iter
                           << " rel. err. " << rnorm/rnorm0 << '\n';
        }

        if (converged(rnorm)) { m_status = 0; break; }

        m_linop->precond(sh, r);
        m_linop->apply(t, sh);

        RT const tt = m_linop->dotProduct(t, t);
        if (tt == RT(0) || !std::isfinite(tt)) { m_status = 4; break; }
        omega = m_linop->dotProduct(t, r) / tt;
        if (omega == RT(0) || !std::isfinite(omega)) { m_status = 5; break; }

        m_linop->increment(a_sol, sh, omega);
        m_linop->increment(r, t, -omega);

        rnorm = m_linop->norm2(r);
        m_res = rnorm;

        if (m_verbose > 1) {
            amrex::Print() << "BiCGStab: Iteration " << std::setw(11) << iter
                           << " rel. err. " << rnorm/rnorm0 << '\n';
        }

        if (converged(rnorm)) { m_status = 0; break; }

        rho_1 = rho;
    }

    if (m_status == -1) { m_status = 1; }

    if (m_verbose > 0) {
        amrex::Print() << "BiCGStab: Final: Iteration " << std::setw(4) << m_its
                       << " rel. err. " << rnorm/rnorm0 << '\n';
    }

    if (m_status != 0 && m_verbose > 0) {
        amrex::Print() << "BiCGStab: Failed to converge, status = " << m_status << '\n';
    }
}

}

#endif
