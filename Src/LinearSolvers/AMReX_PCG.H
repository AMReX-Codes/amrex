#ifndef AMREX_PCG_H_
#define AMREX_PCG_H_
#include <AMReX_Config.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>

#include <cmath>
#include <iomanip>
#include <limits>

namespace amrex {

/**
 * \brief Preconditioned conjugate gradient solver.
 *
 * Requires a symmetric positive definite operator and preconditioner.
 * The operator type `M` must provide the same interface as for
 * amrex::GMRES: `makeVecRHS`, `makeVecLHS`, `norm2`, `dotProduct`, `assign`,
 * `increment`, `linComb`, `scale`, `setToZero`, `apply` and `precond`.
 */
template <typename V, typename M>
class PCG
{
public:

    using RT = typename M::RT; // double or float

    PCG () = default;

    //! Defines with a reference to M. It's the user's responsibility to
    //! keep the M object alive for PCG to be functional. This function
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
     * 0 converged, 1 max iterations reached, 2 p.Ap <= 0 (A not positive
     * definite), 3 r.z <= 0 (preconditioner not positive definite). A NaN
     * or infinity counts as the matching breakdown.
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
void PCG<V,M>::solve (V& a_sol, V const& a_rhs, RT a_tol_rel, RT a_tol_abs, int a_its)
{
    BL_PROFILE("PCG::solve()");

    AMREX_ALWAYS_ASSERT(m_linop != nullptr);

    m_rtol = a_tol_rel;
    m_atol = a_tol_abs;
    if (a_its < 0) { a_its = m_maxiter; }

    m_its = 0;
    m_status = -1;

    V r = m_linop->makeVecRHS();
    V z = m_linop->makeVecLHS();
    V p = m_linop->makeVecLHS();
    V q = m_linop->makeVecRHS();

    if (m_nonzero_guess) {
        // r = b - A x. apply() needs an x made by makeVecLHS.
        m_linop->assign(p, a_sol);
        m_linop->apply(q, p);
        m_linop->linComb(r, RT(1), a_rhs, RT(-1), q);
    } else {
        m_linop->setToZero(a_sol);
        m_linop->assign(r, a_rhs);
    }

    RT rnorm = m_linop->norm2(r);
    RT const rnorm0 = rnorm;
    m_res = rnorm;

    if (m_verbose > 0) {
        amrex::Print() << "PCG: Initial residual (error0) = " << rnorm0 << '\n';
    }

    if (rnorm0 == RT(0) || rnorm0 < m_atol) {
        m_status = 0;
        return;
    }

    auto converged = [&] (RT rn) {
        return rn == RT(0) || rn < m_rtol*rnorm0 || rn < m_atol;
    };

    m_linop->precond(z, r);
    m_linop->assign(p, z);
    RT rz = m_linop->dotProduct(r, z);
    if (rz <= RT(0) || !std::isfinite(rz)) { m_status = 3; }

    for (int iter = 1; iter <= a_its && m_status == -1; ++iter)
    {
        m_linop->apply(q, p);
        RT const pq = m_linop->dotProduct(p, q);
        if (pq <= RT(0) || !std::isfinite(pq)) { m_status = 2; break; }
        RT const alpha = rz/pq;

        m_linop->increment(a_sol, p, alpha);
        m_linop->increment(r, q, -alpha);

        rnorm = m_linop->norm2(r);
        m_its = iter;
        m_res = rnorm;

        if (m_verbose > 1) {
            amrex::Print() << "PCG: Iteration " << std::setw(11) << iter
                           << " rel. err. " << rnorm/rnorm0 << '\n';
        }

        if (converged(rnorm)) { m_status = 0; break; }
        if (!std::isfinite(rnorm)) { m_status = 3; break; }
        if (iter == a_its) { break; }

        m_linop->precond(z, r);
        RT const rz_new = m_linop->dotProduct(r, z);
        if (rz_new <= RT(0) || !std::isfinite(rz_new)) { m_status = 3; break; }
        RT const beta = rz_new/rz;
        rz = rz_new;

        m_linop->scale(p, beta);
        m_linop->increment(p, z, RT(1));         // p = z + beta*p
    }

    if (m_status == -1) { m_status = 1; }

    if (m_verbose > 0) {
        amrex::Print() << "PCG: Final: Iteration " << std::setw(4) << m_its
                       << " rel. err. " << rnorm/rnorm0 << '\n';
    }

    if (m_status != 0 && m_verbose > 0) {
        amrex::Print() << "PCG: Failed to converge, status = " << m_status << '\n';
    }
}

}

#endif
