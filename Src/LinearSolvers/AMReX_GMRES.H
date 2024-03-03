#ifndef AMREX_GMRES_H_
#define AMREX_GMRES_H_
#include <AMReX_Config.H>

#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>
#include <AMReX_TableData.H>
#include <AMReX_Vector.H>
#include <cmath>
#include <limits>
#include <memory>

namespace amrex {

/**
 * \brief GMRES
 *
 * This class implements the GMRES algorithm. The template parameter V is
 * for a linear algebra vector class. For example, it could be
 * amrex::MultiFab. The other template parameter M is for a linear operator
 * class with a number of required member functions. Note that conceptually
 * M contains a matrix. However, it does not mean it needs to have a member
 * variable storing the matrix, because GMRES only needs the matrix vector
 * product, not the matrix itself.
 *
 * \tparam V linear algebra vector. It must be default constructible, move
 *           constructible, and move assignable.
 * \tparam M linear operator. A list of required member functions for M is
 *           shown below. Here RT (typename V::value_type) is either double
 *           or float.
 *             - void apply(V& lhs, V const& rhs)\n
 *               lhs = L(rhs), where L is the linear operator performing matrix
 *               vector product.
 *             - void assign(V& lhs, V const& rhs)\n
 *               lhs = rhs.
 *             - RT dotProduct(V const& v1, V const& v2)\n
 *               returns v1 * v2.
 *             - void increment(V& lhs, V const& rhs, RT a)\n
 *               lhs += a * rhs.
 *             - void linComb(V& lhs, RT a, V const& rhs_a, RT b, V const& rhs_b)\n
 *               lhs = a * rhs_a + b * rhs_b.
 *             - V makeVecRHS()\n
 *               returns a V object that is suitable as RHS in M x = b. The reason
 *               we distinguish between LHS and RHS is M might need the distinction
 *               for efficiency. For example, if V is MultiFab, we might need the x
 *               in the LHS of M x = b to have ghost cells for efficiency, whereas
 *               no ghost cells are needed for the RHS (i.e., b).
 *             - V makeVecLHS()\n
 *               returns a V object that is suitable as LHS in M x = b. See the
 *               description for makeVecRHS for more details.
 *             - RT norm2(V const& v)\n
 *               returns the 2-norm of v.
 *             - void precond(V& lhs, V const& rhs)\n
 *               applies preconditioner to rhs. If there is no preconditioner,
 *               this function should do lhs = rhs.
 *             - void setToZero(V& v)\n
 *               v = 0.
 */
template <typename V, typename M>
class GMRES
{
public:

    using RT = typename M::RT; // double or float

    GMRES ();

    //! Defines with a reference to M. It's the user's responsibility to
    //! keep the M object alive for GMRES to be functional. This function
    //! must be called before solve() can be called.
    void define (M& linop);

    /**
     * \brief Solve the linear system
     *
     * \param a_sol     unknowns, i.e., x in A x = b.
     * \param a_rhs     RHS, i.e., b in A x = b.
     * \param a_tol_rel relative tolerance.
     * \param a_tol_abs absolute tolerance.
     * \patam a_its     optional argument specifying the maximum number of iterations.
     */
    void solve (V& a_sol, V const& a_rhs, RT a_tol_rel, RT a_tol_abs, int a_its=-1);

    //! Sets verbosity.
    void setVerbose (int v) { m_verbose = v; }

    //! Sets restart length. The default is 30.
    void setRestartLength (int rl);

    //! Sets the max number of iterations
    void setMaxIters (int niters) { m_maxiter = niters; }

    //! Gets the number of iterations.
    [[nodiscard]] int getNumIters () const { return m_its; }

    //! Gets the solver status.
    [[nodiscard]] int getStatus () const { return m_status; }

    //! Gets the 2-norm of the residual.
    [[nodiscard]] RT getResidualNorm () const { return m_res; }

private:
    void clear ();
    void allocate_scratch ();
    void cycle (V& a_xx, int& a_status, int& a_itcount, RT& a_rnorm0);
    void build_solution (V& a_xx, int it);
    void compute_residual (V& a_rr, V const& a_xx, V const& a_bb);

    bool converged (RT r0, RT r) const;

    void gram_schmidt_orthogonalization (int it);
    void update_hessenberg (int it, bool happyend, RT& res);

    int m_verbose = 0;
    int m_maxiter = 2000;
    int m_its = 0;
    int m_status = -1;
    int m_restrtlen = 30;
    RT m_res = std::numeric_limits<RT>::max();
    RT m_rtol = RT(0);
    RT m_atol = RT(0);
    Vector<RT> m_hh_1d;
    Vector<RT> m_hes_1d;
    Table2D<RT> m_hh;
    Table2D<RT> m_hes;
    Vector<RT> m_grs;
    Vector<RT> m_cc;
    Vector<RT> m_ss;
    std::unique_ptr<V> m_v_tmp_rhs;
    std::unique_ptr<V> m_v_tmp_lhs;
    Vector<V> m_vv;
    M* m_linop = nullptr;
};

template <typename V, typename M>
GMRES<V,M>::GMRES ()
{
    allocate_scratch();
}

template <typename V, typename M>
void GMRES<V,M>::allocate_scratch ()
{
    int rs = m_restrtlen;

    m_hh_1d.resize(std::size_t(rs + 2) * (rs + 1));
    m_hh = Table2D<RT>(m_hh_1d.data(), {0,0}, {rs+1,rs}); // (0:rs+1,0:rs)

    m_hes_1d.resize(std::size_t(rs + 2) * (rs + 1));
    m_hes = Table2D<RT>(m_hes_1d.data(), {0,0}, {rs+1,rs}); // (0:rs+1,0:rs)

    m_grs.resize(rs + 2);
    m_cc.resize(rs + 1);
    m_ss.resize(rs + 1);
}

template <typename V, typename M>
void GMRES<V,M>::setRestartLength (int rl)
{
    if (m_restrtlen != rl) {
        m_restrtlen = rl;
        allocate_scratch();
        m_vv.clear();
    }
}

template <typename V, typename M>
void GMRES<V,M>::define (M& linop)
{
    clear();
    m_linop = &linop;
}

template <typename V, typename M>
void GMRES<V,M>::clear ()
{
    m_its = 0;
    m_status = -1;
    m_res = std::numeric_limits<RT>::max();
    m_v_tmp_rhs.reset();
    m_v_tmp_lhs.reset();
    m_vv.clear();
    m_linop = nullptr;
}

template <typename V, typename M>
bool GMRES<V,M>::converged (RT r0, RT r) const
{
    return (r < r0*m_rtol) || (r < m_atol);
}

template <typename V, typename M>
void GMRES<V,M>::solve (V& a_sol, V const& a_rhs, RT a_tol_rel, RT a_tol_abs, int a_its)
{
    BL_PROFILE("GMRES::solve()");

    AMREX_ALWAYS_ASSERT(m_linop != nullptr);

    auto t0 = amrex::second();

    if (m_v_tmp_rhs == nullptr) {
        m_v_tmp_rhs = std::make_unique<V>(m_linop->makeVecRHS());
    }
    if (m_v_tmp_lhs == nullptr) {
        m_v_tmp_lhs = std::make_unique<V>(m_linop->makeVecLHS());
    }
    if (m_vv.empty()) {
        m_vv.reserve(m_restrtlen+1);
        for (int i = 0; i < 2; ++i) { // to save space, start with just 2
            m_vv.emplace_back(m_linop->makeVecRHS());
        }
    }

    m_rtol = a_tol_rel;
    m_atol = a_tol_abs;

    if (a_its < 0) { a_its = m_maxiter; }

    auto rnorm0 = RT(0);

    m_linop->assign(m_vv[0], a_rhs);
    m_linop->setToZero(a_sol);

    m_its = 0;
    m_status = -1;
    cycle(a_sol, m_status, m_its, rnorm0);

    while (m_status == -1 && m_its < a_its) {
        compute_residual(m_vv[0], a_sol, a_rhs);
        cycle(a_sol, m_status, m_its, rnorm0);
    }

    if (m_status == -1 && m_its >= a_its) { m_status = 1; }

    m_v_tmp_rhs.reset();
    m_v_tmp_lhs.reset();
    m_vv.clear();

    auto t1 = amrex::second();
    if (m_verbose > 0) {
        amrex::Print() << "GMRES: Solve Time = " << t1-t0 << std::endl;
    }
}

template <typename V, typename M>
void GMRES<V,M>::cycle (V& a_xx, int& a_status, int& a_itcount, RT& a_rnorm0)
{
    BL_PROFILE("GMREA::cycle()");

    m_res = m_linop->norm2(m_vv[0]);
    m_grs[0] = m_res;

    if (m_res == RT(0.0)) {
        a_status = 0;
        return;
    }

    m_linop->scale(m_vv[0], RT(1.0)/m_res);

    if (a_itcount == 0) { a_rnorm0 = m_res; }

    a_status = converged(a_rnorm0,m_res) ? 0 : -1;

    int it = 0;
    while (it < m_restrtlen && a_itcount < m_maxiter)
    {
        if (m_verbose > 1) {
            amrex::Print() << "GMRES: iter = " << a_itcount
                           << ", residual = " << m_res << ", " << m_res/a_rnorm0
                           << " (rel.)\n";
        }

        if (a_status == 0) { break; }

        while (m_vv.size() < it+2) {
            m_vv.emplace_back(m_linop->makeVecRHS());
        }

        auto const& vv_it  = m_vv[it  ];
        auto      & vv_it1 = m_vv[it+1];

        m_linop->precond(*m_v_tmp_lhs, vv_it);
        m_linop->apply(vv_it1, *m_v_tmp_lhs);

        gram_schmidt_orthogonalization(it);

        auto tt = m_linop->norm2(vv_it1);

        auto const small = RT((sizeof(RT) == 8) ? 1.e-99 : 1.e-30);
        bool happyend = (tt < small);
        if (!happyend) {
            m_linop->scale(vv_it1, RT(1.0)/tt);
        }

        m_hh (it+1,it) = tt;
        m_hes(it+1,it) = tt;

        update_hessenberg(it, happyend, m_res);

        ++it;
        ++a_itcount;
        a_status = converged(a_rnorm0, m_res) ? 0 : -1;
        if (happyend) { break; }
    }

    if ((m_verbose > 1) && (a_status != 0 || a_itcount >= m_maxiter)) {
        amrex::Print() << "GMRES: iter = " << a_itcount
                       << ", residual = " << m_res << ", " << m_res/a_rnorm0
                       << " (rel.)\n";
    }

    build_solution(a_xx, it-1);
}

template <typename V, typename M>
void GMRES<V,M>::gram_schmidt_orthogonalization (int const it)
{
    // Two unmodified Gram-Schmidt Orthogonalization

    BL_PROFILE("GMRES::GramSchmidt");

    auto& vv_1 = m_vv[it+1];

    Vector<RT> lhh(it+1);

    for (int j = 0; j <= it; ++j) {
        m_hh (j,it) = RT(0.0);
        m_hes(j,it) = RT(0.0);
    }

    for (int ncnt = 0; ncnt < 2 ; ++ncnt)
    {
        for (int j = 0; j <= it; ++j) {
            lhh[j] = m_linop->dotProduct(vv_1, m_vv[j]);
        }

        for (int j = 0; j <= it; ++j) {
            m_linop->increment(vv_1, m_vv[j], -lhh[j]);
            m_hh (j,it) += lhh[j];
            m_hes(j,it) -= lhh[j];
        }
    }
}

template <typename V, typename M>
void GMRES<V,M>::update_hessenberg (int it, bool happyend, RT& res)
{
    BL_PROFILE("GMRES::update_hessenberg()");

    for (int j = 1; j <= it; ++j) {
        auto tt = m_hh(j-1,it);
        m_hh(j-1,it) = m_cc[j-1] * tt         + m_ss[j-1] * m_hh(j,it);
        m_hh(j  ,it) = m_cc[j-1] * m_hh(j,it) - m_ss[j-1] * tt;
    }

    if (!happyend)
    {
        auto tt = std::sqrt(m_hh(it,it)*m_hh(it,it) + m_hh(it+1,it)*m_hh(it+1,it));
        m_cc[it] = m_hh(it  ,it) / tt;
        m_ss[it] = m_hh(it+1,it) / tt;
        m_grs[it+1] = - (m_ss[it] * m_grs[it]);
        m_grs[it  ] =    m_cc[it] * m_grs[it];
        m_hh(it,it) = m_cc[it] * m_hh(it,it) + m_ss[it] * m_hh(it+1,it);
        res = std::abs(m_grs[it+1]);
    }
    else
    {
        res = RT(0.0);
    }
}

template <typename V, typename M>
void GMRES<V,M>::build_solution (V& a_xx, int const it)
{
    BL_PROFILE("GMRES:build_solution()");

    if (it < 0) { return; }

    if (m_hh(it,it) != RT(0.0)) {
        m_grs[it] /= m_hh(it,it);
    } else {
        m_grs[it] = RT(0.0);
    }

    for (int ii = 1; ii <= it; ++ii) {
        int k = it - ii;
        auto tt = m_grs[k];
        for (int j = k+1; j <= it; ++j) {
            tt -= m_hh(k,j) * m_grs[j];
        }
        m_grs[k] = tt / m_hh(k,k);
    }

    m_linop->setToZero(*m_v_tmp_rhs);
    for (int ii = 0; ii < it+1; ++ii) {
        m_linop->increment(*m_v_tmp_rhs, m_vv[ii], m_grs[ii]);
    }

    m_linop->precond(*m_v_tmp_lhs, *m_v_tmp_rhs);
    m_linop->increment(a_xx, *m_v_tmp_lhs, RT(1.0));
}

template <typename V, typename M>
void GMRES<V,M>::compute_residual (V& a_rr, V const& a_xx, V const& a_bb)
{
    BL_PROFILE("GMRES::compute_residual()");
    m_linop->assign(*m_v_tmp_lhs, a_xx);
    m_linop->apply(*m_v_tmp_rhs, *m_v_tmp_lhs);
    m_linop->linComb(a_rr, RT(1.0), a_bb, RT(-1.0), *m_v_tmp_rhs);
}

}
#endif
