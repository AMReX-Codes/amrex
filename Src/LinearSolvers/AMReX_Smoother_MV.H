#ifndef AMREX_SMOOTHER_MV_H_
#define AMREX_SMOOTHER_MV_H_

#include <AMReX_Algebra.H>
#include <AMReX_Array.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_OpenMP.H>
#include <cmath>
#include <memory>
#include <utility>

namespace amrex {

/**
 * \file AMReX_Smoother_MV.H
 *
 * Jacobi, l1-Jacobi, Chebyshev and l1 hybrid Gauss-Seidel smoothers
 * operating on AlgVector/SpMatrix.
 */

/**
 * \brief Jacobi smoother for AlgVector/SpMatrix linear systems.
 *
 * Performs `x += w * (b - A x) / d` sweeps, where `d` is the diagonal of
 * `A` (weighted Jacobi, weight 2/3 by default) or, with `setL1(true)`,
 * the l1 row norm `sign(a_ii) * sum_j |a_ij|` (l1-Jacobi, weight 1 by
 * default), which converges for any symmetric positive definite matrix
 * [Baker, Falgout, Kolev, Yang, SIAM J. Sci. Comput. 33 (2011), 2864-2887].
 */
template <typename T>
class JacobiSmoother
{
public:
    /**
     * \brief Construct a smoother that operates on \p a_A.
     *
     * \param a_A Sparse matrix supplying diagonal entries and SpMV.
     * \param a_l1 Use l1-Jacobi scaling.
     */
    explicit JacobiSmoother (SpMatrix<T> const* a_A, bool a_l1 = false)
        : m_A(a_A), m_l1(a_l1), m_weight(a_l1 ? T(1) : T(2./3.)) {}

    /**
     * \brief Update how many Jacobi sweeps to perform per apply.
     *
     * \param a_niters Desired number of smoothing iterations (4 by default).
     * \return Previous iteration count.
     */
    int setNumIters (int a_niters) { return std::exchange(m_niters, a_niters); }

    //! Set the relaxation weight. Call this after setL1(), which resets it.
    void setWeight (T w) { m_weight = w; }

    //! Switch between l1-Jacobi and weighted Jacobi, resetting the weight
    //! to the default of the new mode.
    void setL1 (bool b) { m_l1 = b; m_weight = b ? T(1) : T(2./3.); m_scaling.reset(); }

    /**
     * \brief Apply the smoother to solve approximately for \p xvec.
     *
     * \param xvec               Solution vector updated in place.
     * \param bvec               RHS vector.
     * \param with_initial_guess True if \p xvec already holds an initial guess.
     */
    void operator() (AlgVector<T>& xvec, AlgVector<T> const& bvec,
                     bool with_initial_guess = false)
    {
        BL_PROFILE("JacobiSmoother");
        AMREX_ASSERT(xvec.partition() == m_A->partition());
        auto const& diag = scaling();
        if ( ! with_initial_guess) {
            xvec.setValAsync(0);
        }
        T const w = m_weight;
        for (int iter = 0; iter < m_niters; ++iter) {
            if ((iter == 0) && ! with_initial_guess) {
                // A x is zero, so this sweep is just x = w b / d.
                ForEach(xvec, bvec, diag,
                        [=] AMREX_GPU_DEVICE (T& x, T const& b, T const& d)
                        {
                            if (d != T(0)) {
                                x += w * b/d;
                            }
                        });
            } else {
                auto& Axvec = workspace(xvec.partition()); // allocated on first need
                SpMV(Axvec, *m_A, xvec);
                ForEach(xvec, Axvec, bvec, diag,
                        [=] AMREX_GPU_DEVICE (T& x, T const& ax, T const& b, T const& d)
                        {
                            if (d != T(0)) {
                                x += w * (b-ax)/d;
                            }
                        });
            }
        }
        Gpu::streamSynchronize();
    }

private:
    SpMatrix<T> const* m_A;
    bool m_l1 = false;
    T m_weight;
    int m_niters = 4;
    // l1 row norms, built on first use. Shared so that the smoother stays
    // copyable (it is used as a std::function preconditioner).
    std::shared_ptr<AlgVector<T>> m_scaling;
    std::shared_ptr<AlgVector<T>> m_work;

    AlgVector<T>& workspace (AlgPartition const& partition)
    {
        if (!m_work || m_work->partition() != partition) {
            m_work = std::make_shared<AlgVector<T>>(partition);
        }
        return *m_work;
    }

public: // NOLINT  Private, but public for CUDA
    //! Diagonal or l1 scaling vector. The first call is collective.
    AlgVector<T> const& scaling ()
    {
        if (!m_l1) { return m_A->diagonalVector(); }
        if (!m_scaling) {
            // The matrix might not have been split yet (e.g., no SpMV so far).
            const_cast<SpMatrix<T>*>(m_A)->setColumnPartition(m_A->partition());
            m_scaling = std::make_shared<AlgVector<T>>(m_A->partition());
            auto* AMREX_RESTRICT p = m_scaling->data();
            auto const& a = m_A->const_parcsr();
            auto const& d = m_A->diagonalVector();
            auto const* AMREX_RESTRICT pd = d.data();
            ParallelForOMP(m_scaling->numLocalRows(), [=] AMREX_GPU_DEVICE (Long i) noexcept
            {
                T s = T(0);
                for (auto idx = a.csr0.row_offset[i]; idx < a.csr0.row_offset[i+1]; ++idx) {
                    s += std::abs(a.csr0.mat[idx]);
                }
                if (a.csr1.nnz > 0 && a.row_map[i] >= 0) {
                    auto ii = a.row_map[i];
                    for (auto idx = a.csr1.row_offset[ii]; idx < a.csr1.row_offset[ii+1]; ++idx) {
                        s += std::abs(a.csr1.mat[idx]);
                    }
                }
                p[i] = std::copysign(s, pd[i]);
            });
            Gpu::streamSynchronize();
        }
        return *m_scaling;
    }
};

/**
 * \brief Chebyshev polynomial smoother for AlgVector/SpMatrix systems.
 *
 * Applies a Chebyshev polynomial of the diagonally (or l1) scaled operator
 * `D^{-1} A` that damps the spectrum in `[lambda_max / ratio, lambda_max]`
 * [Adams, Brezina, Hu, Tuminaro, J. Comput. Phys. 188 (2003), 593-610].
 * Each application costs `degree` matrix-vector products. `lambda_max` is
 * estimated with power iterations on first use and multiplied by a safety
 * factor of 1.1; with l1 scaling it is at most one.
 */
template <typename T>
class ChebyshevSmoother
{
public:
    /**
     * \param a_A Sparse matrix.
     * \param a_l1 Use l1 scaling instead of the diagonal.
     */
    explicit ChebyshevSmoother (SpMatrix<T> const* a_A, bool a_l1 = false)
        : m_jacobi(a_A, a_l1), m_A(a_A), m_l1(a_l1) {}

    //! Polynomial degree (number of matrix-vector products per application, 2 by default).
    void setDegree (int d) { m_degree = d; }
    //! Ratio between the largest eigenvalue and the lower end of the smoothing
    //! interval (20 by default).
    void setEigRatio (T r) { m_ratio = r; }
    //! Number of power iterations for the eigenvalue estimate (10 by default).
    void setNumEigIters (int n) { m_eig_iters = n; }
    //! Number of polynomial applications per call (1 by default).
    int setNumIters (int a_niters) { return std::exchange(m_niters, a_niters); }
    //! Set the largest eigenvalue of the scaled operator instead of estimating it.
    void setLambdaMax (T l) { m_lambda_max = l; }
    //! Largest eigenvalue of the scaled operator, estimated on first use. The
    //! estimate is collective, so all processes must call this together.
    T lambdaMax ()
    {
        if (m_lambda_max <= T(0)) { estimate(); }
        return m_lambda_max;
    }

    /**
     * \brief Apply the smoother to solve approximately for \p xvec.
     *
     * \param xvec               Solution vector updated in place.
     * \param bvec               RHS vector.
     * \param with_initial_guess True if \p xvec already holds an initial guess.
     */
    void operator() (AlgVector<T>& xvec, AlgVector<T> const& bvec,
                     bool with_initial_guess = false)
    {
        BL_PROFILE("ChebyshevSmoother");
        AMREX_ASSERT(xvec.partition() == m_A->partition());
        T const lmax = lambdaMax();
        T const lmin = lmax / m_ratio;
        T const theta = (lmax + lmin) / T(2);
        T const delta = (lmax - lmin) / T(2);
        T const sigma = theta / delta;

        auto const& diag = m_jacobi.scaling();
        if (!m_work || (*m_work)[0].partition() != xvec.partition()) {
            m_work = std::make_shared<Array<AlgVector<T>,3>>();
            for (auto& v : *m_work) { v.define(xvec.partition()); }
        }
        auto& r  = (*m_work)[0];
        auto& d  = (*m_work)[1];
        auto& ad = (*m_work)[2];
        if ( ! with_initial_guess) { xvec.setValAsync(0); }

        T const c0 = T(1) / theta;
        for (int it = 0; it < m_niters; ++it) {
            // r = b - A x; d = D^{-1} r / theta; x += d
            if (it == 0 && ! with_initial_guess) {
                ForEach(xvec, d, r, bvec, diag,
                        [=] AMREX_GPU_DEVICE (T& xi, T& di, T& ri, T const& bi,
                                              T const& dd)
                        {
                            ri = bi;
                            di = (dd != T(0)) ? c0 * ri / dd : T(0);
                            xi += di;
                        });
            } else {
                SpMV(r, *m_A, xvec);
                ForEach(xvec, d, r, bvec, diag,
                        [=] AMREX_GPU_DEVICE (T& xi, T& di, T& ri, T const& bi,
                                              T const& dd)
                        {
                            ri = bi - ri;
                            di = (dd != T(0)) ? c0 * ri / dd : T(0);
                            xi += di;
                        });
            }
            T rho_old = T(1) / sigma;
            for (int k = 1; k < m_degree; ++k) {
                // r -= A d; rho = 1/(2 sigma - rho_old);
                // d = rho rho_old d + (2 rho / delta) D^{-1} r; x += d
                SpMV(ad, *m_A, d);
                T const rho = T(1) / (T(2)*sigma - rho_old);
                T const c1 = rho * rho_old;
                T const c2 = T(2) * rho / delta;
                ForEach(xvec, d, r, ad, diag,
                        [=] AMREX_GPU_DEVICE (T& xi, T& di, T& ri, T const& adi,
                                              T const& dd)
                        {
                            ri -= adi;
                            di = c1 * di + ((dd != T(0)) ? c2 * ri / dd : T(0));
                            xi += di;
                        });
                rho_old = rho;
            }
        }
        Gpu::streamSynchronize();
    }

private:
    JacobiSmoother<T> m_jacobi; // supplies the scaling vector
    SpMatrix<T> const* m_A;
    bool m_l1 = false;
    int m_degree = 2;
    T m_ratio = T(20);
    int m_eig_iters = 10;
    int m_niters = 1;
    T m_lambda_max = T(-1);
    std::shared_ptr<Array<AlgVector<T>,3>> m_work; // r, d, A d

public: // NOLINT  Private, but public for CUDA
    //! Power iteration for the largest eigenvalue of D^{-1} A.
    void estimate ()
    {
        auto const& diag = m_jacobi.scaling();
        AlgVector<T> v(m_A->partition());
        AlgVector<T> w(m_A->partition());
        {
            auto* pv = v.data();
            ParallelForRNG(v.numLocalRows(), [=] AMREX_GPU_DEVICE (Long i, RandomEngine const& eng) noexcept
            {
                pv[i] = T(Random(eng)) - T(0.5);
            });
        }
        T lambda = T(0);
        for (int it = 0; it < m_eig_iters; ++it) {
            T const nv = v.norm2();
            if (nv == T(0)) { break; }
            v.scaleAsync(T(1)/nv);
            SpMV(w, *m_A, v);
            ForEach(w, diag, [=] AMREX_GPU_DEVICE (T& wi, T const& dd)
                    { wi = (dd != T(0)) ? wi / dd : T(0); });
            lambda = amrex::Dot(v, w);
            std::swap(v, w);
        }
        m_lambda_max = lambda * T(1.1);
        if (m_l1) { m_lambda_max = std::min(m_lambda_max, T(1)); }
        if (m_lambda_max <= T(0)) { m_lambda_max = T(1); }
    }
};

/**
 * \brief l1 hybrid Gauss-Seidel smoother for CPUs.
 *
 * Gauss-Seidel within each process, and within each of a fixed number of
 * row blocks that the OpenMP threads share out; couplings to other blocks
 * and to other processes use the previous iterate and their absolute
 * values are added to the diagonal, which keeps the sweep convergent for
 * any SPD matrix and any partition
 * [Baker, Falgout, Kolev, Yang, SIAM J. Sci. Comput. 33 (2011), 2864-2887].
 * Forward and backward sweeps are available to keep a V-cycle symmetric.
 * Not available in GPU builds.
 */
template <typename T>
class L1GaussSeidelSmoother
{
public:
    /**
     * \brief Construct a smoother that operates on \p a_A.
     *
     * \param a_A Sparse matrix supplying the stencil and SpMV.
     */
    explicit L1GaussSeidelSmoother (SpMatrix<T> const* a_A) : m_A(a_A) {}

    //! Number of sweeps per application (2 by default).
    int setNumIters (int a_niters) { return std::exchange(m_niters, a_niters); }

    /**
     * \brief Apply the smoother to solve approximately for \p xvec.
     *
     * \param xvec               Solution vector updated in place.
     * \param bvec               RHS vector.
     * \param with_initial_guess True if \p xvec already holds an initial guess.
     * \param backward           Sweep the rows in descending order.
     */
    void operator() (AlgVector<T>& xvec, AlgVector<T> const& bvec,
                     bool with_initial_guess = false, bool backward = false)
    {
#ifdef AMREX_USE_GPU
        amrex::ignore_unused(xvec, bvec, with_initial_guess, backward);
        amrex::Abort("L1GaussSeidelSmoother is only available in CPU builds");
#else
        BL_PROFILE("L1GaussSeidelSmoother");
        AMREX_ASSERT(xvec.partition() == m_A->partition());
        auto& A = *const_cast<SpMatrix<T>*>(m_A);
        // The matrix might not have been split yet (e.g., no SpMV so far).
        A.setColumnPartition(A.partition());
        auto const& a = A.const_parcsr();
        Long const nrows = xvec.numLocalRows();
        auto const& dl1 = scaling();
        int const nblocks = m_nblocks;
        T const* AMREX_RESTRICT d = dl1.data();
        T const* AMREX_RESTRICT b = bvec.data();
        T* AMREX_RESTRICT x = xvec.data();

        if (!m_work || (*m_work)[0].partition() != xvec.partition()) {
            m_work = std::make_shared<Array<AlgVector<T>,2>>();
            for (auto& v : *m_work) { v.define(xvec.partition()); }
        }
        auto& xold = (*m_work)[0];
        auto& yrem = (*m_work)[1];
        T* AMREX_RESTRICT px_old = xold.data();
        T* AMREX_RESTRICT py_rem = yrem.data();

        if ( ! with_initial_guess) { xvec.setVal(0); }

        for (int iter = 0; iter < m_niters; ++iter) {
            bool const zero_guess = (iter == 0 && ! with_initial_guess);
            if (zero_guess) {
                xold.setVal(0);
                yrem.setVal(0);
            } else {
                xold.copy(xvec);
                // Remote contribution with the previous iterate.
                yrem.setVal(0);
                A.startComm_mv(xvec);
                A.finishComm_mv(yrem);
            }
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
            for (int iblock = 0; iblock < nblocks; ++iblock) {
                auto const [rlo, rhi] = row_block(nrows, iblock, nblocks);
                auto const& c0 = a.csr0;
                for (Long n = 0; n < rhi-rlo; ++n) {
                    Long const i = backward ? (rhi-1-n) : (rlo+n);
                    T s = b[i] - py_rem[i];
                    for (Long idx = c0.row_offset[i]; idx < c0.row_offset[i+1]; ++idx) {
                        Long const j = c0.col_index[idx];
                        if (j == i) { continue; }
                        // Rows of this block use the current iterate, others the old one.
                        s -= c0.mat[idx] * ((j >= rlo && j < rhi) ? x[j] : px_old[j]);
                    }
                    if (d[i] != T(0)) {
                        // (s - a_ii x_i) / d_i with d_i = a_ii + l1 term.
                        x[i] += (s - diag_of(c0, i) * x[i]) / d[i];
                    }
                }
            }
        }
#endif
    }

private:
    SpMatrix<T> const* m_A;
    int m_niters = 2;
    int m_nblocks = -1;
    std::shared_ptr<AlgVector<T>> m_scaling; // a_ii + l1 term for the current row blocks
    std::shared_ptr<Array<AlgVector<T>,2>> m_work; // previous iterate, remote part

    //! Row range of a block. The blocks, not the threads, define the sweep.
    static std::pair<Long,Long> row_block (Long nrows, int iblock, int nblocks)
    {
        Long const lo = (nrows * iblock) / nblocks;
        Long const hi = (nrows * (iblock+1)) / nblocks;
        return {lo, hi};
    }

    template <typename C>
    static T diag_of (C const& c0, Long i)
    {
        for (Long idx = c0.row_offset[i]; idx < c0.row_offset[i+1]; ++idx) {
            if (c0.col_index[idx] == i) { return c0.mat[idx]; }
        }
        return T(0);
    }

#ifndef AMREX_USE_GPU
    AlgVector<T> const& scaling ()
    {
        int const nb = OpenMP::get_max_threads();
        if (!m_scaling || nb != m_nblocks) {
            m_nblocks = nb;
            m_scaling = std::make_shared<AlgVector<T>>(m_A->partition());
            T* AMREX_RESTRICT p = m_scaling->data();
            auto const& a = m_A->const_parcsr();
            Long const nrows = m_scaling->numLocalRows();
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
            for (int iblock = 0; iblock < nb; ++iblock) {
                auto const [rlo, rhi] = row_block(nrows, iblock, nb);
                auto const& c0 = a.csr0;
                for (Long i = rlo; i < rhi; ++i) {
                    T const aii = diag_of(c0, i); // same entry the sweep uses
                    T l1 = T(0);
                    for (Long idx = c0.row_offset[i]; idx < c0.row_offset[i+1]; ++idx) {
                        Long const j = c0.col_index[idx];
                        if (j != i && (j < rlo || j >= rhi)) {
                            l1 += std::abs(c0.mat[idx]);
                        }
                    }
                    if (a.csr1.nnz > 0 && a.row_map[i] >= 0) {
                        Long const ii = a.row_map[i];
                        for (Long idx = a.csr1.row_offset[ii]; idx < a.csr1.row_offset[ii+1]; ++idx) {
                            l1 += std::abs(a.csr1.mat[idx]);
                        }
                    }
                    p[i] = aii + std::copysign(l1, aii);
                }
            }
        }
        return *m_scaling;
    }
#endif
};

}

#endif
