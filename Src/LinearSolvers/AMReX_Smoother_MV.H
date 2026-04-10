#ifndef AMREX_SMOOTHER_MV_H_
#define AMREX_SMOOTHER_MV_H_

#include <AMReX_Algebra.H>
#include <utility>

namespace amrex {

/**
 * \file AMReX_Smoother_MV.H
 *
 * Simple weighted-Jacobi smoother operating on AlgVector/SpMatrix.
 */

/**
 * \brief Weighted-Jacobi smoother for AlgVector/SpMatrix linear systems.
 *
 * Encapsulates a configurable number of damped Jacobi sweeps that operate
 * on the diagonal extracted from a supplied sparse matrix.
 */
template <typename T>
class JacobiSmoother
{
public:
    /**
     * \brief Construct a weighted-Jacobi smoother that operates on \p a_A.
     *
     * \param a_A Sparse matrix supplying diagonal entries and SpMV.
     */
    explicit JacobiSmoother (SpMatrix<T> const* a_A) : m_A(a_A) {}

    /**
     * \brief Update how many Jacobi sweeps to perform per apply.
     *
     * \param a_niters Desired number of smoothing iterations.
     * \return Previous iteration count.
     */
    int setNumIters (int a_niters) { return std::exchange(m_niters, a_niters); }

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
        auto const& diag = m_A->diagonalVector();
        AlgVector<T> Axvec(xvec.partition());
        if ( ! with_initial_guess) {
            xvec.setVal(0);
        }
        for (int iter = 0; iter < m_niters; ++iter) {
            if ((iter == 0) && ! with_initial_guess) {
                Axvec.setVal(0);
            } else {
                SpMV(Axvec, *m_A, xvec);
            }
            ForEach(xvec, Axvec, bvec, diag,
                    [=] AMREX_GPU_DEVICE (T& x, T const& ax, T const& b, T const& d)
                    {
                        if (d != T(0)) {
                            x += (b-ax)/d * T(2./3.); // weighted Jacobi
                        }
                    });
        }
        Gpu::streamSynchronize();
    }

private:
    SpMatrix<T> const* m_A;
    int m_niters = 4;
};

}

#endif
