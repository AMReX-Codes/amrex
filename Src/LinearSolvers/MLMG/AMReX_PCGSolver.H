#ifndef AMREX_PCG_SOLVER_H_
#define AMREX_PCG_SOLVER_H_
#include <AMReX_Config.H>

#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <cmath>
#include <type_traits>

namespace amrex {

/**
 * \brief Preconditioned conjugate gradient solver
 *
 * \param x       initial guess
 * \param r       initial residual
 * \param mat     matrix
 * \param precond preconditioner
 * \param maxiter max number of iterations
 * \param rel_tol relative tolerance
 */
template <int N, typename T, typename M, typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int pcg_solve (T* AMREX_RESTRICT x, T* AMREX_RESTRICT r,
               M const& mat, P const& precond, int maxiter, T rel_tol)
{
    static_assert(std::is_floating_point_v<T>);

    T rnorm0 = 0;
    for (int i = 0; i < N; ++i) {
        rnorm0 = std::max(rnorm0, std::abs(r[i]));
    }
    if (rnorm0 == 0) { return 0; }

    int iter = 0;
    T rho_prev = T(1.0); // initialized to quiet gcc warning
    T p[N];
    for (iter = 1; iter <= maxiter; ++iter) {
        T z[N];
        precond(z, r);
        T rho = 0;
        for (int i = 0; i < N; ++i) { rho += r[i]*z[i]; }
        if (rho == 0) { break; }
        if (iter == 1) {
            for (int i = 0; i < N; ++i) { p[i] = z[i]; }
        } else {
            auto rr = rho * (T(1.0)/rho_prev);
            for (int i = 0; i < N; ++i) {
                p[i] = z[i] + rr * p[i];
            }
        }
        T q[N];
        mat(q, p);
        T pq = 0;
        for (int i = 0; i < N; ++i) { pq += p[i]*q[i]; }
        if (pq == 0) { break; }
        T alpha = rho * (T(1.0)/pq);
        T rnorm = 0;
        for (int i = 0; i < N; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
            rnorm = std::max(rnorm, std::abs(r[i]));
        }
        if (rnorm <= rnorm0*rel_tol) { break; }
        rho_prev = rho;
    }
    return iter;
}

}

#endif
