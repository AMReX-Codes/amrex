#ifndef AMREX_SMOOTHER_MV_H_
#define AMREX_SMOOTHER_MV_H_

#include <AMReX_Algebra.H>
#include <utility>

namespace amrex {

template <typename T>
class JacobiSmoother
{
public:
    explicit JacobiSmoother (SpMatrix<T> const* a_A) : m_A(a_A) {}

    int setNumIters (int a_niters) { return std::exchange(m_niters, a_niters); }

    void operator() (AlgVector<T>& xvec, AlgVector<T> const& bvec)
    {
        auto const& diag = m_A->diagonalVector();
        AlgVector<T> Axvec(xvec.partition());
        xvec.setVal(0);
        for (int iter = 0; iter < m_niters; ++iter) {
            if (iter == 0) {
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
