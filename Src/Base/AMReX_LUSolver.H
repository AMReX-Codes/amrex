#ifndef AMREX_LU_SOLVER_H_
#define AMREX_LU_SOLVER_H_
#include <AMReX_Config.H>

#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <cmath>
#include <limits>

namespace amrex {

// https://en.wikipedia.org/wiki/LU_decomposition

template <int N, typename T>
class LUSolver
{
public:

    LUSolver () = default;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    LUSolver (Array2D<T, 0, N-1, 0, N-1, Order::C> const& a_mat);

    void define (Array2D<T, 0, N-1, 0, N-1, Order::C> const& a_mat);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (T* AMREX_RESTRICT x, T const* AMREX_RESTRICT b) const
    {
        for (int i = 0; i < N; ++i) {
            x[i] = b[m_piv(i)];
            for (int k = 0; k < i; ++k) {
                x[i] -= m_mat(i,k) * x[k];
            }
        }

        for (int i = N-1; i >= 0; --i) {
            for (int k = i+1; k < N; ++k) {
                x[i] -= m_mat(i,k) * x[k];
            }
            x[i] *= m_mat(i,i);
        }
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Array2D<T,0,N-1,0,N-1,Order::C> invert () const
    {
        Array2D<T,0,N-1,0,N-1,Order::C> IA;
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                IA(i,j) = (m_piv(i) == j) ? T(1.0) : T(0.0);
                for (int k = 0; k < i; ++k) {
                    IA(i,j) -= m_mat(i,k) * IA(k,j);
                }
            }
            for (int i = N-1; i >= 0; --i) {
                for (int k = i+1; k < N; ++k) {
                    IA(i,j) -= m_mat(i,k) * IA(k,j);
                }
                IA(i,j) *= m_mat(i,i);
            }
        }
        return IA;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    T determinant () const
    {
        T det = m_mat(0,0);
        for (int i = 1; i < N; ++i) {
            det *= m_mat(i,i);
        }
        det = T(1.0) / det;
        return (m_npivs % 2 == 0) ? det : -det;
    }

private:

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void define_innard ();

    Array2D<T, 0, N-1, 0, N-1, Order::C> m_mat;
    Array1D<int, 0, N-1> m_piv;
    int m_npivs = 0;
};

template <int N, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
LUSolver<N,T>::LUSolver (Array2D<T, 0, N-1, 0, N-1, Order::C> const& a_mat)
    : m_mat(a_mat)
{
    define_innard();
}

template <int N, typename T>
void LUSolver<N,T>::define (Array2D<T, 0, N-1, 0, N-1, Order::C> const& a_mat)
{
    m_mat = a_mat;
    define_innard();
}

template <int N, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void LUSolver<N,T>::define_innard ()
{
    static_assert(N > 1);
    static_assert(std::is_floating_point_v<T>);

    for (int i = 0; i < N; ++i) { m_piv(i) = i; }
    m_npivs = 0;

    for (int i = 0; i < N; ++i) {
        T maxA = 0;
        int imax = i;

        for (int k = i; k < N; ++k) {
            auto const absA = std::abs(m_mat(k,i));
            if (absA > maxA) {
                maxA = absA;
                imax = k;
            }
        }

        if (maxA < std::numeric_limits<T>::min()) {
            amrex::Abort("LUSolver: matrix is degenerate");
        }

        if (imax != i) {
            amrex::Swap(m_piv(i), m_piv(imax));
            for (int j = 0; j < N; ++j) {
                amrex::Swap(m_mat(i,j), m_mat(imax,j));
            }
            ++m_npivs;
        }

        for (int j = i+1; j < N; ++j) {
            m_mat(j,i) /= m_mat(i,i);
            for (int k = i+1; k < N; ++k) {
                m_mat(j,k) -= m_mat(j,i) * m_mat(i,k);
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        m_mat(i,i) = T(1) / m_mat(i,i);
    }
}

}

#endif
