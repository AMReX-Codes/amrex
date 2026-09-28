#ifndef AMREX_KRYLOV_MV_H_
#define AMREX_KRYLOV_MV_H_

#include <AMReX_Algebra.H>
#include <AMReX_Smoother_MV.H>
#include <functional>
#include <type_traits>

namespace amrex {

/**
 * \file AMReX_KrylovMV.H
 *
 * Krylov solver front-end specialized for AlgVector/SpMatrix combinations.
 */

template <typename V, typename M> class GMRES;

/**
 * \brief Binds a Krylov solver to an SpMatrix.
 *
 * \tparam S Krylov solver class template with the GMRES operator interface,
 *           e.g., GMRES, BiCGStab or PCG.
 * \tparam T Value type of the matrix and vectors (double or float).
 *
 * Use the aliases GMRES_MV, BiCGStab_MV and PCG_MV.
 */
template <template <typename, typename> class S, typename T>
class KrylovMV
{
public:
    using RT = T;
    using VEC = AlgVector<T>;
    using MAT = SpMatrix<T>;
    using Solver = S<VEC,KrylovMV<S,T>>;
    using GM = Solver;
    using PC = std::function<void(VEC&,VEC const&)>;

    /**
     * \brief Bind the solver to a sparse matrix described by \p a_mat.
     *
     * \param a_mat Matrix supplying the SpMV operation and partition.
     */
    KrylovMV (MAT const* a_mat);

    ~KrylovMV () = default;

    // The solver holds a pointer to this object.
    KrylovMV (const KrylovMV&) = delete;
    KrylovMV (KrylovMV&&) = delete;
    KrylovMV& operator= (const KrylovMV&) = delete;
    KrylovMV& operator= (KrylovMV&&) = delete;

    /**
     * \brief Supply an optional right-preconditioner functor.
     *
     * The functor must follow the signature `void(VEC& lhs, VEC const& rhs)`.
     *
     * \param a_pc Callable that computes `lhs = P^{-1}(rhs)`.
     */
    void setPrecond (PC a_pc) { m_pc = std::move(a_pc); }

    /**
     * \brief Solve the linear system
     *
     * \param a_sol     unknowns, i.e., x in A x = b.
     * \param a_rhs     RHS, i.e., b in A x = b.
     * \param a_tol_rel relative tolerance.
     * \param a_tol_abs absolute tolerance.
     */
    void solve (VEC& a_sol, VEC const& a_rhs, T a_tol_rel, T a_tol_abs);

    //! Set verbosity level \p v for the underlying solver.
    void setVerbose (int v) { m_solver.setVerbose(v); }

    //! Access the underlying solver for additional settings and results.
    Solver& getSolver () { return m_solver; }

    //! Same as getSolver(), for GMRES_MV only.
    Solver& getGMRES () requires std::is_same_v<Solver, GMRES<VEC,KrylovMV<S,T>>>
    { return m_solver; }

    //! Return an AlgVector that mirrors the matrix partition for RHS storage.
    [[nodiscard]] VEC makeVecRHS () const;

    //! Return another AlgVector with the same partition for LHS storage.
    [[nodiscard]] VEC makeVecLHS () const;

    //! Euclidean norm of \p vec.
    static T norm2 (VEC const& vec);

    //! Scale \p vec in place by \p scale_factor.
    static void scale (VEC& vec, T scale_factor);

    //! Dot product between \p vec1 and \p vec2.
    static T dotProduct (VEC const& vec1, VEC const& vec2);

    //! Reset \p lhs to zero.
    static void setToZero (VEC& lhs);

    //! Copy \p rhs into \p lhs.
    static void assign (VEC& lhs, VEC const& rhs);

    //! Accumulate `lhs += a * rhs`.
    static void increment (VEC& lhs, VEC const& rhs, T a);

    //! Form the linear combination `lhs = a * rhs_a + b * rhs_b`.
    static void linComb (VEC& lhs, T a, VEC const& rhs_a, T b, VEC const& rhs_b);

    /**
     * \brief Apply the sparse operator: `lhs = A rhs`.
     *
     * \p rhs is non-const to satisfy the solver interface.
     */
    void apply (VEC& lhs, VEC& rhs) const;

    //! Apply the optional preconditioner (or copy if none is provided).
    void precond (VEC& lhs, VEC const& rhs) const;

private:
    Solver m_solver;
    MAT const* m_mat = nullptr;
    PC m_pc;
};

template <template <typename, typename> class S, typename T>
KrylovMV<S,T>::KrylovMV (MAT const* a_mat)
    : m_mat(a_mat)
{
    m_solver.define(*this);
}

template <template <typename, typename> class S, typename T>
void KrylovMV<S,T>::solve (VEC& a_sol, VEC const& a_rhs, T a_tol_rel, T a_tol_abs)
{
    m_solver.solve(a_sol, a_rhs, a_tol_rel, a_tol_abs);
}

template <template <typename, typename> class S, typename T>
auto KrylovMV<S,T>::makeVecRHS () const -> VEC
{
    return VEC(m_mat->partition());
}

template <template <typename, typename> class S, typename T>
auto KrylovMV<S,T>::makeVecLHS () const -> VEC
{
    return VEC(m_mat->partition());
}

template <template <typename, typename> class S, typename T>
T KrylovMV<S,T>::norm2 (VEC const& vec)
{
    return vec.norm2();
}

template <template <typename, typename> class S, typename T>
void KrylovMV<S,T>::scale (VEC& vec, T scale_factor)
{
    vec.scaleAsync(scale_factor);
}

template <template <typename, typename> class S, typename T>
T KrylovMV<S,T>::dotProduct (VEC const& vec1, VEC const& vec2)
{
    return amrex::Dot(vec1,vec2);
}

template <template <typename, typename> class S, typename T>
void KrylovMV<S,T>::setToZero (VEC& lhs)
{
    lhs.setValAsync(0);
}

template <template <typename, typename> class S, typename T>
void KrylovMV<S,T>::assign (VEC& lhs, VEC const& rhs)
{
    lhs.copyAsync(rhs);
}

template <template <typename, typename> class S, typename T>
void KrylovMV<S,T>::increment (VEC& lhs, VEC const& rhs, T a)
{
    amrex::Axpy(lhs, a, rhs);
}

template <template <typename, typename> class S, typename T>
void KrylovMV<S,T>::linComb (VEC& lhs, T a, VEC const& rhs_a, T b, VEC const& rhs_b)
{
    amrex::LinComb(lhs, a, rhs_a, b, rhs_b);
}

template <template <typename, typename> class S, typename T>
void KrylovMV<S,T>::apply (VEC& lhs, VEC& rhs) const
{
    amrex::SpMV(lhs, *m_mat, rhs);
}

template <template <typename, typename> class S, typename T>
void KrylovMV<S,T>::precond (VEC& lhs, VEC const& rhs) const
{
    if (m_pc) {
        m_pc(lhs, rhs);
    } else {
        lhs.copyAsync(rhs);
    }
}

}
#endif
