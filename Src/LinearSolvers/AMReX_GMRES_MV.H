#ifndef AMREX_GMRES_MV_H_
#define AMREX_GMRES_MV_H_

#include <AMReX_GMRES.H>
#include <AMReX_Algebra.H>
#include <AMReX_Smoother_MV.H>
#include <functional>

namespace amrex {

/**
 * \file AMReX_GMRES_MV.H
 *
 * GMRES front-end specialized for AlgVector/SpMatrix combinations.
 */

template <typename T>
class GMRES_MV
{
public:
    using RT = Real; // double or float
    using VEC = AlgVector<T>;
    using MAT = SpMatrix<T>;
    using GM = amrex::GMRES<VEC,GMRES_MV<T>>;
    using PC = std::function<void(VEC&,VEC const&)>;

    /**
     * \brief Bind GMRES to a sparse matrix described by \p a_mat.
     *
     * \param a_mat Matrix supplying the SpMV operation and partition.
     */
    GMRES_MV (MAT const* a_mat);

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

    //! Set verbosity level \p v for the underlying GMRES solver.
    void setVerbose (int v) { m_gmres.setVerbose(v); }

    //! Access the underlying GMRES object for additional tuning.
    GM& getGMRES () { return m_gmres; }

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

    /**
     * \brief Reset \p lhs to zero.
     *
     * \param lhs Vector whose entries are zeroed asynchronously.
     */
    static void setToZero (VEC& lhs);

    /**
     * \brief Copy \p rhs into \p lhs.
     *
     * \param lhs Destination vector updated in place.
     * \param rhs Source vector.
     */
    static void assign (VEC& lhs, VEC const& rhs);

    /**
     * \brief Accumulate \f$\text{lhs} \mathrel{+}= a\,\text{rhs}\f$.
     *
     * \param lhs Accumulation target.
     * \param rhs Source vector.
     * \param a   Scalar factor applied to \p rhs.
     */
    static void increment (VEC& lhs, VEC const& rhs, T a);

    /**
     * \brief Form the linear combination \f$\text{lhs} = a\,\text{rhs}_a + b\,\text{rhs}_b\f$.
     *
     * \param lhs    Destination vector.
     * \param a      Scalar applied to \p rhs_a.
     * \param rhs_a  First source vector.
     * \param b      Scalar applied to \p rhs_b.
     * \param rhs_b  Second source vector.
     */
    static void linComb (VEC& lhs, T a, VEC const& rhs_a, T b, VEC const& rhs_b);

    /**
     * \brief Apply the sparse operator to \p rhs and store the result in \p lhs.
     *
     * \param lhs Destination vector updated in place.
     * \param rhs Source vector produced by the preconditioner (passed by non-const reference to satisfy the GMRES interface).
     */
    void apply (VEC& lhs, VEC& rhs) const;

    //! Apply the optional preconditioner (or copy if none is provided).
    void precond (VEC& lhs, VEC const& rhs) const;

private:
    GM m_gmres;
    MAT const* m_mat = nullptr;
    PC m_pc;
};

template <typename T>
GMRES_MV<T>::GMRES_MV (MAT const* a_mat)
    : m_mat(a_mat)
{
    m_gmres.define(*this);
}

template <typename T>
void GMRES_MV<T>::solve (VEC& a_sol, VEC const& a_rhs, T a_tol_rel, T a_tol_abs)
{
    m_gmres.solve(a_sol, a_rhs, a_tol_rel, a_tol_abs);
}

template <typename T>
auto GMRES_MV<T>::makeVecRHS () const -> VEC
{
    return VEC(m_mat->partition());
}

template <typename T>
auto GMRES_MV<T>::makeVecLHS () const -> VEC
{
    return VEC(m_mat->partition());
}

template <typename T>
T GMRES_MV<T>::norm2 (VEC const& vec)
{
    return vec.norm2();
}

template <typename T>
void GMRES_MV<T>::scale (VEC& vec, T scale_factor)
{
    vec.scaleAsync(scale_factor);
}

template <typename T>
T GMRES_MV<T>::dotProduct (VEC const& vec1, VEC const& vec2)
{
    return amrex::Dot(vec1,vec2);
}

template <typename T>
void GMRES_MV<T>::setToZero (VEC& lhs)
{
    lhs.setValAsync(0);
}

template <typename T>
void GMRES_MV<T>::assign (VEC& lhs, VEC const& rhs)
{
    lhs.copyAsync(rhs);
}

template <typename T>
void GMRES_MV<T>::increment (VEC& lhs, VEC const& rhs, T a)
{
    amrex::Axpy(lhs, a, rhs);
}

template <typename T>
void GMRES_MV<T>::linComb (VEC& lhs, T a, VEC const& rhs_a, T b, VEC const& rhs_b)
{
    amrex::LinComb(lhs, a, rhs_a, b, rhs_b);
}

template <typename T>
void GMRES_MV<T>::apply (VEC& lhs, VEC& rhs) const
{
    amrex::SpMV(lhs, *m_mat, rhs);
}

template <typename T>
void GMRES_MV<T>::precond (VEC& lhs, VEC const& rhs) const
{
    if (m_pc) {
        m_pc(lhs, rhs);
    } else {
        lhs.copyAsync(rhs);
    }
}

}
#endif
