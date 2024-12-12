#ifndef AMREX_GMRES_MV_H_
#define AMREX_GMRES_MV_H_

#include <AMReX_GMRES.H>
#include <AMReX_Algebra.H>
#include <AMReX_Smoother_MV.H>
#include <functional>

namespace amrex {

template <typename T>
class GMRES_MV
{
public:
    using RT = Real; // double or float
    using VEC = AlgVector<T>;
    using MAT = SpMatrix<T>;
    using GM = amrex::GMRES<VEC,GMRES_MV<T>>;
    using PC = std::function<void(VEC&,VEC const&)>;

    GMRES_MV (MAT const* a_mat);

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

    //! Sets verbosity.
    void setVerbose (int v) { m_gmres.setVerbose(v); }

    //! Get the GMRES object.
    GM& getGMRES () { return m_gmres; }

    //! Make MultiFab without ghost cells
    [[nodiscard]] VEC makeVecRHS () const;

    //! Make MultiFab with ghost cells and set ghost cells to zero
    [[nodiscard]] VEC makeVecLHS () const;

    static T norm2 (VEC const& vec);

    static void scale (VEC& vec, T scale_factor);

    static T dotProduct (VEC const& vec1, VEC const& vec2);

    //! lhs = 0
    static void setToZero (VEC& lhs);

    //! lhs = rhs
    static void assign (VEC& lhs, VEC const& rhs);

    //! lhs += a*rhs
    static void increment (VEC& lhs, VEC const& rhs, T a);

    //! lhs = a*rhs_a + b*rhs_b
    static void linComb (VEC& lhs, T a, VEC const& rhs_a, T b, VEC const& rhs_b);

    //! lhs = L(rhs)
    void apply (VEC& lhs, VEC& rhs) const;

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
    amrex::Axpy(lhs, a, rhs, true);
}

template <typename T>
void GMRES_MV<T>::linComb (VEC& lhs, T a, VEC const& rhs_a, T b, VEC const& rhs_b)
{
    amrex::LinComb(lhs, a, rhs_a, b, rhs_b, true);
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
