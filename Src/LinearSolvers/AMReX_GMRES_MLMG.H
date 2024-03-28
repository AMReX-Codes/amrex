#ifndef AMREX_GMRES_MLMG_H_
#define AMREX_GMRES_MLMG_H_
#include <AMReX_Config.H>

#include <AMReX_GMRES.H>
#include <AMReX_MLMG.H>
#include <utility>

namespace amrex {

/**
 * \brief Solve using GMRES with multigrid as preconditioner
 *
 * The linear system to solve is provided by MLMG, which is also being used
 * as the preconditioner.
 *
 */
template <typename MF>
class GMRESMLMGT
{
public:
    using MG = MLMGT<MF>;
    using RT = typename MG::RT; // double or float
    using GM = GMRES<MF,GMRESMLMGT<MF>>;

    explicit GMRESMLMGT (MG& mlmg);

    /**
     * \brief Solve the linear system
     *
     * \param a_sol     unknowns, i.e., x in A x = b.
     * \param a_rhs     RHS, i.e., b in A x = b.
     * \param a_tol_rel relative tolerance.
     * \param a_tol_abs absolute tolerance.
     */
    void solve (MF& a_sol, MF const& a_rhs, RT a_tol_rel, RT a_tol_abs);

    //! Sets verbosity.
    void setVerbose (int v) { m_gmres.setVerbose(v); }

    //! Sets the max number of iterations
    void setMaxIters (int niters) { m_gmres.setMaxIters(niters); }

    //! Gets the number of iterations.
    [[nodiscard]] int getNumIters () const { return m_gmres.getNumIters(); }

    //! Gets the 2-norm of the residual.
    [[nodiscard]] RT getResidualNorm () const { return m_gmres.getResidualNorm(); }

    //! Get the GMRES object.
    GM& getGMRES () { return m_gmres; }

    /**
     * \brief Set MLMG's multiplicative property of zero
     *
     * This should NOT be called unless MLMG has the multiplicative property
     * of zero. MLMG is not a matrix, and usually does not have the
     * properties of a matrix such as the multiplicative property of zero
     * (i.e., M*0=0) because how domain boundary conditions are
     * handled. However, if MLMG has the property of zero, calling this
     * function with true can have a small performance benefit.
     */
    void setPropertyOfZero (bool b) { m_prop_zero = b; }

    //! Make MultiFab without ghost cells
    MF makeVecRHS () const;

    //! Make MultiFab with ghost cells and set ghost cells to zero
    MF makeVecLHS () const;

    RT norm2 (MF const& mf) const;

    static void scale (MF& mf, RT scale_factor);

    RT dotProduct (MF const& mf1, MF const& mf2) const;

    //! lhs = 0
    static void setToZero (MF& lhs);

    //! lhs = rhs
    static void assign (MF& lhs, MF const& rhs);

    //! lhs += a*rhs
    static void increment (MF& lhs, MF const& rhs, RT a);

    //! lhs = a*rhs_a + b*rhs_b
    static void linComb (MF& lhs, RT a, MF const& rhs_a, RT b, MF const& rhs_b);

    //! lhs = L(rhs)
    void apply (MF& lhs, MF const& rhs) const;

    void precond (MF& lhs, MF const& rhs) const;

    bool usePrecond (bool new_flag) { return std::exchange(m_use_precond, new_flag); }

private:
    GM m_gmres;
    MG* m_mlmg;
    MLLinOpT<MF>* m_linop;
    bool m_use_precond = true;
    bool m_prop_zero = false;
};

template <typename MF>
GMRESMLMGT<MF>::GMRESMLMGT (MG& mlmg)
    : m_mlmg(&mlmg), m_linop(&mlmg.getLinOp())
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_linop->NAMRLevels() == 1,
                                     "Only support single level solve");
    m_mlmg->setVerbose(0);
    m_mlmg->setBottomVerbose(0);
    m_mlmg->prepareForGMRES();
    m_gmres.define(*this);
}

template <typename MF>
auto GMRESMLMGT<MF>::makeVecRHS () const -> MF
{
    return m_linop->make(0, 0, IntVect(0));
}

template <typename MF>
auto GMRESMLMGT<MF>::makeVecLHS () const -> MF
{
    auto mf = m_linop->make(0, 0, IntVect(1));
    setBndry(mf, RT(0), 0, nComp(mf));
    return mf;
}

template <typename MF>
auto GMRESMLMGT<MF>::norm2 (MF const& mf) const -> RT
{
    auto r = m_linop->xdoty(0, 0, mf, mf, false);
    return std::sqrt(r);
}

template <typename MF>
void GMRESMLMGT<MF>::scale (MF& mf, RT scale_factor)
{
    Scale(mf, scale_factor, 0, nComp(mf), 0);
}

template <typename MF>
auto GMRESMLMGT<MF>::dotProduct (MF const& mf1, MF const& mf2) const -> RT
{
    return m_linop->xdoty(0, 0, mf1, mf2, false);
}

template <typename MF>
void GMRESMLMGT<MF>::setToZero (MF& lhs)
{
    setVal(lhs, RT(0.0));
}

template <typename MF>
void GMRESMLMGT<MF>::assign (MF& lhs, MF const& rhs)
{
    LocalCopy(lhs, rhs, 0, 0, nComp(lhs), IntVect(0));
}

template <typename MF>
void GMRESMLMGT<MF>::increment (MF& lhs, MF const& rhs, RT a)
{
    Saxpy(lhs, a, rhs, 0, 0, nComp(lhs), IntVect(0));
}

template <typename MF>
void GMRESMLMGT<MF>::linComb (MF& lhs, RT a, MF const& rhs_a, RT b, MF const& rhs_b)
{
    LinComb(lhs, a, rhs_a, 0, b, rhs_b, 0, 0, nComp(lhs), IntVect(0));
}

template <typename MF>
void GMRESMLMGT<MF>::apply (MF& lhs, MF const& rhs) const
{
    m_linop->apply(0, 0, lhs, const_cast<MF&>(rhs),
                   MLLinOpT<MF>::BCMode::Homogeneous,
                   MLLinOpT<MF>::StateMode::Correction);
}

template <typename MF>
void GMRESMLMGT<MF>::precond (MF& lhs, MF const& rhs) const
{
    if (m_use_precond) {
        m_mlmg->prepareMGcycle();

        LocalCopy(m_mlmg->res[0][0], rhs, 0, 0, nComp(rhs), IntVect(0));

        m_mlmg->mgVcycle(0,0);

        LocalCopy(lhs, m_mlmg->cor[0][0], 0, 0, nComp(rhs), IntVect(0));

    } else {
        LocalCopy(lhs, rhs, 0, 0, nComp(lhs), IntVect(0));
    }
}

template <typename MF>
void GMRESMLMGT<MF>::solve (MF& a_sol, MF const& a_rhs, RT a_tol_rel, RT a_tol_abs)
{
    if (m_prop_zero) {
        auto rhs = makeVecRHS();
        assign(rhs, a_rhs);
        m_linop->setDirichletNodesToZero(0,0,rhs);
        m_gmres.solve(a_sol, rhs, a_tol_rel, a_tol_abs);
    } else {
        auto res = makeVecRHS();
        m_mlmg->apply({&res}, {&a_sol}); // res = L(sol)
        increment(res, a_rhs, RT(-1)); // res = L(sol) - rhs
        auto cor = makeVecLHS();
        m_linop->setDirichletNodesToZero(0,0,res);
        m_gmres.solve(cor, res, a_tol_rel, a_tol_abs); // L(cor) = res
        increment(a_sol, cor, RT(-1)); // sol = sol - cor
    }
}

using GMRESMLMG = GMRESMLMGT<MultiFab>;

}

#endif
