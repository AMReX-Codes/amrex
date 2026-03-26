#ifndef AMREX_GMRES_MLMG_H_
#define AMREX_GMRES_MLMG_H_
#include <AMReX_Config.H>

#include <AMReX_GMRES.H>
#include <AMReX_MLMG.H>
#include <AMReX_TypeTraits.H>
#include <utility>

namespace amrex {

/**
 * \file AMReX_GMRES_MLMG.H
 *
 * GMRES wrapper that plugs into MLMG operators/preconditioners.
 */

/**
 * \ingroup amrex_solver_gmres
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
    using VEC = Vector<MF>;
    using MG = MLMGT<MF>;
    using RT = typename MG::RT; // double or float
    using GM = GMRES<Vector<MF>,GMRESMLMGT<MF>>;

    /**
     * \brief Wrap an existing MLMG hierarchy as a GMRES operator/preconditioner.
     *
     * \param mlmg Multigrid object that provides apply/precond hooks.
     */
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

    /**
     * \brief Solve using pointer-to-MultiFab arrays (utility for existing MLMG data).
     *
     * \param a_sol     Per-level solution MultiFabs (single component).
     * \param a_rhs     Per-level RHS MultiFabs (single component).
     * \param a_tol_rel Relative tolerance.
     * \param a_tol_abs Absolute tolerance.
     */
    void solve (Vector<MF*> const& a_sol, Vector<MF const*> const& a_rhs, RT a_tol_rel, RT a_tol_abs);

    //! Set verbosity level \p v for both GMRES and MLMG components.
    void setVerbose (int v) {
        m_verbose = v;
        m_gmres.setVerbose(v);
    }

    //! Cap the number of GMRES iterations executed per solve to \p niters.
    void setMaxIters (int niters) { m_gmres.setMaxIters(niters); }

    //! Number of iterations executed by the last solve().
    [[nodiscard]] int getNumIters () const { return m_gmres.getNumIters(); }

    //! Final residual 2-norm produced by the last solve().
    [[nodiscard]] RT getResidualNorm () const { return m_gmres.getResidualNorm(); }

    //! Direct access to the underlying GMRES driver for fine-grained tuning.
    GM& getGMRES () { return m_gmres; }

    //! Return a MultiFab vector without ghost cells for RHS storage.
    VEC makeVecRHS () const;

    //! Return a MultiFab vector with one layer of ghost cells for LHS storage.
    //! Also set the data in the ghost cells to zero.
    VEC makeVecLHS () const;

    //! Return the 2-norm of \p mf using the multigrid operator's weighting.
    RT norm2 (VEC const& mf) const;

    //! Scale each component of \p mf by \p scale_factor.
    static void scale (VEC& mf, RT scale_factor);

    //! Return the preconditioner-aware dot product between \p mf1 and \p mf2.
    RT dotProduct (VEC const& mf1, VEC const& mf2) const;

    /**
     * \brief Reset \p lhs to zero on every AMR level.
     *
     * \param lhs Per-level vector that will be zeroed in place.
     */
    static void setToZero (VEC& lhs);

    /**
     * \brief Copy \p rhs into \p lhs level by level.
     *
     * \param lhs Destination vector updated in place.
     * \param rhs Source vector.
     */
    static void assign (VEC& lhs, VEC const& rhs);

    /**
     * \brief Add a scaled vector: \f$ \text{lhs} \mathrel{+}= a\,\text{rhs}\f$.
     *
     * \param lhs Accumulation target.
     * \param rhs Source vector.
     * \param a   Scalar multiplier applied to \p rhs.
     */
    static void increment (VEC& lhs, VEC const& rhs, RT a);

    /**
     * \brief Form a linear combination: \f$ \text{lhs} = a\,\text{rhs}_a + b\,\text{rhs}_b\f$.
     *
     * \param lhs    Destination vector written in place.
     * \param a      Scalar applied to \p rhs_a.
     * \param rhs_a  First source vector.
     * \param b      Scalar applied to \p rhs_b.
     * \param rhs_b  Second source vector.
     */
    static void linComb (VEC& lhs, RT a, VEC const& rhs_a, RT b, VEC const& rhs_b);

    /**
     * \brief Apply the operator supplied to the constructor.
     *
     * \param lhs Destination vector receiving \f$L(\text{rhs})\f$.
     * \param rhs Source vector.
     */
    void apply (VEC& lhs, VEC const& rhs) const;

    /**
     * \brief Apply the MLMG-based preconditioner: \f$\text{lhs} = P^{-1}(\text{rhs})\f$.
     *
     * \param lhs Destination vector, zeroed before solving when the preconditioner is enabled.
     * \param rhs Source vector passed to the preconditioner.
     */
    void precond (VEC& lhs, VEC const& rhs) const;

    /**
     * \brief Enable or disable MLMG as a preconditioner.
     *
     * \param new_flag Set true to enable, false to disable.
     * \return Previous flag value.
     */
    bool usePrecond (bool new_flag) { return std::exchange(m_use_precond, new_flag); }

    /**
     * \brief Adjust how many MLMG smoothing iterations run inside each GMRES iteration.
     *
     * \param precond_niters Number of MLMG iterations per GMRES iteration.
     */
    void setPrecondNumIters (int precond_niters) { m_precond_niters = precond_niters; }

private:
    GM m_gmres;
    MG* m_mlmg;
    MLLinOpT<MF>* m_linop;
    int m_verbose = 0;
    int m_nlevels = 0;
    bool m_use_precond = true;
    int m_precond_niters = 1;
};

template <typename MF>
GMRESMLMGT<MF>::GMRESMLMGT (MG& mlmg)
    : m_mlmg(&mlmg), m_linop(&mlmg.getLinOp()), m_nlevels(m_linop->NAMRLevels())
{
    m_mlmg->preparePrecond();
    m_gmres.define(*this);
}

template <typename MF>
auto GMRESMLMGT<MF>::makeVecRHS () const -> VEC
{
    VEC vmf(m_nlevels);
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        vmf[ilev] = m_linop->make(ilev, 0, IntVect(0));
    }
    return vmf;
}

template <typename MF>
auto GMRESMLMGT<MF>::makeVecLHS () const -> VEC
{
    VEC vmf(m_nlevels);
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        vmf[ilev] = m_linop->make(ilev, 0, IntVect(1));
        setBndry(vmf[ilev], RT(0), 0, nComp(vmf[ilev]));
    }
    return vmf;
}

template <typename MF>
auto GMRESMLMGT<MF>::norm2 (VEC const& mf) const -> RT
{
    return m_linop->norm2Precond(GetVecOfConstPtrs(mf));
}

template <typename MF>
void GMRESMLMGT<MF>::scale (VEC& mf, RT scale_factor)
{
    for (auto& xmf : mf) {
        Scale(xmf, scale_factor, 0, nComp(xmf), 0);
    }
}

template <typename MF>
auto GMRESMLMGT<MF>::dotProduct (VEC const& mf1, VEC const& mf2) const -> RT
{
    return m_linop->dotProductPrecond(GetVecOfConstPtrs(mf1), GetVecOfConstPtrs(mf2));
}

template <typename MF>
void GMRESMLMGT<MF>::setToZero (VEC& lhs)
{
    for (auto& xmf : lhs) {
        setVal(xmf, RT(0.0));
    }
}

template <typename MF>
void GMRESMLMGT<MF>::assign (VEC& lhs, VEC const& rhs)
{
    auto nlevels = int(lhs.size());
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        LocalCopy(lhs[ilev], rhs[ilev], 0, 0, nComp(lhs[ilev]), IntVect(0));
    }
}

template <typename MF>
void GMRESMLMGT<MF>::increment (VEC& lhs, VEC const& rhs, RT a)
{
    auto nlevels = int(lhs.size());
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        Saxpy(lhs[ilev], a, rhs[ilev], 0, 0, nComp(lhs[ilev]), IntVect(0));
    }
}

template <typename MF>
void GMRESMLMGT<MF>::linComb (VEC& lhs, RT a, VEC const& rhs_a, RT b, VEC const& rhs_b)
{
    auto nlevels = int(lhs.size());
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        LinComb(lhs[ilev], a, rhs_a[ilev], 0, b, rhs_b[ilev], 0, 0, nComp(lhs[ilev]), IntVect(0));
    }
}

template <typename MF>
void GMRESMLMGT<MF>::apply (VEC& lhs, VEC const& rhs) const
{
    m_mlmg->applyPrecond(GetVecOfPtrs(lhs), GetVecOfPtrs(const_cast<VEC&>(rhs)));
}

template <typename MF>
void GMRESMLMGT<MF>::precond (VEC& lhs, VEC const& rhs) const
{
    if (m_use_precond) {
        m_mlmg->setPrecondIter(m_precond_niters);
        setToZero(lhs);
        m_mlmg->precond(GetVecOfPtrs(lhs), GetVecOfConstPtrs(rhs), 0, 0);
    } else {
        for (int ilev = 0; ilev < m_nlevels; ++ilev) {
            LocalCopy(lhs[ilev], rhs[ilev], 0, 0, nComp(lhs[ilev]), IntVect(0));
        }
    }
}

template <typename MF>
void GMRESMLMGT<MF>::solve (MF& a_sol, MF const& a_rhs, RT a_tol_rel, RT a_tol_abs)
{
    AMREX_ALWAYS_ASSERT(m_nlevels == 1);
    this->solve({&a_sol}, {&a_rhs}, a_tol_rel, a_tol_abs);
}

template <typename MF>
void GMRESMLMGT<MF>::solve (Vector<MF*> const& a_sol, Vector<MF const*> const& a_rhs, RT a_tol_rel, RT a_tol_abs)
{
    m_mlmg->incPrintIdentation();
    auto mlmg_verbose = m_mlmg->getVerbose();
    auto mlmg_bottom_verbose = m_mlmg->getBottomVerbose();
    m_mlmg->setVerbose(m_verbose);
    auto mlmg_bottom_solver = m_mlmg->getBottomSolver();

    if (mlmg_bottom_solver != BottomSolver::smoother &&
        mlmg_bottom_solver != BottomSolver::hypre &&
        mlmg_bottom_solver != BottomSolver::petsc)
    {
        m_mlmg->setBottomSolver(BottomSolver::smoother);
    }

    auto res = makeVecLHS();
    auto cor = makeVecLHS();

    m_mlmg->apply(GetVecOfPtrs(res), a_sol); // res = L(sol)
    // res = L(sol) - rhs
    bool need_to_scale_rhs = m_linop->scaleRHS(0,nullptr);
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        MF const* prhs;
        if (need_to_scale_rhs) {
            LocalCopy(cor[ilev], *a_rhs[ilev], 0, 0, nComp(cor[ilev]), IntVect(0));
            auto r = m_linop->scaleRHS(ilev, &(cor[ilev]));
            amrex::ignore_unused(r);
            prhs = &(cor[ilev]);
        } else {
            prhs = a_rhs[ilev];
        }
        Saxpy(res[ilev], RT(-1), *prhs, 0, 0, nComp(res[ilev]), IntVect(0));
    }
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        m_linop->setDirichletNodesToZero(ilev,0,res[ilev]);
    }
    m_gmres.solve(cor, res, a_tol_rel, a_tol_abs); // L(cor) = res
    // sol = sol - cor
    for (int ilev = 0; ilev < m_nlevels; ++ilev) {
        Saxpy(*a_sol[ilev], RT(-1), cor[ilev], 0, 0, nComp(*a_sol[ilev]), IntVect(0));
    }

    m_mlmg->setBottomSolver(mlmg_bottom_solver);
    m_mlmg->setVerbose(mlmg_verbose);
    m_mlmg->setBottomVerbose(mlmg_bottom_verbose);
    m_mlmg->decPrintIdentation();
}

/// \ingroup amrex_solver_gmres
using GMRESMLMG = GMRESMLMGT<MultiFab>;

}

#endif
