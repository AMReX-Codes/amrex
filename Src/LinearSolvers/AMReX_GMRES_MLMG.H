#ifndef AMREX_GMRES_MLMG_H_
#define AMREX_GMRES_MLMG_H_
#include <AMReX_Config.H>

#include <AMReX_GMRES.H>
#include <AMReX_MLMG.H>
#include <AMReX_TypeTraits.H>
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
    using VEC = Vector<MF>;
    using MG = MLMGT<MF>;
    using RT = typename MG::RT; // double or float
    using GM = GMRES<Vector<MF>,GMRESMLMGT<MF>>;

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

    void solve (Vector<MF*> const& a_sol, Vector<MF const*> const& a_rhs, RT a_tol_rel, RT a_tol_abs);

    //! Sets verbosity.
    void setVerbose (int v) {
        m_verbose = v;
        m_gmres.setVerbose(v);
    }

    //! Sets the max number of iterations
    void setMaxIters (int niters) { m_gmres.setMaxIters(niters); }

    //! Gets the number of iterations.
    [[nodiscard]] int getNumIters () const { return m_gmres.getNumIters(); }

    //! Gets the 2-norm of the residual.
    [[nodiscard]] RT getResidualNorm () const { return m_gmres.getResidualNorm(); }

    //! Get the GMRES object.
    GM& getGMRES () { return m_gmres; }

    //! Make MultiFab without ghost cells
    VEC makeVecRHS () const;

    //! Make MultiFab with ghost cells and set ghost cells to zero
    VEC makeVecLHS () const;

    RT norm2 (VEC const& mf) const;

    static void scale (VEC& mf, RT scale_factor);

    RT dotProduct (VEC const& mf1, VEC const& mf2) const;

    //! lhs = 0
    static void setToZero (VEC& lhs);

    //! lhs = rhs
    static void assign (VEC& lhs, VEC const& rhs);

    //! lhs += a*rhs
    static void increment (VEC& lhs, VEC const& rhs, RT a);

    //! lhs = a*rhs_a + b*rhs_b
    static void linComb (VEC& lhs, RT a, VEC const& rhs_a, RT b, VEC const& rhs_b);

    //! lhs = L(rhs)
    void apply (VEC& lhs, VEC const& rhs) const;

    void precond (VEC& lhs, VEC const& rhs) const;

    //! Control whether or not to use MLMG as preconditioner.
    bool usePrecond (bool new_flag) { return std::exchange(m_use_precond, new_flag); }

    //! Set the number of MLMG preconditioner iterations per GMRES iteration.
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

using GMRESMLMG = GMRESMLMGT<MultiFab>;

}

#endif
