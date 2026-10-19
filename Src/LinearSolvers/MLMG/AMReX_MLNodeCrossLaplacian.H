#ifndef AMREX_MLNODECROSSLAPLACIAN_H_
#define AMREX_MLNODECROSSLAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_MLNodeLinOp.H>

#if (AMREX_SPACEDIM == 2)

namespace amrex {

/**
 * \file AMReX_MLNodeCrossLaplacian.H
 *
 * Nodal operator D sigma G for an exact nodal projection in 2D.
 */

/**
 * \brief Nodal operator `D sigma G` for an exact (idempotent) nodal projection.
 *
 * D is the nodal divergence of a cell-centered velocity and G the
 * cell-centered gradient of a nodal phi, the same D and G as in
 * MLNodeLaplacian::compRHS and MLNodeLaplacian::updateVelocity, with a
 * constant sigma. After solving `D sigma G phi = D u`, the velocity
 * `u - sigma G phi` has zero discrete divergence, and projecting it again
 * leaves it unchanged.
 *
 * Restrictions: 2D, a single AMR level, constant sigma, no EB. Solve with
 * MLMG::setMultigridType(MultigridType::algebraic), and with GMRESMLMG when
 * the cells are not square. Without a Dirichlet boundary, phi is unique up
 * to a constant and a checkerboard mode; the solver removes the checkerboard
 * mode.
 */
class MLNodeCrossLaplacian
    : public MLNodeLinOp
{
public:

    MLNodeCrossLaplacian () = default;
    MLNodeCrossLaplacian (const Vector<Geometry>& a_geom,
                          const Vector<BoxArray>& a_grids,
                          const Vector<DistributionMapping>& a_dmap,
                          const LPInfo& a_info = LPInfo(),
                          Real a_const_sigma = Real(1.0));
    ~MLNodeCrossLaplacian () override = default;

    MLNodeCrossLaplacian (const MLNodeCrossLaplacian&) = delete;
    MLNodeCrossLaplacian (MLNodeCrossLaplacian&&) = delete;
    MLNodeCrossLaplacian& operator= (const MLNodeCrossLaplacian&) = delete;
    MLNodeCrossLaplacian& operator= (MLNodeCrossLaplacian&&) = delete;

    /**
     * \brief Define the operator on a single AMR level.
     *
     * \param a_geom        Geometry (one level).
     * \param a_grids       Grids (one level); nodal or cell-centered.
     * \param a_dmap        Distribution mapping (one level).
     * \param a_info        LPInfo; the operator always has a single MG level.
     * \param a_const_sigma Constant coefficient sigma.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 Real a_const_sigma = Real(1.0));

    [[nodiscard]] std::string name () const override { return std::string("MLNodeCrossLaplacian"); }

    /**
     * \brief Compute the right-hand side D u for the projection.
     *
     * \param rhs Nodal right-hand side (output).
     * \param vel Cell-centered velocity with at least one ghost cell; at
     *            inflow boundaries the ghost cells hold the inflow velocity.
     *            Its ghost cells are modified.
     */
    void compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel);

    //! Subtract sigma G phi from the cell-centered velocity \p vel.
    void updateVelocity (const Vector<MultiFab*>& vel, const Vector<MultiFab const*>& sol) const;

    //! Cell-centered gradient G phi (without sigma) of the nodal \p sol.
    void compGradPhi (const Vector<MultiFab*>& grad, const Vector<MultiFab const*>& sol) const;

    void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;
    void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                 const MultiFab& fine_sol, const MultiFab& fine_rhs) final;
    void reflux (int crse_amrlev,
                 MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                 MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const final;

    void prepareForSolve () final;
    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final;
    void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    //! Removes the checkerboard mode of a singular solution.
    void postSolve (Vector<MultiFab*> const& sol) const final;

    //! Constant and checkerboard components of \p rhs for singular problems.
    Vector<Real> getSolvabilityOffset (int amrlev, int mglev,
                                       MultiFab const& rhs) const final;
    void fixSolvabilityByOffset (int amrlev, int mglev, MultiFab& rhs,
                                 Vector<Real> const& offset) const final;

    [[nodiscard]] bool supportsAlgMG () const override { return true; }
    void fillAlgMatrix (int mglev, MFIter const& mfi,
                        Array4<Long const> const& gid,
                        Array4<int const> const& lid,
                        Long* ncols, Long* cols, Real* mat) const override;
    void fillRHS (int mglev, MFIter const& mfi,
                  Array4<int const> const& lid,
                  Real* rhs,
                  Array4<Real const> const& bfab) const override;

private:

    //! The checkerboard is a null vector: no Dirichlet boundary and an even
    //! number of cells in every periodic direction.
    [[nodiscard]] bool checkerboardIsNull () const;

    Real m_const_sigma = Real(1.0);
};

}

#endif
#endif
