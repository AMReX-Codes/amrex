#ifndef AMREX_ML_TERRAIN_POISSON_H_
#define AMREX_ML_TERRAIN_POISSON_H_
#include <AMReX_Config.H>

#include <AMReX_MLCellLinOp.H>
#include <AMReX_MLTerrainPoisson_K.H>
#include <AMReX_TagParallelFor.H>

namespace amrex {

/**
 * \file AMReX_MLTerrainPoisson.H
 *
 * Cell-centered Poisson operator on a terrain-following mesh.
 */

/**
 * \brief Cell-centered Poisson operator on a terrain-following mesh.
 *
 * The mesh is uniform in computational space, and the physical height of
 * each cell corner is given by a nodal MultiFab (setZPhys).  The operator
 * is
 *
 *     L(phi) = div(A grad phi) / J,
 *
 * where grad is the physical gradient in terrain-following form, A holds
 * the face area factors (setAreas) and J is the cell volume factor
 * (setDetJ).  The flux on a Neumann domain face is zero.  MLMG solves the
 * J-scaled system div(A grad phi) = J rhs, so apply and the residuals
 * compute div(A grad phi), without the division by J.
 *
 * Only 3D, a single AMR level whose grids cover the domain, and
 * homogeneous boundary conditions are supported.  The level BC must be
 * set with setLevelBC(0, nullptr).  The operator is not symmetric; use
 * BiCGStab or the smoother as the bottom solver.  It can be used with MLMG
 * or with GMRESMLMG, which always uses the smoother as the bottom solver.
 * The semicoarsening settings in LPInfo are ignored.
 *
 * Boxes that span the whole domain in z are fastest, but boxes split in z
 * are also supported.
 *
 * \ingroup amrex_solver_mlmg
 */
class MLTerrainPoisson
    : public MLCellLinOpT<MultiFab>
{
public:

    using MF        = MultiFab;
    using RT        = Real;
    using FAB       = FArrayBox;
    using BCMode    = MLCellLinOpT<MF>::BCMode;
    using StateMode = MLCellLinOpT<MF>::StateMode;
    using Location  = MLCellLinOpT<MF>::Location;

    MLTerrainPoisson () = default;

    //! Construct the operator.  See define().
    MLTerrainPoisson (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const LPInfo& a_info = LPInfo());

    ~MLTerrainPoisson () override = default;

    MLTerrainPoisson (const MLTerrainPoisson&) = delete;
    MLTerrainPoisson (MLTerrainPoisson&&) = delete;
    MLTerrainPoisson& operator= (const MLTerrainPoisson&) = delete;
    MLTerrainPoisson& operator= (MLTerrainPoisson&&) = delete;

    /**
     * \brief Define the operator.
     *
     * \param a_geom  Geometry.  Only one AMR level is supported.
     * \param a_grids BoxArray, which must cover the domain.
     * \param a_dmap  DistributionMapping.
     * \param a_info  LPInfo.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

    /**
     * \brief Set the physical height of the cell corners.
     *
     * \param amrlev    AMR level (must be 0).
     * \param z_phys_nd Nodal MultiFab with at least one ghost node filled.
     */
    void setZPhys (int amrlev, MultiFab const& z_phys_nd);

    /**
     * \brief Set the face area factors (required).
     *
     * \param amrlev AMR level (must be 0).
     * \param area   Face MultiFabs in x, y and z.
     */
    void setAreas (int amrlev, Array<MultiFab const*,AMREX_SPACEDIM> const& area);

    /**
     * \brief Set the cell volume factor J.  If not set, J = 1.
     *
     * \param amrlev AMR level (must be 0).
     * \param detJ   Cell-centered MultiFab.
     */
    void setDetJ (int amrlev, MultiFab const& detJ);

    //! Line solver for levels whose boxes do not span z.
    enum struct ZSplitSolver { Spike, Column };

    /**
     * \brief Choose the solver used when boxes are split in z: Spike (the
     * default) or Column.
     */
    void setZSplitSolver (ZSplitSolver a_solver);

    [[nodiscard]] std::string name () const override { return std::string("MLTerrainPoisson"); }

    [[nodiscard]] bool isCrossStencil () const override { return false; }

    //! The smoother if the coarsest level covers the domain with at most
    //! 5x5 columns, else BiCGStab.
    [[nodiscard]] BottomSolver getDefaultBottomSolver () const override {
        Box const& cdomain = m_geom[0].back().Domain();
        bool const few_columns = cdomain.length(0) <= 5 && cdomain.length(1) <= 5;
        return (m_domain_covered[0] && few_columns) ? BottomSolver::smoother
                                                    : BottomSolver::bicgstab;
    }

    void prepareForSolve () override;
    [[nodiscard]] bool needsUpdate () const override {
        return m_needs_update || MLCellLinOpT<MultiFab>::needsUpdate();
    }
    void update () override;
    [[nodiscard]] bool isSingular (int amrlev) const override { return m_is_singular[amrlev]; }
    [[nodiscard]] bool isBottomSingular () const override { return m_is_singular[0]; }

    [[nodiscard]] bool scaleRHS (int amrlev, MultiFab* rhs) const override;

    //! No-op: cell-centered unknowns are never on a Dirichlet boundary.
    void setDirichletNodesToZero (int /*amrlev*/, int /*mglev*/, MultiFab& /*mf*/) const override {}

    void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode, StateMode s_mode,
                  const MLMGBndryT<MultiFab>* bndry=nullptr,
                  bool skip_fillboundary=false) const override;

    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const override;
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                  int redblack) const override;
    void FFlux (int amrlev, const MFIter& mfi, const Array<FAB*,AMREX_SPACEDIM>& flux,
                const FAB& sol, Location loc, int face_only=0) const override;

    //! Compute the terrain-following gradient, i.e., the negative of the flux.
    void compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                   MultiFab& sol, Location loc) const override;

    struct BCTag {
        Box bx;
        int local_index;

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Box const& box () const noexcept { return bx; }
    };

    [[nodiscard]] mlterrain::BCInfo bcInfo (int mglev) const;
    [[nodiscard]] IntVect coarsenRatio (int mglev) const;
    void fillZPhysGhost (MultiFab& zp, int mglev) const;
    void updateCoefs ();

private:

    Vector<int> m_is_singular;
    ZSplitSolver m_zsplit = ZSplitSolver::Spike;
    void defineLineSolvers ();
    bool m_needs_update = true;

    MultiFab m_zphys; // finest MG level only

    // Per MG level of AMR level 0
    Vector<Array<MultiFab,AMREX_SPACEDIM>> m_area;
    Vector<MultiFab> m_column_lu;
    Vector<MultiFab> m_rx;
    Vector<MultiFab> m_ry;
    Vector<MultiFab> m_zf;
    MultiFab m_detJ;

    mutable Vector<TagVector<BCTag>> m_bc_tags;

    // GPU smoother: residual scratch and unallocated column slabs.
    mutable Vector<MultiFab> m_smooth_res;
    mutable Vector<MultiFab> m_column_slab;

    // Levels whose boxes do not span z, with the same z cuts in every
    // column: block solves coupled through spikes W, V (on this layout).
    // Block end values go to x-y column slabs (k = block index), where the
    // interface values are solved and sent back.
    Vector<MultiFab> m_spike;
    Vector<MultiFab> m_spike_end;         // per column: W, V at block ends
    mutable Vector<MultiFab> m_spike_sum; // per column: y at block ends, then l, f
    mutable Vector<MultiFab> m_spike_buf; // per box: the same
    mutable Vector<MultiFab> m_spike_lu;  // per column: eliminated reduced rows
    Vector<Vector<int>> m_zcut;           // z cuts, p+1 values

    // Other split levels: whole-column layout for the line solve, with its
    // LU, residual and correction.
    Vector<MultiFab> m_col_lu;
    mutable Vector<MultiFab> m_col_res;
    mutable Vector<MultiFab> m_col_cor;
};

}

#endif
