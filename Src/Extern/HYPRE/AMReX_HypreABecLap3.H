#ifndef AMREX_HYPREABECLAP3_H_
#define AMREX_HYPREABECLAP3_H_
#include <AMReX_Config.H>

#include <AMReX_Hypre.H>

#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_HypreIJIface.H>

#include <algorithm>

/**
 * \file AMReX_HypreABecLap3.H
 *
 * IJ-interface ABecLaplacian that supports overset and EB-aware matrices.
 */

namespace amrex
{

/**
 * \brief HYPRE IJ solver for ABecLaplacian (optionally with overset masks).
 */
class HypreABecLap3 final
    : public Hypre
{
public:

    /**
     * \brief Construct the IJ-based solver and capture overset/EB context.
     *
     * \param grids         Grid layout for the AMR level.
     * \param dmap          Distribution map for \p grids.
     * \param geom_         Geometry describing the domain.
     * \param comm_         MPI communicator forwarded to HYPRE.
     * \param overset_mask  Optional overset mask used to drop covered cells.
     */
    HypreABecLap3 (const BoxArray& grids, const DistributionMapping& dmap,
                   const Geometry& geom_, MPI_Comm comm_,
                   const iMultiFab* overset_mask = nullptr);
    ~HypreABecLap3 () override;

    HypreABecLap3 (HypreABecLap3 const&) = delete;
    HypreABecLap3 (HypreABecLap3 &&) = delete;
    HypreABecLap3& operator= (HypreABecLap3 const&) = delete;
    HypreABecLap3& operator= (HypreABecLap3 &&) = delete;

    /**
     * \brief Solve with the IJ wrapper (supports overset/EB contributions).
     *
     * \param soln    Solution MultiFab (updated in place).
     * \param rhs     Right-hand side data.
     * \param rel_tol Relative tolerance.
     * \param abs_tol Absolute tolerance.
     * \param max_iter Maximum iterations.
     * \param bndry   Boundary data.
     * \param max_bndry_order Maximum order at the boundary.
     */
    void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol,
                int max_iter, const BndryData& bndry, int max_bndry_order) override;

#ifdef AMREX_USE_EB
    //! Provide `b` coefficients stored in \p beb for embedded-boundary faces.
    void setEBDirichlet (MultiFab const* beb) { m_eb_b_coeffs = beb; }
#endif

private :
    std::unique_ptr<HypreIJIface> hypre_ij;

    // Non-owning references to HYPRE matrix, rhs, and solution data
    HYPRE_IJMatrix A = nullptr;
    HYPRE_IJVector b = nullptr;
    HYPRE_IJVector x = nullptr;

    LayoutData<HYPRE_Int> ncells_grid;
    FabArray<BaseFab<HYPRE_Int> > cell_id;
    FabArray<BaseFab<HYPRE_Int> > cell_id_vec;

    MultiFab const* m_eb_b_coeffs = nullptr;

    iMultiFab const* m_overset_mask = nullptr;

public: // for CUDA
    //! Build IJ structures and configure the HYPRE solver/preconditioner.
    void prepareSolver ();
    //! Copy HYPRE solution back into AMReX storage (\p soln).
    void getSolution (MultiFab& soln);
    //! Load RHS \p rhs and initial guess \p soln into HYPRE vectors.
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
};

}

#endif
