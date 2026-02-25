#ifndef AMREX_HypreABecLap_H_
#define AMREX_HypreABecLap_H_
#include <AMReX_Config.H>

#include <AMReX_Hypre.H>
#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BndryData.H>
#include <AMReX_LO_BCTYPES.H>

#include <HYPRE_struct_ls.h>

#include <memory>

/**
 * \file AMReX_HypreABecLap.H
 *
 * Struct interface wrapper that solves cell-centered ABecLaplacian with HYPRE.
 */

namespace amrex {

/**
 * \brief Struct-interface HYPRE solver for `(a - b div beta grad) phi = rhs`.
 */
class HypreABecLap final
    : public Hypre
{
public:

    /**
     * \brief Bind the solver to a single AMR level described by \p grids/\p dmap.
     *
     * \param grids  Grid layout for the level.
     * \param dmap   DistributionMapping for \p grids.
     * \param geom_  Geometry describing the level.
     * \param comm_  MPI communicator handed to HYPRE.
     */
    HypreABecLap (const BoxArray& grids, const DistributionMapping& dmap,
                  const Geometry& geom_, MPI_Comm comm_);
    ~HypreABecLap () override;

    HypreABecLap (HypreABecLap const&) = delete;
    HypreABecLap (HypreABecLap &&) = delete;
    HypreABecLap& operator= (HypreABecLap const&) = delete;
    HypreABecLap& operator= (HypreABecLap &&) = delete;

    /**
     * \brief Solve the struct-interface linear system with the configured coefficients.
     *
     * \param soln    Solution MultiFab (filled on return).
     * \param rhs     Right-hand side MultiFab.
     * \param reltol  Relative tolerance.
     * \param abstol  Absolute tolerance.
     * \param maxiter Maximum number of iterations.
     * \param bndry   Boundary  data.
     * \param max_bndry_order Maximum order at the boundary.
     */
    void solve (MultiFab& soln, const MultiFab& rhs, Real reltol, Real abstol,
                int maxiter, const BndryData& bndry, int max_bndry_order) override;

private:

    HYPRE_StructGrid grid = nullptr;
    HYPRE_StructMatrix A = nullptr;
    HYPRE_StructVector b = nullptr;
    HYPRE_StructVector x = nullptr;
    HYPRE_StructSolver solver = nullptr;

    void getSolution (MultiFab& soln);

public: // for cuda
    //! Build the HYPRE structures and apply runtime parameters.
    void prepareSolver ();
    //! Copy AMReX solution guess \p soln and RHS \p rhs into HYPRE vectors prior to a solve.
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
};

}

#endif
