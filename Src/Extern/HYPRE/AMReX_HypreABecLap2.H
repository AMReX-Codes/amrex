#ifndef AMREX_HYPREABECLAP2_H_
#define AMREX_HYPREABECLAP2_H_
#include <AMReX_Config.H>

#include <AMReX_Hypre.H>

#include "_hypre_sstruct_mv.h"
#include "_hypre_sstruct_ls.h"

/**
 * \file AMReX_HypreABecLap2.H
 *
 * Semi-structured HYPRE solver for cell-centered ABecLaplacian problems.
 */

namespace amrex
{

/**
 * \brief Semi-structured ABecLaplacian operator built on HYPRE's SStruct API.
 */
class HypreABecLap2 final
    : public Hypre
{
public:

    /**
     * \brief Bind the semi-structured solver to a single AMR level.
     *
     * \param grids  Grid layout for the level.
     * \param dmap   Distribution map for \p grids.
     * \param geom_  Level geometry.
     * \param comm_  MPI communicator passed to HYPRE.
     */
    HypreABecLap2 (const BoxArray& grids, const DistributionMapping& dmap,
                   const Geometry& geom_, MPI_Comm comm_);
    ~HypreABecLap2 () override;

    HypreABecLap2 (HypreABecLap2 const&) = delete;
    HypreABecLap2 (HypreABecLap2 &&) = delete;
    HypreABecLap2& operator= (HypreABecLap2 const&) = delete;
    HypreABecLap2& operator= (HypreABecLap2 &&) = delete;

    /**
     * \brief Solve `(a - b div beta grad) phi = rhs` via HYPRE's SStruct interface.
     *
     * \param soln    Solution MultiFab written on return.
     * \param rhs     SStruct right-hand side.
     * \param reltol  Relative tolerance.
     * \param abstol  Absolute tolerance.
     * \param maxiter Maximum iterations.
     * \param bndry   Boundary data.
     * \param max_bndry_order Maximum order at the boundary.
     */
    void solve (MultiFab& soln, const MultiFab& rhs, Real reltol, Real abstol,
                int maxiter, const BndryData& bndry, int max_bndry_order) override;

private :

    HYPRE_SStructGrid     hgrid = nullptr;
    HYPRE_SStructStencil  stencil = nullptr;
    HYPRE_SStructGraph    graph = nullptr;
    HYPRE_SStructMatrix   A = nullptr;
    HYPRE_SStructVector   b = nullptr;
    HYPRE_SStructVector   x = nullptr;
    HYPRE_Solver          solver = nullptr;

    void getSolution (MultiFab& soln);

public: // for cuda
    //! Build the SStruct graph/matrix given the current coefficients.
    void prepareSolver ();
    //! Copy AMReX solution guess \p soln and RHS \p rhs into the HYPRE vectors before solving.
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
};

}

#endif
