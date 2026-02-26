#ifndef AMREX_PETSC_H_
#define AMREX_PETSC_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_BndryData.H>

// Do NOT include any other PETSc headers here.  They should be included in
// .cpp files before AMReX_PETSc.H is included.
#include <petscsystypes.h>

/**
 * \file AMReX_PETSc.H
 *
 * PETSc-backed ABecLaplacian solver hooks (constructor, coefficients, and helpers).
 */

namespace amrex {

class amrex_KSP;
class amrex_Mat;
class amrex_Vec;

/**
 * \brief PETSc wrapper that assembles/solves `(a alpha - b div beta grad) phi = rhs`.
 */
class PETScABecLap
{
public:

    /**
     * \brief Construct a PETSc solver bound to the supplied AMR level.
     *
     * \param grids  Grid layout for the level.
     * \param dmap   Distribution map associated with \p grids.
     * \param geom_  Geometry describing the domain.
     * \param comm_  MPI communicator passed to PETSc.
     */
    PETScABecLap (const BoxArray& grids, const DistributionMapping& dmap,
                  const Geometry& geom_, MPI_Comm comm_);

    ~PETScABecLap ();

    PETScABecLap (PETScABecLap const&) = delete;
    PETScABecLap (PETScABecLap &&) noexcept = delete;
    PETScABecLap& operator= (PETScABecLap const&) = delete;
    PETScABecLap& operator= (PETScABecLap &&) noexcept = delete;

    /**
     * \brief Store constant scalars for the diagonal and diffusive terms.
     *
     * \param sa Scalar applied to the `a` term.
     * \param sb Scalar applied to the `b` term.
     */
    void setScalars (Real sa, Real sb);
    /**
     * \brief Install cell-centered `a` coefficients (copied as needed).
     *
     * \param alpha MultiFab containing the `a` coefficients.
     */
    void setACoeffs (const MultiFab& alpha);
    /**
     * \brief Install per-direction face-centered `b` coefficients.
     *
     * \param beta Array of per-direction coefficient MultiFabs.
     */
    void setBCoeffs (const Array<const MultiFab*,BL_SPACEDIM>& beta);
    //! Set PETSc verbosity level forwarded to KSP monitors (\p _verbose).
    void setVerbose (int _verbose);
    /**
     * \brief Solve the assembled ABec system using PETSc/KSP.
     *
     * \param soln    Solution MultiFab (updated in place).
     * \param rhs     Right-hand side.
     * \param rel_tol Relative tolerance.
     * \param abs_tol Absolute tolerance.
     * \param max_iter Maximum allowed iterations.
     * \param bndry   Boundary data.
     * \param max_bndry_order Maximum order at the boundary.
     */
    void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol,
                int max_iter, const BndryData& bndry, int max_bndry_order);

#ifdef AMREX_USE_EB
    //! Provide `b` coefficients stored in \p beb for embedded-boundary faces.
    void setEBDirichlet (MultiFab const* beb) { m_eb_b_coeffs = beb; }
#endif

private:

    static constexpr PetscInt regular_stencil_size = 2*AMREX_SPACEDIM + 1;
    static constexpr PetscInt eb_stencil_size = AMREX_D_TERM(3, *3, *3);

    MPI_Comm comm = MPI_COMM_NULL;
    Geometry geom;

    int verbose = 0;

    MultiFab acoefs;
    Array<MultiFab,AMREX_SPACEDIM> bcoefs;
    Real scalar_a, scalar_b;

    MultiFab diaginv;

    FabFactory<FArrayBox> const* m_factory = nullptr;
    BndryData const* m_bndry = nullptr;
    int m_maxorder = -1;

    std::unique_ptr<amrex_KSP> solver;
    std::unique_ptr<amrex_Mat> A;
    std::unique_ptr<amrex_Vec> b;
    std::unique_ptr<amrex_Vec> x;

    LayoutData<PetscInt> ncells_grid;
    FabArray<BaseFab<PetscInt> > cell_id;
    FabArray<BaseFab<PetscInt> > cell_id_vec;

    MultiFab const* m_eb_b_coeffs = nullptr;

public: // for cuda
    //! Build PETSc matrices/vectors and apply runtime parameters.
    void prepareSolver ();
    //! Copy AMReX solution guess \p soln and RHS \p rhs into PETSc vectors.
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    //! Copy PETSc's solved solution back into AMReX storage \p soln.
    void getSolution (MultiFab& soln);
};

/**
 * \brief Factory helper that instantiates a PETSc ABec Laplacian on one level.
 *
 * \param grids  Grid layout for the level.
 * \param dmap   Distribution mapping.
 * \param geom   Level geometry.
 * \param comm_  MPI communicator passed to PETSc.
 * \return Owning pointer to the configured PETSc solver.
 */
[[nodiscard]] std::unique_ptr<PETScABecLap>
makePetsc (const BoxArray& grids, const DistributionMapping& dmap,
           const Geometry& geom, MPI_Comm comm_);

}

#endif
