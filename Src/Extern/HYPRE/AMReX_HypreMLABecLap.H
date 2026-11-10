#ifndef AMREX_HYPRE_ML_ABECLAP_H_
#define AMREX_HYPRE_ML_ABECLAP_H_
#include <AMReX_Config.H>

#include <AMReX_CompositeSystemUtil.H>
#include <AMReX_FabFactory.H>
#include <AMReX_Geometry.H>
#include <AMReX_HypreIJIface.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MLMGBndry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiMask.H>

#include <HYPRE_sstruct_ls.h>

#include <limits>
#include <unordered_map>
#include <utility>

/**
 * \file AMReX_HypreMLABecLap.H
 *
 * Multilevel ABecLaplacian assembled into HYPRE's SStruct interface, including
 * coarse/fine coupling and optional SSAMG GPU solver support.
 *
 * In EB builds, given EB factories (EBSupport::full, at least one ghost
 * cell) with cut cells on any level, the composite operator of MLMG with
 * MLEBABecLap is assembled in IJ mode on all levels, except that the fluxes
 * across coarse/fine faces are the area-weighted fine fluxes, so that the
 * composite operator is conservative also where the EB crosses the
 * coarse/fine boundary (MLMG's reflux treats the cells there as regular).
 * The rows are weighted by the volume fraction, as is the right-hand side;
 * cells covered by a finer level are the EB average of their fine cells
 * (recomputed exactly after the solve, as MLMG averages down), and b on
 * coarse faces covered by a finer level is the area-weighted average of the
 * fine b, as in MLEBABecLap.  The right-hand side and the initial guess are
 * not used in EB-covered cells and in cells covered by a finer level.  The
 * factories must have the grids of the constructor.
 *
 * Like MLMG with maxorder 3 and an interpolation half width of 1 (which
 * matters in 2D only), both with and without EB; Dirichlet or Neumann
 * domain boundaries, refinement ratios 2 and 4; b at the face centers and
 * phi at the cell centers (MLEBABecLap's default locations, not its face
 * centroid and cell centroid options).  A ghost cell outside two faces
 * whose inner cells are off the level's grids, which MLEBABecLap::applyBC
 * does not fill (MLMG uses the solution's ghost value there, e.g. in the EB
 * Dirichlet stencil of a cut cell in that corner), is filled at the lowest
 * order: outside two coarse/fine faces (a convex corner of the grids) with
 * the coarse interpolation at the face; outside a physical boundary and the
 * coarse/fine boundary (or in a domain corner) with the boundary value, or
 * for Neumann the inner cell's.
 * Untested with EB: level 0 with coarse boundary data (a_coarse_bc).
 *
 * `<prefix>.eb_scaling` (default 1): with EB, scale the unknowns by
 * 1/sqrt(volume fraction), so that the rows of small cut cells count in the
 * residual norm (as in HypreMLNodeLap).  The dumped system is the scaled one.
 */

namespace amrex {

enum struct HypreSolverID {
    BoomerAMG, SSAMG
};

// single component only, cell centered only

/**
 * \brief Multilevel HYPRE solver that mirrors MLMG's ABecLaplacian structure.
 */
class HypreMLABecLap
{
public:
    HypreMLABecLap (HypreMLABecLap const&) = delete;
    HypreMLABecLap (HypreMLABecLap &&) = delete;
    HypreMLABecLap& operator= (HypreMLABecLap const&) = delete;
    HypreMLABecLap& operator= (HypreMLABecLap &&) = delete;

    /**
     * \brief Construct the multilevel HYPRE wrapper and copy the hierarchy metadata.
     *
     * \param a_geom       Per-level geometries (fine to coarse).
     * \param a_grids      Per-level BoxArrays.
     * \param a_dmap       Per-level DistributionMappings.
     * \param a_hypre_solver_id Solver type (BoomerAMG or SSAMG).
     * \param a_parmparse_prefix ParmParse prefix for runtime options.
     */
    HypreMLABecLap (Vector<Geometry> a_geom,
                    Vector<BoxArray> a_grids,
                    Vector<DistributionMapping> a_dmap,
                    HypreSolverID a_hypre_solver_id,
                    std::string a_parmparse_prefix = "hypre_mlabeclap");

    /**
     * \brief As above, with per-level factories: EB factories with cut cells
     * select the EB operator (IJ mode only).
     */
    HypreMLABecLap (Vector<Geometry> a_geom,
                    Vector<BoxArray> a_grids,
                    Vector<DistributionMapping> a_dmap,
                    HypreSolverID a_hypre_solver_id,
                    Vector<FabFactory<FArrayBox> const*> a_factory,
                    std::string a_parmparse_prefix = "hypre_mlabeclap");

    ~HypreMLABecLap ();

    //! Set HYPRE's logging verbosity (\p v).
    void setVerbose (int v) { m_verbose = v; }
    //! Override the maximum number of HYPRE iterations (\p v).
    void setMaxIter (int v) { m_maxiter = v; }
    //! Flag whether the assembled operator is singular (\p v).
    void setIsSingular (bool v) { m_is_singular = v; }

    /**
     * \brief EB operator: Dirichlet boundary condition on the EB for the
     * following setups (the data are copied; the argument order is
     * MLEBABecLap::setEBDirichlet's).
     *
     * \param a_phi  Per-level EB values; missing or null entries mean zero.
     * \param a_beta Per-level b coefficient on the EB (cell-centered).
     */
    void setEBDirichlet (Vector<MultiFab const*> const& a_phi,
                         Vector<MultiFab const*> const& a_beta);
    //! EB operator: homogeneous Dirichlet boundary condition on the EB
    void setEBHomogDirichlet (Vector<MultiFab const*> const& a_beta);

    /**
     * \brief IJ mode only: let setup() be called again with new coefficients
     * (same grids), updating the matrix values in place and keeping the
     * solver/preconditioner setup, which then only approximates the new
     * operator (see resetSolverSetup()).  With EB, the boundary condition
     * types and whether setEBDirichlet() was called must not change.
     */
    void setReuseSolverSetup (bool v) { m_reuse_setup = v; }
    //! IJ mode with setup reuse (hypre or AlgMG backend): redo the setup at
    //! the next solve when a solve with a reused setup needs more than \p n
    //! iterations (0: never; this is also `<prefix>.reuse_setup_max_iter`)
    void setReuseSetupMaxIter (int n) noexcept { m_reuse_max_iter = n; }
    //! IJ mode only: redo the solver/preconditioner setup at the next solve
    void resetSolverSetup () {
        if (m_ij_solver) { m_ij_solver->setNeedSetup(); }
        m_algmg_reset = true;
    }
    //! Throw std::runtime_error instead of aborting when the solve fails
    void setThrowException (bool t) noexcept { m_throw_exception = t; }
    //! Iterations and final relative residual of the last solve
    [[nodiscard]] int getNumIters () const noexcept { return m_num_iters; }
    [[nodiscard]] Real getFinalResidualNorm () const noexcept { return m_final_res; }

    /**
     * \brief Assemble the HYPRE operators/storage for the provided coefficients.
     *
     * \param a_ascalar   Scalar applied to the diagonal term.
     * \param a_bscalar   Scalar applied to the diffusive term.
     * \param a_acoefs    Level-by-level `a` coefficients.
     * \param a_bcoefs    Level-by-level arrays of face-centered `b` coefficients;
     *                    a missing or null entry means 1 on that level (with EB,
     *                    the coarse faces under a finer level then hold 1 too).
     * \param a_levelbcdata BC data per level; a missing or null entry means zero.
     *                    On a level with a physical boundary the data must have
     *                    a ghost cell (the boundary values are read from it).
     * \param a_lobc      Low-side boundary-condition types.
     * \param a_hibc      High-side boundary-condition types.
     * \param a_coarse_bc Optional BC data from even coarser AMR level.
     */
    void setup (Real a_ascalar, Real a_bscalar,
                Vector<MultiFab const*> const& a_acoefs,
                Vector<Array<MultiFab const*,AMREX_SPACEDIM>> const& a_bcoefs,
                Array<LinOpBCType,AMREX_SPACEDIM> const& a_lobc,
                Array<LinOpBCType,AMREX_SPACEDIM> const& a_hibc,
                Vector<MultiFab const*> const& a_levelbcdata,
                std::pair<MultiFab const*, IntVect> const& a_coarse_bc = {nullptr, IntVect(0)});

    /**
     * \brief Solve the assembled system with HYPRE (single-component).
     *
     * \param a_sol Vector of level solution MultiFabs (filled on return).
     * \param a_rhs Vector of level right-hand sides.
     * \param a_reltol Relative tolerance.
     * \param a_abstol Absolute tolerance.
     */
    void solve (Vector<MultiFab*> const& a_sol, Vector<MultiFab const*> const& a_rhs,
                Real a_reltol, Real a_abstol);

    // public for cuda

    /**
     * \brief Communicate `b` coefficients needed by level \p flev at coarse/fine interfaces.
     *
     * \param flev     Finest AMR level whose coarse/fine stencils are being prepared.
     * \param a_bcoefs Face-centered coefficients for that level (per direction).
     */
    void commBCoefs (int flev, Array<MultiFab const*,AMREX_SPACEDIM> const& a_bcoefs);
    /**
     * \brief Local portion of commBCoefs used to issue device-friendly copies.
     *
     * \param flev     AMR level whose coarse/fine data are exchanged.
     * \param a_bcoefs Face arrays for that level.
     * \param tags     Cached copy tags extracted from FabArray communications.
     */
    void commBCoefs_local (int flev, Array<MultiFab const*,AMREX_SPACEDIM> const& a_bcoefs,
                           Vector<FabArrayBase::CopyComTag> const& tags);

    //! Build the device maps of the IJ assembly; public for cuda
    void buildIJMaps ();

    //! Map the assembled triplets to the hypre matrix's value arrays
    //! (m_ij_gather); public for cuda
    void buildIJSlots (Gpu::DeviceVector<HYPRE_BigInt> const& rows,
                       Gpu::DeviceVector<HYPRE_BigInt> const& cols);
    //! Overwrite the hypre matrix's values with the triplet values
    //! (m_ij_gather); public for cuda
    void updateIJValues (Gpu::DeviceVector<Real> const& vals);

    //! Solve with the IJ system (see m_use_ij); public for cuda
    void solveIJ (Vector<MultiFab*> const& a_sol, Vector<MultiFab const*> const& a_rhs,
                  Real a_reltol, Real a_abstol);

    //! Per-box (row, column, value) triplets of the IJ assembly
    struct Triplets {
        Gpu::DeviceVector<HYPRE_BigInt> rows, cols;
        Gpu::DeviceVector<Real> vals;
    };

    //! EB operator: the triplets of level ilev's cells (see
    //! hypmlabeclap_eb_row), with b from ebBCoefs; public for cuda
    void buildEBTriplets (int ilev, Real a_ascalar, Real a_bscalar, MultiFab const* a_acoef,
                          Array<MultiFab,AMREX_SPACEDIM> const& a_bcoef,
                          Vector<Triplets>& chunks);

private:

    void addNonStencilEntriesToGraph ();

    //! EB operator: copies of the caller's EB data, checked against the grids
    void copyEBData (Vector<MultiFab const*> const& src, Vector<MultiFab>& dst) const;
    //! EB operator: b with ghost faces, averaged down as MLEBABecLap does
    [[nodiscard]] Vector<Array<MultiFab,AMREX_SPACEDIM>>
    ebBCoefs (Vector<Array<MultiFab const*,AMREX_SPACEDIM>> const& a_bcoefs) const;
    //! Scaling of the local rows' unknowns (one without EB scaling)
    [[nodiscard]] Vector<Real> rowScaling () const;

    //! Global row id of cell iv on level lev in the IJ system (-1 if the
    //! cell is not on the level's grids); periodic images are wrapped.
    [[nodiscard]] Long cellId (int lev, IntVect const& iv) const;

    int m_verbose = 0;
    int m_maxiter = 200;
    bool m_is_singular = false;

    Vector<Geometry> m_geom;
    Vector<BoxArray> m_grids;
    Vector<DistributionMapping> m_dmap;
    std::string m_parmparse_prefix;
    Vector<FabFactory<FArrayBox> const*> m_factory;
    int m_nlevels = 0;
    MPI_Comm m_comm = MPI_COMM_NULL;

    Vector<IntVect> m_ref_ratio;

    Real m_ascalar = std::numeric_limits<Real>::max();
    Real m_bscalar = std::numeric_limits<Real>::max();
    Array<LinOpBCType,AMREX_SPACEDIM> m_lobc{AMREX_D_DECL(LinOpBCType::bogus,
                                                          LinOpBCType::bogus,
                                                          LinOpBCType::bogus)};
    Array<LinOpBCType,AMREX_SPACEDIM> m_hibc{AMREX_D_DECL(LinOpBCType::bogus,
                                                          LinOpBCType::bogus,
                                                          LinOpBCType::bogus)};

    Vector<std::unique_ptr<MLMGBndry>> m_bndry;
    Vector<std::unique_ptr<BndryRegister>> m_bndry_rhs;
    Vector<iMultiFab> m_fine_masks;
    Vector<iMultiFab> m_crse_masks;

    // For coarse cells at coarse/fine interface. The vector is for AMR
    // levels.
    Vector<iMultiFab> m_c2f_offset_from; // offset for sparse coarse from-cells
    Vector<LayoutData<int>> m_c2f_total_from; // # of coarse from-cells w/ c2f entries
    Vector<iMultiFab> m_c2f_nentries; // # of non-stencil entries
    Vector<iMultiFab> m_c2f_offset_to; // offset for sparse to-cells, including fine (and coarse in 3d) cells
    Vector<LayoutData<int>> m_c2f_total_to; // total sum of non-stencil entries in a Box

    // B coefficients at coarse/fine interface
    Vector<Array<iMultiFab,AMREX_SPACEDIM>> m_offset_cf_bcoefs;
    Vector<Array<LayoutData<std::unique_ptr<Gpu::DeviceVector<Real>>>,AMREX_SPACEDIM>> m_cf_bcoefs;

#ifdef AMREX_USE_GPU
    template <class T> using HostVector = Gpu::PinnedVector<T>;
#else
    template <class T> using HostVector = Vector<T>;
#endif

    // For fine cells at coarse/fine interface. The non-stencil entries are
    // from fine to coarse. The outer vector is for AMR levels.
    Vector<HostVector<int>> m_f2c_bno; // local box number
    Vector<HostVector<IntVect>> m_f2c_cell; // fine cell
    Vector<Vector<HYPRE_Int>> m_f2c_nentries; // # of non-stencil entries
    Vector<HostVector<std::size_t>> m_f2c_offset; // offset into m_f2c_values
    Vector<HostVector<Real>> m_f2c_values; // values for non-stencil entries

    HYPRE_SStructGrid    m_ss_grid = nullptr;
    HYPRE_SStructStencil m_ss_stencil = nullptr;
    HYPRE_SStructGraph   m_ss_graph = nullptr;
    HYPRE_SStructSolver  m_ss_solver = nullptr;
    HYPRE_SStructSolver  m_ss_precond = nullptr;
    HYPRE_SStructMatrix  m_ss_A = nullptr;
    HYPRE_SStructVector  m_ss_x = nullptr;
    HYPRE_SStructVector  m_ss_b = nullptr;

    HYPRE_Solver m_solver = nullptr;

    //! Configurable ParCSR solver (e.g. FlexGMRES + BoomerAMG), used when
    //! `<prefix>.hypre_solver` is given; otherwise the built-in BoomerAMG
    std::unique_ptr<HypreIJIface> m_ij_solver;

    //! IJ mode (`<prefix>.hypre_solver` given or `backend = algmg`): one
    //! assembled matrix, no SStruct matrix or vectors
    bool m_use_ij = false;
    bool m_reuse_setup = false;
    bool m_throw_exception = false;
    //! The hypre vectors were assembled once; later solves write into them
    bool m_ij_vectors_assembled = false;
    bool m_ij_gather_failed = false; // hypre dropped an entry: no in-place updates
    //! AlgMG backend: build a new hierarchy at the next setup or solve
    bool m_algmg_reset = false;
    //! AlgMG backend: the last setup kept the hierarchy
    bool m_algmg_kept = false;
    int m_reuse_max_iter = 0;
    Vector<Vector<Long>> m_box_offset;  // first row id of each box, per level
    Long m_ij_lower = 0;
    Long m_ij_upper = -1;
    //! Targets of the non-stencil entries of each local row (until buildIJMaps)
    std::unordered_map<Long,Vector<Long>> m_ns_targets;
    //! Device maps for the IJ assembly: row id of each cell (ghost cells:
    //! same-level neighbours, -1 off the grids), and the targets of the
    //! coarse-to-fine and fine-to-coarse entries aligned with their values
    Vector<std::unique_ptr<FabArray<BaseFab<Long>>>> m_cell_id;
    Vector<LayoutData<Gpu::DeviceVector<HYPRE_BigInt>>> m_c2f_targets;
    Vector<Gpu::DeviceVector<HYPRE_BigInt>> m_f2c_rows;
    Vector<Gpu::DeviceVector<HYPRE_BigInt>> m_f2c_cols;
    Gpu::DeviceVector<HYPRE_BigInt> m_ij_rowids; // the local rows

    //! EB operator (IJ mode): the scaling of the unknowns (with one ghost
    //! cell), the coarser level's row ids and scaling on the coarsened grids
    //! of each level, the right-hand side terms of the boundary values from
    //! the last setup, and the EB Dirichlet data for the next one
    bool m_eb = false;
    bool m_eb_scaling = true;
    Vector<MultiFab> m_eb_scale;
    Vector<std::unique_ptr<FabArray<BaseFab<Long>>>> m_crse_id;
    Vector<MultiFab> m_crse_scale;
    Vector<MultiFab> m_eb_bcrhs;
    Vector<MultiFab> m_eb_beta;
    Vector<MultiFab> m_eb_phi;
    //! With setup reuse, the triplets of later setups are summed directly
    //! into a matrix's value array, each entry adding its triplets in
    //! increasing triplet order so that the sums do not depend on thread
    //! scheduling: the triplets of entry e are trip[offset[e]:offset[e+1]].
    struct EntryGather
    {
        Gpu::DeviceVector<Long> offset;
        Gpu::DeviceVector<Long> trip;
        [[nodiscard]] bool empty () const { return offset.empty(); }
        [[nodiscard]] Long numTriplets () const { return Long(trip.size()); }
        void clear () { offset.clear(); trip.clear(); }
    };
    //! Hypre backend: entries of the assembled hypre matrix, diag values
    //! first, then offd values; empty if not built
    EntryGather m_ij_gather;

    //! `<prefix>.backend = hypre` (default) or `algmg`: in IJ mode, solve the
    //! assembled composite system with hypre or with AMReX's AlgMG
    //! (`<prefix>.algmg.*`).  `backend = algmg` implies IJ mode.
    enum class Backend { hypre, algmg };
    Backend m_backend = Backend::hypre;
    int m_num_iters = 0;                     // last solve (IJ mode)
    Real m_final_res = 0.0;
    Vector<Long> m_rank_offset;              // first row of each rank, plus the total
    composite::OwnedCSR m_csr;               // owned rows (AlgMG backend, dumps)
    std::unique_ptr<SpMatrix<Real>> m_spA;
    std::unique_ptr<AlgMG<Real>> m_algmg;
    //! AlgMG backend with setup reuse: the sorted CSR pattern of the local
    //! rows and the triplets of each of its entries (empty if not built)
    Gpu::DeviceVector<Long> m_alg_row_offset;
    Gpu::DeviceVector<Long> m_alg_col_index;
    EntryGather m_alg_gather;

    //! `<prefix>.dump_matrix = N`: Matrix Market dumps of the first N solves
    //! with this dump prefix, counted over all objects of the run: for solve
    //! k, `<dump_prefix>_<k>_A.mtx`, `_b.mtx`, `_x.mtx` and `_rows.txt`
    //! (level and cell of each row).  `<prefix>.dump_prefix` (default: the
    //! ParmParse prefix).  Files with more than `<prefix>.dump_gather_max`
    //! values (matrix entries, or vector and row-map values; default 5e7)
    //! are written one per rank, as `<file>.<rank>`.
    int m_dump_max = 0;
    std::string m_dump_prefix;
    Long m_dump_gather_max = 50000000;
    void dumpSystem (int k, Gpu::DeviceVector<Real> const& v, bool solution);
    //! Throw (setThrowException) or abort with \p msg
    void fail (std::string const& msg) const;

    HypreSolverID m_hypre_solver_id = HypreSolverID::BoomerAMG;
    HYPRE_Int m_hypre_object_type = HYPRE_PARCSR;
};

}

#endif
