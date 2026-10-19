#ifndef AMREX_HYPRE_ML_NODELAP_H_
#define AMREX_HYPRE_ML_NODELAP_H_
#include <AMReX_Config.H>

#include <AMReX_CompositeSystemUtil.H>
#include <AMReX_HypreIJIface.H>
#include <AMReX_MultiFab.H>

#include <memory>
#include <string>

/**
 * \file AMReX_HypreMLNodeLap.H
 *
 * Composite (all AMR levels) nodal Laplacian assembled into a single hypre
 * ParCSR matrix, reproducing the composite system that MLMG solves with
 * MLNodeLaplacian.
 */

namespace amrex {

class MLNodeLaplacian;

/**
 * \brief Composite multi-level solve of an MLNodeLaplacian with hypre.
 *
 * The composite operator of MLMG + MLNodeLaplacian is (up to a lumped-mass
 * row scaling and extra fine-equation terms in the coarse/fine rows) the
 * bilinear/trilinear finite-element operator with hanging-node constraints.
 * This class assembles the symmetric form S = Q^T K Q, where K sums the Q1
 * element matrices of all cells not covered by a finer level and Q maps
 * coarse/fine fine-level nodes to the coarse unknowns by multilinear
 * interpolation.  The right-hand side is transformed so that the solution
 * equals MLMG's composite solution (up to a constant for singular problems,
 * which are made consistent with nodal-volume weights; with an inconsistent
 * composite right-hand side this differs slightly from MLMG's correction).
 *
 * Supports the Sigma coarsening strategy with scalar, constant or mapped
 * (per-direction) sigma, Cartesian coordinates, and refinement ratios 2
 * and 4; and, in EB builds, the RAP strategy with embedded boundaries, where
 * the element matrices of cut cells are integrated over their fluid part
 * (MLNodeLaplacian's EB integrals), covered cells are skipped, and nodes
 * without fluid are not unknowns.  At coarse/fine nodes next to cut cells
 * this composite operator is the consistent finite-element one, whereas
 * MLMG's reflux treats the cells there as regular.
 *
 * The solver is configured with the HypreIJIface keys under the given
 * ParmParse prefix (default "hypre_mlnodelap").  Further keys:
 *
 * - `<prefix>.backend = hypre | algmg` (default hypre): solve the
 *   assembled system with hypre or with AMReX's AlgMG, configured by
 *   `<prefix>.algmg.*` (see composite::algmg_parse_options).
 * - `<prefix>.dump_matrix = N` (default 0): write the systems of the first
 *   N solves with this dump prefix (counted over all objects of the run)
 *   in Matrix Market format for external analysis: for solve k,
 *   `<dump_prefix>_<k>_A.mtx`, `_b.mtx`, `_x.mtx` (matrix, right-hand
 *   side, solution) and `_rows.txt` (per row: level, node index, scaling).
 *   The dumped system is the one solved, S' x' = b' with S' = s S s, b' = s b,
 *   x = s x'.  `<prefix>.dump_prefix` defaults to the ParmParse prefix.
 *   Files with more than `<prefix>.dump_gather_max` values (matrix
 *   entries, or vector and row-map values; default 5e7) are written one
 *   per rank, as `<file>.<rank>`.
 * - `<prefix>.eb_scaling` (default 1): in EB builds, scale the unknowns by
 *   1/sqrt(nodal fluid fraction) so that rows next to small cut cells count
 *   in the residual norm.
 */
class HypreMLNodeLap
{
public:
    using BigInt = HYPRE_BigInt;

    /**
     * \brief Assemble the composite matrix from a prepared linear operator.
     *
     * \param a_linop  MLNodeLaplacian with sigma set and masks built (e.g. after
     *                 an MLMG solve or MLMG::prepareForSolve).
     * \param a_prefix ParmParse prefix for the hypre solver options.
     * \param a_verbose Verbosity.
     */
    HypreMLNodeLap (MLNodeLaplacian const& a_linop,
                    std::string a_prefix = "hypre_mlnodelap",
                    int a_verbose = 0);

    ~HypreMLNodeLap ();

    HypreMLNodeLap (HypreMLNodeLap const&) = delete;
    HypreMLNodeLap (HypreMLNodeLap &&) = delete;
    HypreMLNodeLap& operator= (HypreMLNodeLap const&) = delete;
    HypreMLNodeLap& operator= (HypreMLNodeLap &&) = delete;

    /**
     * \brief Solve the composite system.
     *
     * \param a_sol Per-level nodal solution.  On input, its values on
     *              Dirichlet domain boundaries are the boundary values (as
     *              with MLMG); on output it is filled everywhere, including
     *              covered and coarse/fine nodes, like MLMG.
     * \param a_rhs Per-level nodal right-hand side as produced by
     *              MLNodeLaplacian::compRHS (MLMG's composite RHS).
     * \param a_reltol Relative tolerance.
     * \param a_abstol Absolute tolerance (<= 0: unused).
     * \param a_maxiter Iteration limit.
     */
    void solve (Vector<MultiFab*> const& a_sol, Vector<MultiFab const*> const& a_rhs,
                Real a_reltol, Real a_abstol, int a_maxiter);

    //! Throw std::runtime_error instead of aborting when the solve fails
    //! (no convergence, a hypre error, a non-finite residual, AlgMG diverging)
    void setThrowException (bool t) noexcept { m_throw_exception = t; }

    //! Iterations and final relative residual of the last solve
    [[nodiscard]] int getNumIters () const { return m_num_iters; }
    [[nodiscard]] Real getFinalResidualNorm () const { return m_final_res; }

    //! Number of composite unknowns
    [[nodiscard]] BigInt numUnknowns () const { return m_nglobal; }

    //! Internal: called by the constructor (public only for CUDA extended lambdas)
    void buildMatrix ();

private:
    using IdFab = BaseFab<Long>;
    using IdFabArray = FabArray<IdFab>;

    void numberUnknowns ();
    void loadRHS (Vector<MultiFab const*> const& a_rhs, Vector<MultiFab*> const& a_sol);
    void getSolution (Vector<MultiFab*> const& a_sol);
    void dumpSystem (int k, bool solution);
    //! Throw (setThrowException) or abort with \p msg
    void fail (std::string const& msg) const;

    //! Neumann scale (2 per Neumann/inflow domain face) of node iv on level lev
    [[nodiscard]] int neumannScale (int lev, IntVect const& iv) const;
    //! Is node iv on a Dirichlet domain face of level lev?
    [[nodiscard]] bool onDirichletFace (int lev, IntVect const& iv) const;

    [[nodiscard]] FabFactory<FArrayBox> const* factory (int lev) const;

    MLNodeLaplacian const* m_linop;
    std::string m_prefix;
    int m_verbose = 0;
    int m_nlevels = 0;
    MPI_Comm m_comm = MPI_COMM_NULL;
    int m_nprocs = 1;
    int m_myproc = 0;

    Vector<Geometry> m_geom;
    Vector<BoxArray> m_nd_grids;  // nodal BoxArrays
    Vector<DistributionMapping> m_dmap;
    Vector<int> m_ref_ratio;      // m_ref_ratio[lev]: ratio between lev and lev+1

    // Host-accessible (pinned in GPU builds) per-level data
    Vector<std::unique_ptr<IdFabArray>> m_gid;         // global id or -1
    Vector<std::unique_ptr<iMultiFab>> m_owner;        // node ownership
    Vector<std::unique_ptr<iMultiFab>> m_dmsk;         // Dirichlet (incl. c/f on fine levels)
    Vector<std::unique_ptr<iMultiFab>> m_ndfine;       // nodal fine mask (numbering only)
    Vector<std::unique_ptr<IdFabArray>> m_cgid_on_f;   // level lev-1 gid on coarsened level lev nodes
    Vector<std::unique_ptr<iMultiFab>> m_cndfine_on_f; // level lev-1 nodal fine mask, same layout
    Vector<std::unique_ptr<MultiFab>> m_ndvol;         // EB: nodal fluid volume fraction (or null)
    Vector<std::unique_ptr<MultiFab>> m_scale;         // EB: unknown scaling 1/sqrt(ndvol) (or null)
    Vector<std::unique_ptr<MultiFab>> m_cscale_on_f;   // level lev-1 scaling on coarsened level lev nodes
    bool m_eb_scaling = true;     // <prefix>.eb_scaling
    bool m_has_dirichlet = false; // any non-periodic Dirichlet domain face
    bool m_throw_exception = false;

    //! Coupling of an unknown to a Dirichlet domain node with a value
    struct DirichletCoupling { Long row; Real coef; int lev; int li; IntVect node; };
    Vector<DirichletCoupling> m_dbc;

    Vector<Long> m_rank_offset;   // first gid of each rank, plus total at the end
    BigInt m_ilower = 0;
    BigInt m_iupper = -1;
    BigInt m_nglobal = 0;
    Long m_pin_gid = -1;          // pinned unknown for singular problems

    enum class Backend { hypre, algmg };
    Backend m_backend = Backend::hypre;

    // Matrix Market dumps
    int m_dump_max = 0;
    std::string m_dump_prefix;
    Long m_dump_gather_max = 50000000;

    // Owned rows of the assembled matrix (AlgMG backend and dumps)
    composite::OwnedCSR m_csr;
    std::unique_ptr<SpMatrix<Real>> m_spA;
    std::unique_ptr<AlgMG<Real>> m_algmg;

    // Right-hand side and solution of this rank's rows
    Vector<Real> m_bvec;
    Vector<Real> m_xvec;

    std::unique_ptr<HypreIJIface> m_ij;
    HYPRE_IJMatrix m_A = nullptr;
    HYPRE_IJVector m_b = nullptr;
    HYPRE_IJVector m_x = nullptr;

    int m_num_iters = 0;
    Real m_final_res = 0.0;
    Real m_bnorm = 0.0;           // 2-norm of the assembled right-hand side
};

}

#endif
