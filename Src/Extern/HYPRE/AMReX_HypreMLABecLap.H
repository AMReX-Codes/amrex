#ifndef AMREX_HYPRE_ML_ABECLAP_H_
#define AMREX_HYPRE_ML_ABECLAP_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MLMGBndry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiMask.H>

#include <HYPRE_sstruct_ls.h>

#include <limits>
#include <utility>

/**
 * \file AMReX_HypreMLABecLap.H
 *
 * Multilevel ABecLaplacian assembled into HYPRE's SStruct interface, including
 * coarse/fine coupling and optional SSAMG GPU solver support.
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

    ~HypreMLABecLap ();

    //! Set HYPRE's logging verbosity (\p v).
    void setVerbose (int v) { m_verbose = v; }
    //! Override the maximum number of HYPRE iterations (\p v).
    void setMaxIter (int v) { m_maxiter = v; }
    //! Flag whether the assembled operator is singular (\p v).
    void setIsSingular (bool v) { m_is_singular = v; }

    /**
     * \brief Assemble the HYPRE operators/storage for the provided coefficients.
     *
     * \param a_ascalar   Scalar applied to the diagonal term.
     * \param a_bscalar   Scalar applied to the diffusive term.
     * \param a_acoefs    Level-by-level `a` coefficients.
     * \param a_bcoefs    Level-by-level arrays of face-centered `b` coefficients.
     * \param a_lobc      Low-side boundary-condition types.
     * \param a_hibc      High-side boundary-condition types.
     * \param a_levelbcdata BC data per level.
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

    // update? updateDirichleBC? or updateCoeffs?

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

private:

    void addNonStencilEntriesToGraph ();

    int m_verbose = 0;
    int m_maxiter = 200;
    bool m_is_singular = false;

    Vector<Geometry> m_geom;
    Vector<BoxArray> m_grids;
    Vector<DistributionMapping> m_dmap;
    std::string m_parmparse_prefix;
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

    HypreSolverID m_hypre_solver_id = HypreSolverID::BoomerAMG;
    HYPRE_Int m_hypre_object_type = HYPRE_PARCSR;
};

}

#endif
