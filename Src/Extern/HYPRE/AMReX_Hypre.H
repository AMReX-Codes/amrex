#ifndef AMREX_HYPRE_H_
#define AMREX_HYPRE_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BndryData.H>

#include "HYPRE.h"
#include "_hypre_utilities.h"

#include <memory>

/**
 * \file AMReX_Hypre.H
 *
 * Base class and factory helpers that wrap HYPRE solvers for AMReX
 * (struct, semi-struct, and IJ interfaces).
 */

namespace amrex
{

/**
 * \brief Lightweight base class shared by HYPRE-backed linear solvers,
 * HypreABecLap, HypreABecLap2, and HypreABecLap3.
 *
 * Stores geometry, coefficients, and solver-wide options that concrete
 * implementations (Struct, SStruct, IJ) consume in their respective
 * `solve()` routines.
 */
class Hypre
{
public:

    //! HYPRE interface modes supported.
    enum class Interface : int { structed, semi_structed, ij };

    /**
     * \brief Construct a HYPRE wrapper bound to the supplied grids.
     *
     * \param grids Grid layout on the AMR level being solved.
     * \param dmap  Distribution map for \p grids.
     * \param geom  Geometry describing the domain.
     * \param comm_ MPI communicator handed to HYPRE.
     */
    Hypre (const BoxArray& grids,
           const DistributionMapping& dmap,
           const Geometry& geom,
           MPI_Comm comm_);

    virtual ~Hypre ();
    Hypre (Hypre &&) noexcept = default;
    Hypre& operator= (Hypre &&) noexcept = default;

    Hypre (Hypre const&) = delete;
    Hypre& operator= (Hypre const&) = delete;

    /**
     * \brief Store constant scalars applied to the `a` (diagonal) and `b` (Laplacian) terms.
     */
    void setScalars (Real sa, Real sb);
    /**
     * \brief Install per-cell `a` coefficients.
     */
    void setACoeffs (const MultiFab& alpha);
    /**
     * \brief Install per-face `b` coefficients for each spatial direction.
     */
    void setBCoeffs (const Array<const MultiFab*,BL_SPACEDIM>& beta);
    //! Set HYPRE verbosity level (\p _verbose).
    void setVerbose (int _verbose);
    //! Flag a singular operator (enables null-space handling in derived classes).
    void setIsMatrixSingular(bool flag) { is_matrix_singular = flag; }
    /**
     * \brief Solve the linear system using the concrete interface implementation.
     *
     * \param soln             Solution MultiFab updated in place.
     * \param rhs              Right-hand side.
     * \param rel_tol          Relative convergence tolerance.
     * \param abs_tol          Absolute convergence tolerance.
     * \param max_iter         Maximum iterations allowed.
     * \param bndry            Boundary condition data.
     * \param max_bndry_order  Maximum order at the boundary.
     */
    virtual void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol,
                        int max_iter, const BndryData& bndry, int max_bndry_order) = 0;

    /**
     * \brief Return 1 if \p i is a power of two, otherwise 0.
     */
    static HYPRE_Int ispow2 (HYPRE_Int i)
    {
        if (i == 1) {
            return 1;
        } else if ((i <= 0) || (i & 1)) {
            return 0;
        } else {
            return ispow2(i/2);
        }
    }

    //! Convenience helper that casts a box's low corner to HYPRE_Int.
    static Array<HYPRE_Int,AMREX_SPACEDIM> loV (const Box& b) {
        const auto& v = b.loVect();
        return {AMREX_D_DECL(static_cast<HYPRE_Int>(v[0]),
                             static_cast<HYPRE_Int>(v[1]),
                             static_cast<HYPRE_Int>(v[2]))};
    }

    //! Convenience helper that casts a box's high corner to HYPRE_Int.
    static Array<HYPRE_Int,AMREX_SPACEDIM> hiV (const Box& b) {
        const auto& v = b.hiVect();
        return {AMREX_D_DECL(static_cast<HYPRE_Int>(v[0]),
                             static_cast<HYPRE_Int>(v[1]),
                             static_cast<HYPRE_Int>(v[2]))};
    }


    //! Override the ParmParse namespace consulted for HYPRE runtime options (\p ns).
    void setHypreOptionsNamespace (const std::string& ns) noexcept
    { options_namespace = ns; }
    //! Toggle HYPRE's "old default" AMG setup (Falgout + modified classical interpolation) via \p l.
    void setHypreOldDefault (bool l) noexcept {old_default = l;}
    //! Choose the relaxation type handed to HYPRE (see HYPRE docs for valid codes) via \p n.
    void setHypreRelaxType (int n) noexcept {relax_type = n;}
    //! Choose the relaxation ordering (C/F, lexicographic, etc.) via \p n.
    void setHypreRelaxOrder (int n) noexcept {relax_order = n;}
    //! Set the number of pre/post smoothing sweeps HYPRE executes on each level via \p n.
    void setHypreNumSweeps (int n) noexcept {num_sweeps = n;}
    //! Adjust HYPRE's AMG strong-connection threshold used during coarsening via \p t.
    void setHypreStrongThreshold (Real t) noexcept {strong_threshold = t;}

    //! Stencil size used on regular cells.
    static constexpr HYPRE_Int regular_stencil_size = 2*AMREX_SPACEDIM + 1;
    //! Stencil size used on EB cells (3^D).
    static constexpr HYPRE_Int eb_stencil_size = AMREX_D_TERM(3, *3, *3);

protected:

    MPI_Comm comm = MPI_COMM_NULL;
    Geometry geom;

    int verbose = 0;
    bool old_default = true; // Falgout coarsening with modified classical interpolation
    int relax_type = 6;  // G-S/Jacobi hybrid relaxation
    int relax_order = 1; // uses C/F relaxation
    int num_sweeps = 2;  // Sweeeps on each level
    Real strong_threshold = Real(0.25); // HYPRE default is 0.25

    std::string options_namespace{"hypre"};

    MultiFab acoefs;
    Array<MultiFab,AMREX_SPACEDIM> bcoefs;
    Real scalar_a, scalar_b;

    MultiFab diaginv;

    FabFactory<FArrayBox> const* m_factory = nullptr;
    BndryData const* m_bndry = nullptr;
    int m_maxorder = -1;

    bool is_matrix_singular { false };
};

/**
 * \brief Factory that instantiates the requested HYPRE interface.
 *
 * \param grids         Grid layout for the solve.
 * \param dmap          Distribution map for \p grids.
 * \param geom          Geometry describing the domain.
 * \param comm_         MPI communicator used by HYPRE.
 * \param interface     Interface flavor (Struct/SStruct/IJ).
 * \param overset_mask  Optional overset mask used by IJ solvers.
 * \return Owning pointer to the chosen HYPRE implementation.
 */
[[nodiscard]] std::unique_ptr<Hypre>
makeHypre (const BoxArray& grids, const DistributionMapping& dmap,
           const Geometry& geom, MPI_Comm comm_, Hypre::Interface interface,
           const iMultiFab* overset_mask = nullptr);

}

#endif
