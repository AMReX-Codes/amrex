#ifndef BL_AMRCORE_H_
#define BL_AMRCORE_H_
#include <AMReX_Config.H>

#include <AMReX_AmrMesh.H>

/**
 * \file AMReX_AmrCore.H
 * \brief Base class for AMR applications that manages mesh hierarchy but not state data.
 */

#include <iosfwd>
#include <memory>

namespace amrex {

#ifdef AMREX_PARTICLES
class AmrParGDB;
#endif

/**
 * \ingroup amrex_amrcore
 * \brief Provide basic functionalities to set up an AMR hierarchy
 *
 * Because it does not own any floating point data, AmrCore used
 * virtual functions to allocate, initialize and delete data.  It also
 * requires the derived class to tag cells for refinement.
 */
class AmrCore
    : public AmrMesh
{
public:

    //! Construct an empty AmrCore.
    AmrCore ();

    /**
     * \brief Construct using legacy inputs that mirror AmrMesh.
     *
     * \param rb          Optional physical domain for level 0 (nullptr uses defaults).
     * \param max_level_in Maximum AMR level allowed.
     * \param n_cell_in   Number of cells at level 0 in each direction.
     * \param coord     Coordinate system index (defaults to Cartesian if not
     *                  specified in ParmParse parameter geometry.coord_sys).
     * \param ref_ratios  Per-level refinement ratios (level 0 onward).
     * \param is_per    Optional periodicity flags (defaults to non-periodic if
     *                  not specified in ParmParse parameter geometry.is_periodic).
     */
    AmrCore (const RealBox* rb, int max_level_in,
             const Vector<int>& n_cell_in, int coord=-1,
             Vector<IntVect> ref_ratios = Vector<IntVect>(),
             const int* is_per = nullptr);

    /**
     * \brief Construct with an explicit domain and periodicity array.
     *
     * \param rb          Physical domain for level 0.
     * \param max_level_in Maximum AMR level allowed.
     * \param n_cell_in   Number of cells at level 0 in each direction.
     * \param coord       Coordinate system index.
     * \param ref_ratios  Per-level refinement ratios (level 0 onward).
     * \param is_per      Array of periodicity flags (0 = nonperiodic, 1 = periodic).
     */
    AmrCore (const RealBox& rb, int max_level_in,
             const Vector<int>& n_cell_in, int coord,
             Vector<IntVect> const& ref_ratios,
             Array<int,AMREX_SPACEDIM> const& is_per);

    /**
     * \brief Construct from a ready-made Geometry and AmrInfo bundle.
     *
     * \param level_0_geom Geometry describing level 0.
     * \param amr_info     Configuration struct controlling level counts and limits.
     */
    AmrCore (Geometry const& level_0_geom, AmrInfo const& amr_info);

    AmrCore (AmrCore&& rhs) noexcept;
    AmrCore& operator= (AmrCore&& rhs) noexcept;

    AmrCore (const AmrCore& rhs) = delete;
    AmrCore& operator= (const AmrCore& rhs) = delete;

    ~AmrCore () override;

#ifdef AMREX_PARTICLES
    //! Return AmrParGDB pointer for particle containers.
    [[nodiscard]] AmrParGDB* GetParGDB () const noexcept { return m_gdb.get(); }
#endif

    /**
     * \brief Initialize BoxArray, DistributionMapping, and level data from scratch.
     *
     * Calling this function requires the derived class to implement MakeNewLevelFromScratch()
     * so it can allocate and initialize data on every level. After initialization completes,
     * most applications average the fine data back down to the coarse levels.
     *
     * \param time Simulation time associated with the freshly initialized hierarchy.
     */
    void InitFromScratch (Real time);

    /**
     * \brief Rebuild all levels finer than \p lbase at simulation time \p time.
     *
     * \param lbase    Coarse level after which to regenerate hierarchy.
     * \param time     Physical time used when calling user callbacks.
     * \param initial  Whether this regrid corresponds to initial creation.
     */
    virtual void regrid (int lbase, Real time, bool initial=false);

    /**
     * \brief Print a human-readable summary of grid sizes between \p min_lev and \p max_lev.
     *
     * \param os      Destination stream.
     * \param min_lev First level to report.
     * \param max_lev Last level to report (inclusive).
     */
    void printGridSummary (std::ostream& os, int min_lev, int max_lev) const noexcept;

protected:

    /**
     * \brief Tag cells for refinement on level \p lev.
     *
     * \param lev   Level whose valid region is being tagged.
     * \param tags  Tag array built on level \p lev grids.
     * \param time  Simulation time used by tagging logic.
     * \param ngrow Number of ghost layers in \p tags that must be filled.
     */
    void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override = 0;

    /**
     * \brief Create and fill a new level from scratch.
     *
     * Only invoked during initialization; derived classes must allocate and populate
     * all data for level \p lev using the supplied metadata.
     */
    void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override = 0;

    /**
     * \brief Make a new level using provided metadata and populate it from the next coarser level.
     *
     * \param lev Level index being created.
     * \param time Simulation time when the level is initialized.
     * \param ba  BoxArray describing the new grids.
     * \param dm  Distribution mapping describing ownership of those grids.
     */
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    /**
     * \brief Remake an existing level using provided metadata and refill from old data.
     *
     * \param lev Level index being regenerated.
     * \param time Simulation time to stamp on the data.
     * \param ba  New BoxArray layout for that level.
     * \param dm  Distribution mapping describing ownership.
     */
    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    /**
     * \brief Delete level \p lev data structures owned by the derived class.
     *
     * \param lev Level index whose data should be released.
     */
    virtual void ClearLevel (int lev) = 0;

#ifdef AMREX_PARTICLES
    //! Geometry, DistributionMapping and BoxArray for particle containers.
    std::unique_ptr<AmrParGDB> m_gdb;
#endif

private:
    //! Common helper invoked by constructors to initialize geometry and bookkeeping.
    void InitAmrCore ();
};

}

#endif
