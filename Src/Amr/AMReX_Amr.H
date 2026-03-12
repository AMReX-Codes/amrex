
#ifndef AMREX_Amr_H_
#define AMREX_Amr_H_
#include <AMReX_Config.H>

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BCRec.H>
#include <AMReX_AmrCore.H>

/**
 * \file AMReX_Amr.H
 * \brief Driver that coordinates AmrLevel instances, time stepping, and output schedules.
 */

#include <iosfwd>
#include <list>
#include <memory>

namespace amrex {

class AmrLevel;
class LevelBld;
#if defined(AMREX_USE_SENSEI_INSITU) && !defined(AMREX_NO_SENSEI_AMR_INST)
class AmrInSituBridge;
#endif

/**
* \ingroup amrex_amrlevel
* \brief Manage hierarchy of levels for time-dependent AMR computations.
*
* The Amr class is designed to manage parts of the computation  which do
* not belong on a single level, like establishing and updating the hierarchy
* of levels, global timestepping, and managing the different AmrLevels
*/
class Amr
    : public AmrCore
{
    using BoundaryPointList = std::multimap< std::pair<int, int>, double >;

public:
    /**
     * \brief Construct an AMR driver bound to a factory of AmrLevel implementations.
     *
     * \param a_levelbld Factory responsible for creating level objects at runtime (must be non-null).
     */
    Amr (LevelBld* a_levelbld /* One must pass LevelBld* as an argument now*/);

    /**
     * \brief Construct with explicit problem domain metadata and a level factory.
     *
     * \param rb          Optional physical domain for level 0.
     * \param max_level_in Maximum AMR level allowed.
     * \param n_cell_in   Level-0 cell counts per direction.
     * \param coord       Coordinate system index.
     * \param a_levelbld  Factory responsible for creating AmrLevel instances.
     */
    Amr (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord,
        LevelBld* a_levelbld /* One must pass LevelBld* as an argument now*/);

    //! Amr owns global driver state—disable copying to avoid duplicate hierarchies.
    Amr (const Amr& rhs) = delete;
    //! Amr instances are not movable; they rely on pointers registered with static callbacks.
    Amr (Amr&& rhs) = delete;
    //! Copy assignment is likewise disabled.
    Amr& operator= (const Amr& rhs) = delete;
    //! Move assignment is disabled to keep ownership contracts simple.
    Amr& operator= (Amr&& rhs) = delete;

    /**
     * \brief Initialize Amr-wide defaults and parse the amr.* ParmParse options.
     *
     * This drives LevelBld::variableSetUp(), resets bookkeeping (plot/checkpoint cadence,
     * logging, etc.), and seeds shared resources before grids are created.
     */
    void InitAmr ();

    //! The destructor.
    ~Amr () override;

    /**
     * \brief Initialize hierarchy state (grid creation, initial conditions, etc.).
     *
     * \param strt_time Simulation start time.
     * \param stop_time Optional stop time for the initial solve (used by restart logic).
     */
    virtual void init (Real strt_time, Real stop_time);

    /**
     * \brief First half of the two-phase initialization used by initialInit().
     *
     * \param strt_time Simulation start time.
     * \param stop_time Optional stop time for the initial solve.
     * \param lev0_grids Optional user-supplied level-0 grids.
     * \param pmap       Optional processor map matching \p lev0_grids.
     */
    void InitializeInit (Real strt_time, Real stop_time,
                         const BoxArray* lev0_grids = nullptr, const Vector<int>* pmap = nullptr);

    //! Second half of initial initialization; pairs with InitializeInit() using times \p strt_time and \p stop_time.
    void FinalizeInit (Real strt_time, Real stop_time);

    //! Set per-level time steps to the values in \p dt_lev.
    void setDtLevel (const Vector<Real>& dt_lev) noexcept;

    //! Set the time step stored at AMR level \p lev to \p dt.
    void setDtLevel (Real dt, int lev) noexcept;

    //! Override the maximum allowable dt (typically CFL-limited) for each level using \p dt_min_in.
    void setDtMin (const Vector<Real>& dt_min_in) noexcept;

    //! Provide the number of time subcycles per level via \p ns.
    void setNCycle (const Vector<int>& ns) noexcept;

    //! Subcycle in time?
    int subCycle () const noexcept { return sub_cycle; }

    //! How are we subcycling?
    const std::string& subcyclingMode() const noexcept { return subcycling_mode; }

    /**
    * \brief What is "level" in Amr::timeStep?  This is only relevant if we are still in Amr::timeStep;
    *      it is set back to -1 on leaving Amr::timeStep.
    */
    int level_being_advanced () const noexcept { return which_level_being_advanced; }
    //! Physical time (time advanced so far).
    Real cumTime () const noexcept { return cumtime; }
    //! Set the accumulated physical time to \p t.
    void setCumTime (Real t) noexcept { cumtime = t; }
    //! Physical time this simulation started.
    Real startTime () const noexcept { return start_time; }
    //! Set the recorded start time to \p t.
    void setStartTime (Real t) noexcept { start_time = t; }
    //! Time step currently stored for level \p level.
    Real dtLevel (int level) const noexcept { return dt_level[level]; }
    //! Max time step (typically physics-limited) for level \p level.
    Real dtMin (int level) const noexcept { return dt_min[level]; }
    //! Vector of time steps (dt) for all levels.
    const Vector<Real>& dtLevel () const noexcept { return dt_level; }
    //! Number of subcycled time steps taken per coarse step on level \p level.
    int nCycle (int level) const noexcept { return n_cycle[level]; }
    //! Number of time steps completed so far on level \p lev.
    int levelSteps (int lev) const noexcept { return level_steps[lev]; }
    /**
     * \brief Override the recorded time-step count for level \p lev.
     *
     * \param lev Level index to update.
     * \param n   New time-step count.
     */
    void setLevelSteps (int lev, int n) noexcept { level_steps[lev] = n; }
    //! Which cycle number are we on for level \p lev?
    int levelCount (int lev) const noexcept { return level_count[lev]; }
    /**
     * \brief Set the current cycle count for level \p lev.
     *
     * \param lev Level index to update.
     * \param n   New cycle number.
     */
    void setLevelCount (int lev, int n) noexcept { level_count[lev] = n; }
    //! Whether to regrid right after restart
    static bool RegridOnRestart () noexcept;
    //! Interval (in steps) between regrids on level \p lev.
    int regridInt (int lev) const noexcept { return regrid_int[lev]; }
    //! Number of time steps between checkpoint files.
    int checkInt () const noexcept { return check_int; }
    //! Time between checkpoint files.
    Real checkPer() const noexcept { return check_per; }
    //! Number of time steps between plot files.
    int plotInt () const noexcept { return plot_int; }
    //! Time between plot files.
    Real plotPer () const noexcept { return plot_per; }
    //! Spacing in log10(time) of logarithmically spaced plot files
    Real plotLogPer () const noexcept { return plot_log_per; }
    //! Maximum level index to include when writing plot files.
    int plotMaxLevel () const noexcept { return plot_max_level; }
    //! Number of time steps between small plot files.
    int smallplotInt () const noexcept { return small_plot_int; }
    //! Time between plot files.
    Real smallplotPer () const noexcept { return small_plot_per; }
    //! Spacing in log10(time) of logarithmically spaced small plot files
    Real smallplotLogPer () const noexcept { return small_plot_log_per; }
    /**
    * \brief The names of state variables to output in the
    * plotfile.  They can be set using the amr.plot_vars variable
    * in a ParmParse inputs file.
    */
    static const std::list<std::string>& statePlotVars () noexcept { return state_plot_vars; }
    //! List of state variables written to the reduced “small” plot files.
    static const std::list<std::string>& stateSmallPlotVars () noexcept { return state_small_plot_vars; }
    /**
     * \brief Return true if \p name is already registered for full plotfiles.
     *
     * \param name State variable identifier.
     */
    static bool isStatePlotVar (const std::string& name);
    /**
     * \brief Return true if \p name is already registered for small plotfiles.
     *
     * \param name State variable identifier.
     */
    static bool isStateSmallPlotVar (const std::string& name);
    /**
    * \brief Register \p name as a state variable to include in full plotfiles (no-op if already present).
    *
    * \param name State variable identifier.
    */
    static void addStatePlotVar (const std::string& name);
    /**
     * \brief Append \p name to the small-plot variable list if not already present.
     *
     * \param name State variable identifier.
     */
    static void addStateSmallPlotVar (const std::string& name);
    /**
     * \brief Remove \p name from the full-plot list if present.
     *
     * \param name State variable identifier.
     */
    static void deleteStatePlotVar (const std::string& name);
    //! Clear the list of state_plot_vars.
    static void clearStatePlotVarList ();
    //! Clear the list of variables written to small plot files.
    static void clearStateSmallPlotVarList ();
    //!  Fill the list of state_plot_vars with all of the state quantities.
    static void fillStatePlotVarList ();
    //! Populate state_small_plot_vars with every state quantity registered on AmrLevel 0.
    static void fillStateSmallPlotVarList ();
    //!  Write out plotfiles (True/False)?
    static bool Plot_Files_Output ();
    /**
    * \brief The names of derived variables to output in the
    * plotfile.  They can be set using the amr.derive_plot_vars
    * variable in a ParmParse inputs file.
    */
    static const std::list<std::string>& derivePlotVars () noexcept { return derive_plot_vars; }
    //! List of derived quantities included when writing small plot files.
    static const std::list<std::string>& deriveSmallPlotVars () noexcept { return derive_small_plot_vars; }
    /**
     * \brief Return true if \p name is already registered for full derived-output.
     *
     * \param name Derived quantity identifier.
     */
    static bool isDerivePlotVar (const std::string& name) noexcept;
    /**
     * \brief Return true if \p name is already registered for small derived-output.
     *
     * \param name Derived quantity identifier.
     */
    static bool isDeriveSmallPlotVar (const std::string& name) noexcept;
    /**
    * \brief Register \p name as a derived quantity to include in full plotfiles (no-op if already present).
    *
    * \param name Derived quantity identifier.
    */
    static void addDerivePlotVar (const std::string& name);
    /**
     * \brief Append \p name to the small-plot derived list if not already present.
     *
     * \param name Derived quantity identifier.
     */
    static void addDeriveSmallPlotVar (const std::string& name);
    /**
     * \brief Remove \p name from the full derived list.
     *
     * \param name Derived quantity identifier.
     */
    static void deleteDerivePlotVar (const std::string& name);
    /**
     * \brief Remove \p name from the small derived list.
     *
     * \param name Derived quantity identifier.
     */
    static void deleteDeriveSmallPlotVar (const std::string& name);
    //! Clear the list of derive_plot_vars.
    static void clearDerivePlotVarList ();
    //! Clear the list of derived quantities written to small plot files.
    static void clearDeriveSmallPlotVarList ();
    //!  Fill the list of derive_plot_vars with all derived quantities.
    static void fillDerivePlotVarList ();
    //! Populate derive_small_plot_vars with every registered derived quantity.
    static void fillDeriveSmallPlotVarList ();

    //! Toggle whether new time steps are recomputed immediately after each regrid (\p flag=true recomputes).
    static void setComputeNewDtOnRegrid (bool flag) { compute_new_dt_on_regrid = flag; }

    //! Initialize global, Amr-wide static state (e.g., ParmParse defaults, logs).
    static void Initialize ();
    //! Release global Amr static resources allocated by Initialize().
    static void Finalize ();
    /**
     * \brief Access the AmrLevel object stored at index \p lev.
     *
     * \param lev Level index (0 = coarsest).
     */
    AmrLevel& getLevel (int lev) noexcept { return *amr_level[lev]; }
    //! Return the container of level objects (index 0..finest).
    Vector<std::unique_ptr<AmrLevel> >& getAmrLevels () noexcept;
    //! \return Total number of cells owned by all levels.
    Long cellCount () noexcept;
    //! \return Number of cells at level \p lev.
    Long cellCount (int lev) noexcept;
    //! \return Total number of grids currently allocated across all levels.
    int numGrids () noexcept;
    //! \return Number of grids at level \p lev.
    int numGrids (int lev) noexcept;
    //! \return Nonzero if the driver should continue advancing.
    int okToContinue () noexcept;
    /**
     * \brief Rebuild grid hierarchy finer than \p lbase.
     *
     * \param lbase  Base level (levels <= \p lbase are unchanged).
     * \param time   Physical time used when tagging/refining.
     * \param initial Whether this call happens during initialization (affects logging and I/O).
     */
    void regrid (int  lbase,
                 Real time,
                 bool initial = false) override;
    /**
     * \brief Regrid the hierarchy at physical time \p time without advancing any levels.
     *
     * \param time  Physical time used during tagging/refinement.
     * \param do_io When true, emit log messages describing the regrid.
     */
    void RegridOnly (Real time, bool do_io = true);
    /**
     * \brief Return true if level \p level is permitted to regrid on this step.
     *
     * \param level Level index being queried.
     */
    bool okToRegrid (int level) noexcept;
    /**
     * \brief Return the user-specified BoxArray that initializes level \p level (1-based).
     *
     * \param level Logical level index (1-based).
     */
    static const BoxArray& initialBa (int level) noexcept
        { BL_ASSERT(level-1 < initial_ba.size()); return initial_ba[level-1]; }
    //! \return Number of levels for which grids were explicitly provided via inputs.
    static int initialBaLevels () noexcept { return static_cast<int>(initial_ba.size()); }
    /**
     * \brief Advance the hierarchy by one coarse step.
     *
     * \param stop_time Optional time limit; time stepping will stop once this value is reached.
     */
    virtual void coarseTimeStep (Real stop_time);

    /**
     * \brief Advance once and report the coarse level dt that was used.
     *
     * \param stop_time Optional time limit (same semantics as coarseTimeStep()).
     * \return The dt associated with level 0 for that step.
     */
    Real coarseTimeStepDt (Real stop_time);
    /**
     * \brief Build a MultiFab containing the derived quantity named \p name.
     *
     * \param name Identifier registered with AmrLevel::derive().
     * \param time Physical time at which the quantity should be evaluated.
     * \param lev  AMR level on which to compute.
     * \param ngrow Number of ghost cells to fill in the new MultiFab.
     * \return Owning pointer to a MultiFab populated with the requested data.
     */
    std::unique_ptr<MultiFab> derive (const std::string& name,
                                      Real           time,
                                      int            lev,
                                      int            ngrow);
    /**
     * \brief Compute derived data on all levels and store it in the supplied MultiFabs.
     *
     * \param name  Identifier registered with AmrLevel::derive().
     * \param time  Physical time at which to evaluate the derived quantity.
     * \param mf    Per-level buffers that will receive the results.
     * \param dcomp Destination component inside each MultiFab \p mf entry.
     */
    void derive (const std::string&       name,
                 Real                     time,
                 const Vector<MultiFab*>& mf,
                 int                      dcomp);
    /**
     * \brief Compute derived data on all levels and return as owning MultiFabs.
     *
     * \param name  Identifier registered with AmrLevel::derive().
     * \param time  Physical time at which the derived results should be evaluated.
     * \param ngrow Number of ghost cells to fill (applied uniformly to all levels).
     * \return Vector of MultiFabs, slot \p lev containing the data for that level.
     */
    Vector<std::unique_ptr<MultiFab>> derive (const std::string& name,
                                              amrex::Real        time,
                                              int                ngrow);
    //! Name of the most recent restart checkpoint file.
    const std::string& theRestartFile () const noexcept { return restart_chkfile; }
    //! Name of the restart plotfile.
    const std::string& theRestartPlotFile () const noexcept { return restart_pltfile; }
    /**
     * \brief Return the \p i th datalog stream (opened lazily) for writing diagnostics.
     *
     * \param i Datalog index.
     */
    std::ostream& DataLog (int i);
    //! The filename backing the \p i th datalog stream.
    std::string DataLogName (int i) const noexcept { return datalogname[i]; }
    //! \return Number of datalog streams that have been opened so far.
    int NumDataLogs () noexcept;
    /**
    * \brief Compute the optimal subcycling pattern.
    * This assumes that anything less than cycle_max[i] is a valid
    * number of subcycles at level[i]. For example
    * if ref_ratio[i] = cycle_max[i] = 4, then 1,2,3,4 are all valid
    * values for n_cycles[i]
    *
    * \param n         Number of refinement levels under consideration.
    * \param best      Output array receiving cycle counts per level (size \p n).
    * \param dt_max    Maximum dt allowed at each level.
    * \param est_work  Estimated per-level work weights (same length as \p dt_max).
    * \param cycle_max Maximum allowed subcycles per level (matches \p n).
    * \return The coarse-level dt associated with the best pattern.
    */
    static Real computeOptimalSubcycling (int   n,
                                          int*  best,
                                          const Real* dt_max,
                                          const Real* est_work,
                                          const int*  cycle_max);

    /// \ingroup amrex_io
    /**
     * \brief Write a full-resolution plot file to disk using the current hierarchy state.
     */
    virtual void writePlotFile ();
    //! Return the AMR step index associated with the most recent plot file.
    int stepOfLastPlotFile () const noexcept {return last_plotfile;}
    /**
     * \brief Write the reduced “small plot” file (subset of variables/levels) to disk.
     */
    virtual void writeSmallPlotFile ();
    //! Return the AMR step index associated with the most recent small plot file.
    int stepOfLastSmallPlotFile () const noexcept {return last_smallplotfile;}
    /// \ingroup amrex_io
    /**
     * \brief Write current state into a chk* checkpoint directory.
     */
    virtual void checkPoint ();
    //! Return the AMR step index associated with the most recent checkpoint.
    int stepOfLastCheckPoint () const noexcept {return last_checkpoint;}

    //! \return Collection of BoxArrays that define the initial grid hierarchy.
    static const Vector<BoxArray>& getInitialBA() noexcept;

    /**
    * \brief Specialized version: register intersection points for faces orthogonal to X/Y.
    *
    * Define BoundaryPointLists that give the intersections of the external geometry
    * with constant (i,k) and (j,k) planes on the coarsest level.
    *
    * \param IntersectLoX Intersection list for the low-X face.
    * \param IntersectHiX Intersection list for the high-X face.
    * \param IntersectLoY Intersection list for the low-Y face.
    * \param IntersectHiY Intersection list for the high-Y face.
    */
    void setBoundaryGeometry(BoundaryPointList& IntersectLoX,
                             BoundaryPointList& IntersectHiX,
                             BoundaryPointList& IntersectLoY,
                             BoundaryPointList& IntersectHiY) noexcept
    {
        intersect_lox = IntersectLoX;
        intersect_hix = IntersectHiX;
        intersect_loy = IntersectLoY;
        intersect_hiy = IntersectHiY;
    }

    /**
    * \brief General version: register intersections for all Cartesian faces.
    *
    * Define BoundaryPointLists that give the intersections of the external geometry
    * with constant (i,k), (j,k), and (i,j) planes on the coarsest level.
    *
    * \param IntersectLoX Intersection list for the low-X face.
    * \param IntersectHiX Intersection list for the high-X face.
    * \param IntersectLoY Intersection list for the low-Y face.
    * \param IntersectHiY Intersection list for the high-Y face.
    * \param IntersectLoZ Intersection list for the low-Z face.
    * \param IntersectHiZ Intersection list for the high-Z face.
    */
    void setBoundaryGeometry(BoundaryPointList& IntersectLoX,
                             BoundaryPointList& IntersectHiX,
                             BoundaryPointList& IntersectLoY,
                             BoundaryPointList& IntersectHiY,
                             BoundaryPointList& IntersectLoZ,
                             BoundaryPointList& IntersectHiZ) noexcept
    {
        intersect_lox = IntersectLoX;
        intersect_hix = IntersectHiX;
        intersect_loy = IntersectLoY;
        intersect_hiy = IntersectHiY;
        intersect_loz = IntersectLoZ;
        intersect_hiz = IntersectHiZ;
    }

    //! \return Intersection points between the domain boundary and the low-X plane.
    BoundaryPointList& getIntersectLoX() noexcept
    {
        return intersect_lox;
    }
    //! \return Intersection points between the boundary and the high-X plane.
    BoundaryPointList& getIntersectHiX() noexcept
    {
        return intersect_hix;
    }
    //! \return Intersection points between the boundary and the low-Y plane.
    BoundaryPointList& getIntersectLoY() noexcept
    {
        return intersect_loy;
    }
    //! \return Intersection points between the boundary and the high-Y plane.
    BoundaryPointList& getIntersectHiY() noexcept
    {
        return intersect_hiy;
    }
    //! \return Intersection points between the boundary and the low-Z plane.
    BoundaryPointList& getIntersectLoZ() noexcept
    {
        return intersect_loz;
    }
    //! \return Intersection points between the boundary and the high-Z plane.
    BoundaryPointList& getIntersectHiZ() noexcept
    {
        return intersect_hiz;
    }

#ifdef AMREX_PARTICLES
    /**
     * \brief Trigger a particle redistribution pass across all levels.
     *
     * Mirrors AmrLevel::particle_redistribute() but drives the full hierarchy.
     */
    void RedistributeParticles ();
#endif

    /**
     * \brief Replace the distribution map on level \p lev with \p newdm and update metadata.
     *
     * \param lev   Level index whose data will be remapped.
     * \param newdm DistributionMapping describing the new ownership of level \p lev grids.
     */
    void InstallNewDistributionMap (int lev, const DistributionMapping& newdm);

    //! Whether plot/checkpoint writers pre-create directories (helps avoid races on some filesystems).
    static bool UsingPrecreateDirectories () noexcept;

protected:

    /**
     * \brief Initialize the grid hierarchy (called by Amr::init()).
     *
     * \param strt_time Simulation start time.
     * \param stop_time Optional time limit for the initialization phase.
     * \param lev0_grids Optional user-defined BoxArray for level 0.
     * \param pmap       Optional processor map that pairs with \p lev0_grids.
     */
    void initialInit (Real strt_time, Real stop_time,
                      const BoxArray* lev0_grids = nullptr, const Vector<int>* pmap = nullptr);
#ifndef AMREX_NO_PROBINIT
    //! Read the probin file, updating the flag \p init that tracks initialization iterations.
    void readProbinFile (int& init);
#endif
    //! Check for valid input.
    void checkInput ();
    //! Restart from a checkpoint file named \p filename.
    void restart (const std::string& filename);
    //! Define and initialize the coarsest level using optional level-0 grids \p lev0_grids and map \p pmap.
    void defBaseLevel (Real strt_time, const BoxArray* lev0_grids = nullptr, const Vector<int>* pmap = nullptr);
    //! Define and initialize refined levels at simulation time \p strt_time.
    void bldFineLevels (Real strt_time);
    //! Regrid level 0 on restart.
    virtual void regrid_level_0_on_restart ();
    /**
     * \brief Define new grid locations (called from regrid) and store them in \p new_grids.
     *
     * \param lbase       Coarsest level that may change.
     * \param time        Physical time used when tagging/refining.
     * \param new_finest  Output: index of the resulting finest level.
     * \param new_grids   Output: BoxArrays describing each new level.
     */
    void grid_places (int               lbase,
                      Real              time,
                      int&              new_finest,
                      Vector<BoxArray>& new_grids);

    /**
     * \brief Build a DistributionMapping for level \p lev at time \p time using grid layout \p ba.
     *
     * \return Load-balanced DistributionMapping suitable for the supplied BoxArray.
     */
    DistributionMapping makeLoadBalanceDistributionMap (int lev, Real time, const BoxArray& ba) const;
    /**
     * \brief Compute and install a new distribution map for level 0.
     *
     * \param time Simulation time at which load-balancing decisions are made.
     */
    void LoadBalanceLevel0 (Real time);

    /**
     * \brief Default error-estimation driver; forwards to the levelbld-defined callbacks.
     *
     * \param lev   Level being tagged.
     * \param tags  Tag arrays populated by derived routines.
     * \param time  Physical time passed to tagging functors.
     * \param ngrow Number of ghost cells available in \p tags.
     */
    void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override;
    //! Return cached regions that must not be tagged on level \p lev.
    BoxArray GetAreaNotToTag (int lev) override;
    /**
     * \brief Allow derived classes to override tagging manually on level \p lev.
     *
     * \param lev     Level being tagged.
     * \param tags    Tag arrays that may be modified in-place.
     * \param bf_lev  Buffer sizes supplied by Amr::regrid().
     */
    void ManualTagsPlacement (int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev) override;

    /**
     * \brief Do a single time step on level \p level.
     *
     * \param level     Level index to advance.
     * \param time      Physical time at the beginning of the step.
     * \param iteration Which sub-iteration of this coarse step we are in.
     * \param niter     Total number of iterations that will run on this level.
     * \param stop_time Optional physical time limit that halts advancement when reached.
     */
    virtual void timeStep (int  level,
                           Real time,
                           int  iteration,
                           int  niter,
                           Real stop_time);

    // pure virtual function in AmrCore
    void MakeNewLevelFromScratch (int /*lev*/, Real /*time*/, const BoxArray& /*ba*/, const DistributionMapping& /*dm*/) override
        { amrex::Abort("How did we get here!"); }
    void MakeNewLevelFromCoarse (int /*lev*/, Real /*time*/, const BoxArray& /*ba*/, const DistributionMapping& /*dm*/) override
        { amrex::Abort("How did we get here!"); }
    void RemakeLevel (int /*lev*/, Real /*time*/, const BoxArray& /*ba*/, const DistributionMapping& /*dm*/) override
        { amrex::Abort("How did we get here!"); }
    void ClearLevel (int /*lev*/) override
        { amrex::Abort("How did we get here!"); }

    //! Whether checkpoint criteria are satisfied on this step.
    bool checkPointNow () noexcept;
    //! Whether to write a plotfile now
    bool writePlotNow () noexcept;
    //! Whether to emit a small plot file on this step.
    bool writeSmallPlotNow () noexcept;

    /**
     * \brief Print per-level grid metadata to \p os for levels [\p min_lev, \p max_lev].
     *
     * \param os      Output stream receiving the report.
     * \param min_lev Smallest AMR level to describe.
     * \param max_lev Largest AMR level to describe.
     */
    void printGridInfo (std::ostream& os,
                        int           min_lev,
                        int           max_lev);

    //! Record the log filename \p filename used when dumping grid info metadata.
    void setRecordGridInfo (const std::string& filename);

    //! Record the log filename \p filename used when writing detailed run information.
    void setRecordRunInfo (const std::string& filename);

    //! Record the log filename \p filename used for terse run information.
    void setRecordRunInfoTerse (const std::string& filename);

    //! Configure the \p i th datalog stream to mirror file \p filename.
    void setRecordDataInfo (int i, const std::string& filename);

    //! Parse ParmParse state and initialize subcycling-related arrays/flags.
    void initSubcycle();
    //! Initialize plot/checkpoint cadence, names, and counters from inputs.
    void initPltAndChk();

    //! Set up in-situ analysis hooks; returns 0 on success.
    static int initInSitu();
    //! Drive any per-step in-situ actions; returns 0 on success.
    int updateInSitu();
    //! Tear down in-situ infrastructure; returns 0 on success.
    static int finalizeInSitu();

    //
    // The data ...
    //
    std::string      regrid_grids_file;   //!< Grids file that will bypass regridding.
    std::string      initial_grids_file;  //!< Grids file that will bypass regridding only at initialization.
    Vector<std::unique_ptr<AmrLevel> > amr_level;    //!< Vector of levels
    Real             cumtime = std::numeric_limits<Real>::lowest();      //!< Physical time variable.
    Real             start_time = std::numeric_limits<Real>::lowest();   //!< Physical time this simulation started.
    Vector<Real>      dt_level;     //!< Timestep at this level.
    Vector<int>       level_steps;  //!< Number of time steps at this level.
    Vector<int>       level_count;
    Vector<int>       n_cycle;
    std::string      subcycling_mode; //!<Type of subcycling to use.
    Vector<Real>      dt_min;
    Vector<int>       regrid_int;      //!< Interval between regridding.
    int              last_checkpoint; //!< Step number of previous checkpoint.
    int              check_int;       //!< How often checkpoint (# time steps).
    Real             check_per;       //!< How often checkpoint (units of time).
    std::string      check_file_root; //!< Root name of checkpoint file.
    int              last_plotfile;   //!< Step number of previous plotfile.
    int              last_smallplotfile;   //!< Step number of previous small plotfile.
    int              plot_int;        //!< How often plotfile (# of time steps)
    Real             plot_per;        //!< How often plotfile (in units of time)
    Real             plot_log_per;    //!< How often plotfile (in units of log10(time))
    int              plot_max_level;  //!< Maximum AMR level to write to a plotfile
    int              small_plot_int;  //!< How often small plotfile (# of time steps)
    Real             small_plot_per;  //!< How often small plotfile (in units of time)
    Real             small_plot_log_per;  //!< How often small plotfile (in units of log10(time))
    bool             write_plotfile_with_checkpoint;  //!< Write out a plotfile whenever we checkpoint
    int              file_name_digits; //!< How many digits to use in the plotfile and checkpoint names
    int              message_int;     //!< How often checking messages touched by user, such as "stop_run"
    std::string      plot_file_root;  //!< Root name of plotfile.
    std::string      small_plot_file_root;  //!< Root name of small plotfile.

    int              which_level_being_advanced = -1; //!< Only >=0 if we are in Amr::timeStep(level,...)

    int              record_grid_info;
    int              record_run_info;
    int              record_run_info_terse;
    std::ofstream    gridlog;
    std::ofstream    runlog;
    std::ofstream    runlog_terse;
    Vector<std::unique_ptr<std::fstream> > datalog;
    Vector<std::string> datalogname;
    int              sub_cycle;
    std::string      restart_chkfile;
    std::string      restart_pltfile;
#ifndef AMREX_NO_PROBINIT
    std::string      probin_file;
#endif
    LevelBld*        levelbld;
    bool             abort_on_stream_retry_failure;
    int              stream_max_tries;
    int              loadbalance_with_workestimates;
    int              loadbalance_level0_int;
    Real             loadbalance_max_fac;

    bool             bUserStopRequest;

    //
    // The static data ...
    //
    static std::list<std::string> state_plot_vars;  //!< State Vars to dump to plotfile
    static std::list<std::string> state_small_plot_vars;  //!< State Vars to dump to small plotfile
    static std::list<std::string> derive_plot_vars; //!< Derived Vars to dump to plotfile
    static std::list<std::string> derive_small_plot_vars; //!< Derived Vars to dump to small plotfile
    static bool                   first_plotfile;
    //! Array of BoxArrays read in to initially define grid hierarchy
    static Vector<BoxArray> initial_ba;
    //! Array of BoxArrays read in to externally define grid hierarchy at each regrid
    static Vector<BoxArray> regrid_ba;
    static bool compute_new_dt_on_regrid;

#if defined(AMREX_USE_SENSEI_INSITU) && !defined(AMREX_NO_SENSEI_AMR_INST)
    static AmrInSituBridge *insitu_bridge;
#endif

public:
    BoundaryPointList intersect_lox;
    BoundaryPointList intersect_loy;
    BoundaryPointList intersect_loz;
    BoundaryPointList intersect_hix;
    BoundaryPointList intersect_hiy;
    BoundaryPointList intersect_hiz;

    static bool first_smallplotfile;

private:
    void writePlotFileDoit (std::string const& pltfile, bool regular);
};

}

#endif /*_Amr_H_*/
