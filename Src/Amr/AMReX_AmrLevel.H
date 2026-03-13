
#ifndef AMREX_AmrLevel_H_
#define AMREX_AmrLevel_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_LayoutData.H>
#include <AMReX_Derive.H>
#include <AMReX_BCRec.H>
#include <AMReX_Amr.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_StateData.H>
#include <AMReX_VisMF.H>
#include <AMReX_RungeKutta.H>
#include <AMReX_FillPatcher.H>

/**
 * \file AMReX_AmrLevel.H
 * \brief Virtual base class that owns data and time-stepping logic for a single AMR level.
 */
#ifdef AMREX_USE_EB
#include <AMReX_EBSupport.H>
#endif

#include <memory>
#include <map>

namespace amrex {

class TagBox;
class TagBoxArray;

/**
* \ingroup amrex_amrlevel
* \brief Virtual base class for managing individual levels.
* AmrLevel functions both as a container for state data on a level
* and also manages the advancement of data in time.
*/
class AmrLevel
{
    friend class Amr;
    friend class FillPatchIterator;
    friend class FillPatchIteratorHelper;

public:
    //! What time are we at?
    enum TimeLevel { AmrOldTime,
                     AmrHalfTime,
                     AmrNewTime,
                     Amr1QtrTime,
                     Amr3QtrTime,
                     AmrOtherTime };
    //! The destructor.
    virtual ~AmrLevel ();

    //! Disable copying because AmrLevel owns state data and geometry caches.
    AmrLevel (const AmrLevel&) = delete;
    //! Disable moving; derived classes manage non-movable FillPatch caches.
    AmrLevel (AmrLevel&&) = delete;
    //! Copy assignment is undefined for the same ownership reasons.
    AmrLevel& operator= (const AmrLevel&) = delete;
    //! Move assignment is also disabled—use factory methods instead.
    AmrLevel& operator= (AmrLevel&&) = delete;

    /**
     * \brief Derive the directory name for this level inside a plot/checkpoint tree.
     *
     * \param dir       Root directory (typically plot or checkpoint folder).
     * \param LevelDir  Output: relative level directory (e.g., Level_3).
     * \param FullPath  Output: absolute path combining \p dir and \p LevelDir.
     */
    void LevelDirectoryNames (const std::string &dir,
                              std::string &LevelDir,
                              std::string &FullPath) const;
    /**
     * \brief Create the Level_* directory inside \p dir (if needed) before writing plots/checkpoints.
     *
     * \param dir Base output directory (plot or checkpoint root).
     */
    virtual void CreateLevelDirectory (const std::string &dir);
    /**
     * \brief Record whether CreateLevelDirectory() has already created the folder.
     *
     * \param ldc True if the level directory exists; false to clear the flag.
     */
    void SetLevelDirectoryCreated(bool ldc) noexcept { levelDirectoryCreated = ldc; }
    /**
    * \brief A string written as the first item in writePlotFile() at
    * level zero.
    * It is so we can distinguish between different types of
    * plot files.
    * This default "HyperCLaw-V1.1" is for VisIt software and some of our
    * internal postprocessing routines
    */
    virtual std::string thePlotFileType () const
    {
        static const std::string the_plot_file_type("HyperCLaw-V1.1"); //!< Default identifier consumed by VisIt and AMReX tools.
        return the_plot_file_type;
    }
    /**
     * \brief Write plot file data for this level into directory \p dir.
     *
     * \param dir Output root (e.g., plt00010).
     * \param os  Metadata stream shared across levels.
     * \param how VisMF write mode (NFiles or OneFilePerCPU).
     */
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how = VisMF::NFiles);

    /**
     * \brief Perform pre-plotfile operations (e.g., gather metadata) for directory \p dir using \p os.
     *
     * \param dir Plotfile directory about to be written.
     * \param os  Stream receiving metadata (often the Header).
     */
    virtual void writePlotFilePre (const std::string& dir,
                                   std::ostream&      os);

    /**
     * \brief Perform post-plotfile cleanup for directory \p dir once writes complete.
     *
     * \param dir Plotfile directory that was just written.
     * \param os  Metadata stream (same one passed to writePlotFilePre()).
     */
    virtual void writePlotFilePost (const std::string& dir,
                                    std::ostream&      os);

    /**
     * \brief Write a reduced plot file (same parameters as writePlotFile()) to \p dir.
     *
     * \param dir Target directory.
     * \param os  Metadata stream.
     * \param how VisMF write mode.
     */
    virtual void writeSmallPlotFile (const std::string& dir,
                                     std::ostream&      os,
                                     VisMF::How         how = VisMF::NFiles)
    {
        (void)dir;
        (void)os;
        (void)how;
    }
    /**
     * \brief Write current state to a checkpoint directory.
     *
     * \param dir       Destination directory (e.g., chk00010).
     * \param os        Metadata stream shared across levels.
     * \param how       VisMF write mode.
     * \param dump_old  Whether to include the "old" state data.
     */
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             VisMF::How         how = VisMF::NFiles,
                             bool               dump_old = true);
    /**
     * \brief Do pre-checkpoint work (e.g., directory creation) for \p dir using \p os.
     *
     * \param dir Checkpoint directory being created.
     * \param os  Metadata stream for header output.
     */
    virtual void checkPointPre (const std::string& dir,
                                std::ostream&      os);
    /**
     * \brief Do post-checkpoint cleanup for \p dir once writes complete.
     *
     * \param dir Checkpoint directory that was just written.
     * \param os  Metadata stream used during writing.
     */
    virtual void checkPointPost (const std::string& dir,
                                 std::ostream&      os);
    /**
     * \brief Populate this level from checkpoint input stream \p is.
     *
     * \param papa          Owning Amr driver (provides geometry and config).
     * \param is            Input stream pointing at the level's checkpoint data.
     * \param bReadSpecial  Whether to parse additional application-specific fields.
     */
    virtual void restart (Amr&          papa,
                          std::istream& is,
                          bool          bReadSpecial = false);

    /**
     * \brief Hook to adjust which state descriptors are expected in checkpoint streams.
     *
     * \param state_in_checkpoint Array of flags (one per descriptor) that this level may modify.
     */
    virtual void set_state_in_checkpoint (Vector<int>& state_in_checkpoint);

    /**
     * \brief Return true if \p name matches a registered state variable, also providing descriptor info.
     *
     * \param name        State variable identifier.
     * \param state_indx  Output: descriptor index when found.
     * \param ncomp       Output: number of components associated with that state.
     */
    static bool isStateVariable (const std::string& name,
                                int&               state_indx,
                                int&               ncomp);

    //! Resize/clear cached FillPatchIterator state in preparation for new fill operations.
    static void FlushFPICache ();
    /**
     * \brief Compute the initial time step size for each level.
     *
     * \param finest_level Finest AMR level existing so far.
     * \param sub_cycle    Whether subcycling is enabled.
     * \param n_cycle      In/out: number of subcycles per level.
     * \param ref_ratio    Level-to-level refinement ratios.
     * \param dt_level     Output per-level dt.
     * \param stop_time    Simulation stop time (used to clamp dt if needed).
     */
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Vector<int>&           n_cycle,
                                   const Vector<IntVect>& ref_ratio,
                                   Vector<Real>&          dt_level,
                                   Real                  stop_time) = 0;
    /**
     * \brief Compute the next time step sizes after an advance.
     *
     * \param finest_level     Finest AMR level existing.
     * \param sub_cycle        Whether subcycling is enabled.
     * \param n_cycle          In/out: number of subcycles per level.
     * \param ref_ratio        Level-to-level refinement ratios.
     * \param dt_min           Output per-level dt upper bounds (physics-limited).
     * \param dt_level         Output per-level dt to use.
     * \param stop_time        Simulation stop time (used to clamp dt).
     * \param post_regrid_flag Flag describing regrid state (0 normal, 1 immediate post-regrid).
     */
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_min,
                               Vector<Real>&          dt_level,
                               Real                  stop_time,
                               int                   post_regrid_flag) = 0;
    /**
     * \brief Advance solution on this level by \p dt and return the next stable dt.
     *
     * \param time      Start time of the advance.
     * \param dt        Step size to take.
     * \param iteration Sub-step counter within the coarse step.
     * \param ncycle    Number of fine subcycles per coarse step.
     * \return Recommended dt for the next call (based on CFL, etc.).
     */
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) = 0;

    /**
     * \brief Perform level-specific updates after each subcycle.
     *
     * \param iteration Sub-step index within the coarse step.
     *
     * \note Override implementations should reset any FillPatcher caches they own.
     */
    virtual void post_timestep (int iteration);
    /**
     * \brief Called once per coarse step after fine levels finish refluxing.
     *
     * \param time Current coarse-level time.
     */
    virtual void postCoarseTimeStep (Real time);
    //! Override to perform any fixups immediately after restart (default: no-op).
    virtual void post_restart () {}
    /**
     * \brief Operations to be performed immediately after regridding completes.
     *
     * \param lbase      Coarsest level that triggered the regrid.
     * \param new_finest New finest level index after regridding.
     */
    virtual  void post_regrid (int lbase,
                               int new_finest) = 0;
    /**
     * \brief Operations to run after initialization finishes.
     *
     * \param stop_time Simulation stop time (useful for multi-stage inits).
     */
    virtual  void post_init (Real stop_time) = 0;
    /**
     * \brief Return nonzero to continue advancing; zero to halt the run.
     */
    virtual  int okToContinue () { return 1; }
    /**
     * \brief Return nonzero if refinement starting from this level is allowed this cycle.
     *
     * Only evaluated when `regrid_int > 0` and `level_count >= regrid_int`. Default returns true.
     */
    virtual  int okToRegrid ();
    /**
    * \brief Init grid data at problem start-up.
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual void initData () = 0;
    /**
     * \brief Set the time metadata stored in the state data arrays (old/new).
     *
     * \param time    New-time value to record.
     * \param dt_old  Dt associated with the old slot.
     * \param dt_new  Dt associated with the new slot.
     */
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new);
    //! Allocate storage for "old" time slices before an advance.
    virtual void allocOldData ();
    //! Delete the old-time slots once they are no longer needed.
    virtual void removeOldData ();
    /**
    * \brief Init data on this level from another AmrLevel (during regrid).
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    *
    * \param old Source level to interpolate from.
    */
    virtual void init (AmrLevel &old) = 0;
    /**
    * Init data on this level after regridding if old AmrLevel
    * did not previously exist. This is a pure virtual function
    * and hence MUST be implemented by derived classes.
    */
    virtual void init () = 0;
    //! Reset data to initial time by swapping new and old time slots.
    void reset ();
    //! Returns this AmrLevel index (0 = coarsest).
    int Level () const noexcept { return level; }
    //! List of grids at this level.
    const BoxArray& boxArray () const noexcept { return grids; }
    //! Return the edge-centered grid layout for direction \p dir (cached).
    const BoxArray& getEdgeBoxArray (int dir) const noexcept;
    //! Return the nodal grid layout for this level (cached).
    const BoxArray& getNodalBoxArray () const noexcept;
    //
    const DistributionMapping& DistributionMap () const noexcept { return dmap; }
    //
    const FabFactory<FArrayBox>& Factory () const noexcept { return *m_factory; }
    //
#ifdef AMREX_USE_EB
    const EBFArrayBoxFactory&
    EBFactory () const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*m_factory);
    }
#endif
    //
    //! Number of grids at this level.
    int numGrids () const noexcept { return static_cast<int>(grids.size()); }
    //! Number of states at this level.
    int numStates () const noexcept { return static_cast<int>(state.size()); }
    //! Returns the indices defining physical domain.
    const Box& Domain () const noexcept { return geom.Domain(); }
    //! Timestep n at this level.
    int nStep () const noexcept { return parent->levelSteps(level); }
    //! Returns the geometry object.
    const Geometry& Geom () const noexcept { return geom; }
    //
    const IntVect& fineRatio () const noexcept { return fine_ratio; }
    //! Returns number of valid cells on this level.
    Long countCells () const noexcept;

    //! Get the cached BoxArray describing regions that should not be tagged.
    const BoxArray& getAreaNotToTag () noexcept;
    //! Get the Box (single region) that is permitted to be tagged.
    const Box& getAreaToTag () noexcept;
    //! Build internal data structures defining the area not to tag.
    void constructAreaNotToTag ();
    //! Override the area-not-to-tag cache with user-supplied BoxArray \p ba.
    void setAreaNotToTag (BoxArray& ba) noexcept;

    //! Clear cached FillPatcher objects (call after reflux/averaging).
    void resetFillPatcher ();

    /**
    * \brief Error estimation for regridding. This is a pure virtual
    * function and hence MUST be implemented by derived classes.
    *
    * \param tb          Tag array to mark for refinement (output).
    * \param clearval    Value that clears a tag entry.
    * \param tagval      Value that marks a cell for refinement.
    * \param time        Physical time at which tagging occurs.
    * \param n_error_buf Number of grow cells added around tagged features (default 0).
    * \param ngrow       Number of ghost layers available in \p tb (default 0).
    */
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
                           int          n_error_buf = 0,
                           int          ngrow = 0) = 0;
    //! Interpolate from coarse level to the valid area in mf.
    /**
     * \brief Fill \p mf on this level using data interpolated from the coarser level.
     *
     * \param mf        Destination MultiFab.
     * \param dcomp     Destination component offset in \p mf.
     * \param time      Time at which to sample coarse data.
     * \param state_idx State descriptor index selecting which state to fill.
     * \param scomp     Source component index within the coarse state.
     * \param ncomp     Number of components to fill.
     * \param nghost    Number of ghost cells to fill (default 0).
     */
    void FillCoarsePatch (MultiFab& mf,
                          int       dcomp,
                          Real      time,
                          int       state_idx,
                          int       scomp,
                          int       ncomp,
                          int       nghost = 0);
    /**
     * \brief Apply physical boundary conditions to cell-centered data.
     *
     * \param dest        Destination FArrayBox.
     * \param state_indx  State descriptor index associated with \p dest.
     * \param time        Physical time at which BCs should be evaluated.
     * \param dest_comp   Destination component offset in \p dest.
     * \param src_comp    Component offset in the state data.
     * \param num_comp    Number of components to update.
     */
    virtual void setPhysBoundaryValues (FArrayBox& dest,
                                        int        state_indx,
                                        Real       time,
                                        int        dest_comp,
                                        int        src_comp,
                                        int        num_comp);
    /**
    * \brief Returns a MultiFab containing the derived data for this level.
    * If ngrow>0 the MultiFab is built on the appropriately
    * grown BoxArray.
    *
    * \param name  Derived quantity identifier.
    * \param time  Physical time at which to evaluate the quantity.
    * \param ngrow How many ghost cells to include.
    * \return Owning pointer to a freshly allocated MultiFab.
    */
    virtual std::unique_ptr<MultiFab> derive (const std::string& name,
                                              Real               time,
                                              int                ngrow);
    /**
    * \brief This version of derive() fills the dcomp'th component of \p mf
    * with the derived quantity.
    *
    * \param name Derived quantity identifier.
    * \param time Physical time at which to evaluate the quantity.
    * \param mf   Destination MultiFab.
    * \param dcomp Destination component index.
    */
    virtual void derive (const std::string& name,
                         Real               time,
                         MultiFab&          mf,
                         int                dcomp);
    //! State data object accessor for descriptor \p state_indx.
    StateData& get_state_data (int state_indx) noexcept { return state[state_indx]; }
    //! Mutable access to the MultiFab at the old time slot for state \p state_indx.
    MultiFab& get_old_data (int state_indx) noexcept { return state[state_indx].oldData(); }
    //! Const access to the MultiFab at the old time slot for state \p state_indx.
    const MultiFab& get_old_data (int state_indx) const noexcept { return state[state_indx].oldData(); }
    //! Mutable access to the MultiFab at the new time slot for state \p state_indx.
    MultiFab& get_new_data (int state_indx) noexcept { return state[state_indx].newData(); }
    //! Const access to the MultiFab at the new time slot for state \p state_indx.
    const MultiFab& get_new_data (int state_indx) const noexcept { return state[state_indx].newData(); }
    //! Returns list of Descriptors.
    static const DescriptorList& get_desc_lst () noexcept { return desc_lst; }
    //! Returns list of derived variables.
    static DeriveList& get_derive_lst () noexcept;
    //! Returns whether a post-timestep regrid is currently requested (nonzero for true).
    int postStepRegrid () const noexcept { return post_step_regrid; }
    /**
     * \brief Set the post-timestep regrid trigger (nonzero requests a regrid).
     *
     * \param new_val Trigger value applied after each step.
     */
    void setPostStepRegrid (int new_val) noexcept { post_step_regrid = new_val; }

    /**
     * \brief Update every StateData's DistributionMapping to match \p dmap (e.g., after rebalancing).
     *
     * \param dmap Distribution mapping describing the new ownership.
     */
    void UpdateDistributionMaps ( DistributionMapping& dmap );

    /**
     * \brief Return an array describing boundary types for components [strt_comp,strt_comp+ncomp).
     *
     * \param State_Type State descriptor index.
     * \param gridno     Grid index on this level.
     * \param strt_comp  Starting component.
     * \param ncomp      Number of components.
     */
    Vector<int> getBCArray (int State_Type,
                            int gridno,
                            int strt_comp,
                            int ncomp);
    /**
     * \brief Fetch the MultiFab (old or new) that best matches \p time for state \p state_indx.
     *
     * \param state_indx State descriptor index.
     * \param time       Target physical time.
     * \return Reference to either the old or new MultiFab, whichever brackets \p time.
     */
    MultiFab& get_data (int  state_indx, Real time) noexcept;
    /**
     * \brief Hook to override fillpatched boundary values before applying physical BCs (default no-op).
     *
     * \param S            MultiFab to edit.
     * \param state_index  State descriptor index.
     * \param scomp        Source component offset.
     * \param dcomp        Destination component offset.
     * \param ncomp        Number of components to adjust.
     * \param time         Physical time at which BCs are evaluated.
     */
    virtual void set_preferred_boundary_values (MultiFab& S,
                                                int       state_index,
                                                int       scomp,
                                                int       dcomp,
                                                int       ncomp,
                                                Real      time) const;
    /**
     * \brief Adjust tags after core tagging routines run (default: no-op).
     *
     * \param tags   Tag arrays built on this level.
     * \param bf_lev Buffer sizes per level from Amr.
     */
    virtual void manual_tags_placement (TagBoxArray&    tags,
                                        const Vector<IntVect>& bf_lev);
    //! Override to customize which variables appear in full plotfiles on this level (default uses global list).
    virtual void setPlotVariables ();
    //! Override to customize which variables appear in small plotfiles on this level.
    virtual void setSmallPlotVariables ();
    /**
     * \brief Estimate the amount of work required to advance just this level, typically based on cell counts.
     *
     * \return Work units proportional to the expected compute effort (higher means more expensive).
     */
    virtual Real estimateWork();

    /**
     * \brief Select which state descriptor drives work estimates.
     *
     * \return State descriptor index used for estimateWork() sampling, or -1 to disable.
     */
    virtual int WorkEstType () { return -1; }

    /**
    * \brief Returns one the TimeLevel enums.
    * Asserts that time is between AmrOldTime and AmrNewTime.
    *
    * \param state_indx State descriptor index providing time metadata.
    * \param time       Target physical time.
    */
    TimeLevel which_time (int  state_indx, Real time) const noexcept;

    /**
     * \brief Return true if this level requests an immediate checkpoint (default false).
     *
     * Overrides can trigger checkpointing after conditional events.
     */
    virtual bool checkPointNow ();

    /**
     * \brief Return true if this level requests an immediate plotfile (default false).
     */
    virtual bool writePlotNow ();

    /**
     * \brief Return true if this level requests an immediate small-plot file (default false).
     */
    virtual bool writeSmallPlotNow ();

#ifdef AMREX_PARTICLES
    /**
     * \brief Optional particle redistribution hook called by the parent Amr object.
     *
     * \param lbase Coarsest level whose particles should be considered.
     * \param a_init Whether the call happens during initialization.
     */
    virtual void particle_redistribute (int lbase = 0, bool a_init = false)
    {
        (void)lbase;
        (void)a_init;
    }
#endif

    /**
     * \brief Fill ghost cells using cached FillPatcher (level>0) or FillPatch (level 0).
     *
     * \param mf          Destination MultiFab.
     * \param dcomp       Destination component offset.
     * \param ncomp       Number of components to fill.
     * \param nghost      Number of ghost cells to fill.
     * \param time        Time at which to sample data.
     * \param state_index State descriptor index to pull from.
     * \param scomp       Source component offset within the state.
     */
    void FillPatcherFill (amrex::MultiFab& mf, int dcomp, int ncomp, int nghost,
                          amrex::Real time, int state_index, int scomp);

    /**
     * \brief Convenience wrapper around FillPatchSingleLevel/TwoLevels.
     *
     * \param amrlevel  Level providing state data.
     * \param leveldata Destination MultiFab to fill.
     * \param boxGrow   Number of grow cells to request.
     * \param time      Time at which to sample.
     * \param index     State descriptor index.
     * \param scomp     Source component.
     * \param ncomp     Number of components.
     * \param dcomp     Destination component (default 0).
     */
    static void FillPatch (AmrLevel& amrlevel,
                           MultiFab& leveldata,
                           int       boxGrow,
                           Real      time,
                           int       index,
                           int       scomp,
                           int       ncomp,
                           int       dcomp=0);

    /**
     * \brief Same as FillPatch(), but adds into \p leveldata instead of overwriting.
     *
     * \param amrlevel  Level providing state data.
     * \param leveldata Destination MultiFab to accumulate into.
     * \param boxGrow   Number of grow cells to request.
     * \param time      Time at which to sample.
     * \param index     State descriptor index.
     * \param scomp     Source component.
     * \param ncomp     Number of components.
     * \param dcomp     Destination component (default 0).
     */
    static void FillPatchAdd (AmrLevel& amrlevel,
                              MultiFab& leveldata,
                              int       boxGrow,
                              Real      time,
                              int       index,
                              int       scomp,
                              int       ncomp,
                              int       dcomp=0);

    /**
     * \brief Evolve one step with Runge-Kutta (2, 3, or 4)
     *
     * To use RK, the StateData must have all the ghost cells needed.  See
     * namespace RungeKutta for expected function signatures of the callable
     * parameters.
     *
     * \param order      order of RK
     * \param state_type index of StateData
     * \param time       time at the beginning of the step.
     * \param dt         time step
     * \param iteration  iteration number on fine level during a coarse time
     *                   step.  For an AMR simulation with subcycling and a
     *                   refinement ratio of 2, the number is either 1 or 2,
     *                   denoting the first and second substep, respectively.
     * \param ncycle     number of subcyling steps.  It's usually 2 or 4.
     *                   Without subcycling, this will be 1.
     * \param f          computing right-hand side for evolving the StateData.
     *                   One can also register data for flux registers in this.
     * \param p          optionally post-processing RK stage results
     */
    template <typename F, typename P = RungeKutta::PostStageNoOp>
    void RK (int order, int state_type, Real time, Real dt, int iteration,
             int ncycle, F&& f, P&& p = RungeKutta::PostStageNoOp());

#ifdef AMREX_USE_EB
    //! Configure the ghost-cell requirements (basic/volume/full) for EB FillPatch operations.
    static void SetEBMaxGrowCells (int nbasic, int nvolume, int nfull) noexcept {
        m_eb_basic_grow_cells = nbasic;
        m_eb_volume_grow_cells = nvolume;
        m_eb_full_grow_cells = nfull;
    }
    static int            m_eb_basic_grow_cells;
    static int            m_eb_volume_grow_cells;
    static int            m_eb_full_grow_cells;
    //! Select which level of EB data to build (basic, volume, full) via \p ebs.
    static void SetEBSupportLevel (EBSupport ebs) { m_eb_support_level = ebs; }
    static EBSupport      m_eb_support_level;
#endif

    /**
     * \brief Recommendation of a proper blocking factor for filling.
     *
     * \param amr_level Owning level supplying geometry.
     * \param boxGrow   Grow cells requested by the caller.
     * \param boxType   Index type (cell-centered, nodal, etc.).
     * \param desc      State descriptor tied to the data.
     * \param SComp     Starting component index.
     */
    static IntVect ProperBlockingFactor (AmrLevel const& amr_level, int boxGrow,
                                         IndexType const& boxType,
                                         StateDescriptor const& desc, int SComp);

protected:
    //! The constructors -- for derived classes.
    AmrLevel () noexcept {} // NOLINT

    AmrLevel (Amr&            papa,
              int             lev,
              const Geometry& level_geom,
              const BoxArray& ba,
              const DistributionMapping& dm,
              Real            time);

    //! Common code used by all constructors.
    void finishConstructor ();

    //
    // The Data.
    //
    int                   level{-1};    // AMR level (0 is coarsest).
    Geometry              geom;         // Geom at this level.
    BoxArray              grids;        // Cell-centered locations of grids.
    DistributionMapping   dmap;         // Distribution of grids among processes
    Amr*                  parent{nullptr};// Pointer to parent AMR structure.
    IntVect               crse_ratio;   // Refinement ratio to coarser level.
    IntVect               fine_ratio;   // Refinement ratio to finer level.
    static DeriveList     derive_lst;   // List of derived quantities.
    static DescriptorList desc_lst;     // List of state variables.
    Vector<StateData>     state;        // Array of state data.

    BoxArray m_AreaNotToTag; //Area which shouldn't be tagged on this level.
    Box      m_AreaToTag;    //Area which is allowed to be tagged on this level.

    int post_step_regrid{0}; // Whether or not to do a regrid after the timestep.

    bool levelDirectoryCreated{false}; // for checkpoints and plotfiles

    std::unique_ptr<FabFactory<FArrayBox> > m_factory;

    Vector<std::unique_ptr<FillPatcher<MultiFab>>> m_fillpatcher;

private:

    template <std::size_t order>
    void storeRKCoarseData (int state_type, Real time, Real dt,
                            MultiFab const& S_old,
                            Array<MultiFab,order> const& rkk);

    void FillRKPatch (int state_index, MultiFab& S, Real time,
                      int stage, int iteration, int ncycle);

    mutable BoxArray      edge_grids[AMREX_SPACEDIM];  // face-centered grids
    mutable BoxArray      nodal_grids;              // all nodal grids
};

//
// Forward declaration.
//
class FillPatchIteratorHelper;

class FillPatchIterator
    :
    public MFIter
{
  public:

    friend class AmrLevel;

    /**
     * \brief Construct an iterator that fills the default state (new-time data) from scratch.
     *
     * \param amrlevel Owning level that provides state data.
     * \param leveldata Destination MultiFab that will hold the filled tiles.
     */
    FillPatchIterator (AmrLevel& amrlevel,
                       MultiFab& leveldata);

    /**
     * \brief Construct and immediately fill ghost cells using specified state components.
     *
     * \param amrlevel Owning level.
     * \param leveldata Destination MultiFab.
     * \param boxGrow  Number of grow cells to request.
     * \param time     Time at which to sample data.
     * \param idx      State descriptor index.
     * \param scomp    Source component offset.
     * \param ncomp    Number of components.
     */
    FillPatchIterator (AmrLevel& amrlevel,
                       MultiFab& leveldata,
                       int       boxGrow,
                       Real      time,
                       int       idx,
                       int       scomp,
                       int       ncomp);

    /**
     * \brief Prepare another round of fill operations using the cached buffers.
     *
     * \param boxGrow Number of grow cells to request.
     * \param time    Time at which to sample state data.
     * \param idx     State descriptor index.
     * \param scomp   Source component offset.
     * \param ncomp   Number of components to fill.
     */
    void Initialize (int  boxGrow,
                     Real time,
                     int  idx,
                     int  scomp,
                     int  ncomp);

    //! Access the filled FAB for the current tile.
    FArrayBox& operator() () noexcept { return m_fabs[MFIter::index()]; }

    //! Return the un-grown valid box for the current tile.
    Box UngrownBox () const noexcept { return MFIter::validbox(); }

    //! Mutable reference to the internal MultiFab holding filled tiles.
    MultiFab& get_mf() noexcept { return m_fabs; }

private:

    void FillFromLevel0 (Real time, int idx, int scomp, int dcomp, int ncomp);
    void FillFromTwoLevels (Real time, int idx, int scomp, int dcomp, int ncomp);

    //
    // The data.
    //
    AmrLevel*                         m_amrlevel;
    MultiFab*                         m_leveldata;
    std::vector< std::pair<int,int> > m_range;
    MultiFab                          m_fabs;
    int                               m_ncomp;
};

class FillPatchIteratorHelper
{
public:

    friend class FillPatchIterator;

    //! Build a helper bound to \p amrlevel and destination MultiFab \p leveldata.
    FillPatchIteratorHelper (AmrLevel& amrlevel,
                             MultiFab& leveldata);

    /**
     * \brief Construct a helper and immediately populate ghost cells.
     *
     * \param amrlevel Source level providing state data.
     * \param leveldata Destination MultiFab.
     * \param boxGrow  Number of grow cells to request.
     * \param time     Time at which to sample data.
     * \param index    State descriptor index.
     * \param scomp    Source component offset.
     * \param ncomp    Number of components to fill.
     * \param mapper   Interpolater override (defaults to descriptor's mapper when null).
     */
    FillPatchIteratorHelper (AmrLevel&     amrlevel,
                             MultiFab&     leveldata,
                             int           boxGrow,
                             Real          time,
                             int           index,
                             int           scomp,
                             int           ncomp,
                             InterpBase*   mapper);

    /**
     * \brief Prepare another fill using cached metadata.
     *
     * \param boxGrow Number of grow cells to request.
     * \param time    Time at which to sample data.
     * \param idx     State descriptor index.
     * \param scomp   Source component offset.
     * \param ncomp   Number of components to fill.
     * \param mapper  Interpolater override (may be null).
     */
    void Initialize (int           boxGrow,
                     Real          time,
                     int           idx,
                     int           scomp,
                     int           ncomp,
                     InterpBase*   mapper);

    /**
     * \brief Fill the provided FAB using cached data.
     *
     * \param fab   Destination FAB.
     * \param dcomp Destination component offset.
     * \param idx   State descriptor index.
     */
    void fill (FArrayBox& fab, int dcomp, int idx);

private:
    //
    // The data.
    //
    AmrLevel*                  m_amrlevel;
    MultiFab*                  m_leveldata;
    MultiFabCopyDescriptor     m_mfcd;
    Vector< Vector<MultiFabId> > m_mfid;     // [level][oldnew]
    Interpolater*              m_map = nullptr;
    std::map<int,Box>          m_ba;
    Real                       m_time = std::numeric_limits<Real>::lowest();
    int                        m_growsize = std::numeric_limits<int>::lowest();
    int                        m_index = -1;
    int                        m_scomp = -1;
    int                        m_ncomp = -1;
    bool                       m_FixUpCorners = false;

    std::map< int,Vector< Vector<Box> > >                m_fbox; // [grid][level][validregion]
    std::map< int,Vector< Vector<Box> > >                m_cbox; // [grid][level][fillablesubbox]
    std::map< int,Vector< Vector< Vector<FillBoxId> > > > m_fbid; // [grid][level][fillablesubbox][oldnew]
};

template <typename F, typename P>
void AmrLevel::RK (int order, int state_type, Real time, Real dt, int iteration,
                   int ncycle, F&& f, P&& p)
{
    BL_PROFILE("AmrLevel::RK()");

    AMREX_ASSERT(AmrLevel::desc_lst[state_type].nExtra() > 0); // Need ghost cells in StateData

    MultiFab& S_old = get_old_data(state_type);
    MultiFab& S_new = get_new_data(state_type);

    if (order == 2) {
        RungeKutta::RK2(S_old, S_new, time, dt, std::forward<F>(f),
                        [&] (int /*stage*/, MultiFab& mf, Real t) {
                            FillPatcherFill(mf, 0, mf.nComp(), mf.nGrow(), t,
                                            state_type, 0); },
                        std::forward<P>(p));
    } else if (order == 3) {
        RungeKutta::RK3(S_old, S_new, time, dt, std::forward<F>(f),
                        [&] (int stage, MultiFab& mf, Real t) {
                            FillRKPatch(state_type, mf, t, stage, iteration, ncycle);
                        },
                        [&] (Array<MultiFab,3> const& rkk) {
                            if (level < parent->finestLevel()) {
                                storeRKCoarseData(state_type, time, dt, S_old, rkk);
                            }
                        },
                        std::forward<P>(p));
    } else if (order == 4) {
        RungeKutta::RK4(S_old, S_new, time, dt, std::forward<F>(f),
                        [&] (int stage, MultiFab& mf, Real t) {
                            FillRKPatch(state_type, mf, t, stage, iteration, ncycle);
                        },
                        [&] (Array<MultiFab,4> const& rkk) {
                            if (level < parent->finestLevel()) {
                                storeRKCoarseData(state_type, time, dt, S_old, rkk);
                            }
                        },
                        std::forward<P>(p));
    } else {
        amrex::Abort("AmrLevel::RK: order = "+std::to_string(order)+" is not supported");
    }
}

template <std::size_t order>
void AmrLevel::storeRKCoarseData (int state_type, Real time, Real dt,
                                  MultiFab const& S_old,
                                  Array<MultiFab,order> const& rkk)
{
    BL_PROFILE("AmrLevel::storeRKCoarseData()");
    if (level == parent->finestLevel()) { return; }

    const StateDescriptor& desc = AmrLevel::desc_lst[state_type];

    auto& fillpatcher = parent->getLevel(level+1).m_fillpatcher[state_type];
    fillpatcher = std::make_unique<FillPatcher<MultiFab>>
        (parent->boxArray(level+1), parent->DistributionMap(level+1),
         parent->Geom(level+1),
         parent->boxArray(level), parent->DistributionMap(level),
         parent->Geom(level),
         IntVect(desc.nExtra()), desc.nComp(), desc.interp(0));

    fillpatcher->storeRKCoarseData(time, dt, S_old, rkk);
}


}

#endif /*_AmrLevel_H_*/
