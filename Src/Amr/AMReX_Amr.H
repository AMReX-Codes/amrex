
#ifndef AMREX_Amr_H_
#define AMREX_Amr_H_

#include <fstream>
#include <memory>
#include <list>

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BCRec.H>

#include <AMReX_AmrCore.H>

#ifdef USE_PERILLA
#include <RegionGraph.H>
#include <Perilla.H>
#endif

namespace amrex {

class AmrLevel;
class LevelBld;
class BoxDomain;
template <class T>
class MFGraph;
class AmrTask;
#if defined(BL_USE_SENSEI_INSITU)
class AmrInSituBridge;
#endif

/**
* \brief Manage hierarchy of levels for time-dependent AMR computations.
*
* The Amr class is designed to manage parts of the computation  which do
* not belong on a single level, like establishing and updating the hierarchy
* of levels, global timestepping, and managing the different AmrLevels
*/

class Amr
    : public AmrCore
{
  template <class T>
  friend class MFGraph;
  friend class AmrTask;
  typedef std::multimap< std::pair<int, int>, double >  BoundaryPointList;

public:
    //! The constructor.
    Amr ();

    Amr (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord);

    Amr (const Amr& rhs) = delete;
    Amr& operator= (const Amr& rhs) = delete;

    void InitAmr ();

    //! The destructor.
    virtual ~Amr ();

    //! Init data after construction. Must be called before timestepping.
    virtual void init (Real strt_time, Real stop_time);

    //! First part of initialInit
    void InitializeInit (Real strt_time, Real stop_time,
                         const BoxArray* lev0_grids = 0, const Vector<int>* pmap = 0);

    //! Second part of initialInit
    void FinalizeInit (Real strt_time, Real stop_time);

    //! Set the timestep on each level.
    void setDtLevel (const Vector<Real>& dt_lev) noexcept;

    //! Set the timestep at one level.
    void setDtLevel (Real dt, int lev) noexcept;

    //! Set the dtmin on each level.
    void setDtMin (const Vector<Real>& dt_lev) noexcept;

    //! Set the cycle count on each level.
    void setNCycle (const Vector<int>& mss) noexcept;

    //! Subcycle in time?
    int subCycle () const noexcept { return sub_cycle; }

    //! How are we subcycling?
    const std::string& subcyclingMode() const noexcept { return subcycling_mode; }

    /**
    * \brief What is "level" in Amr::timeStep?  This is only relevant if we are still in Amr::timeStep;
    *      it is set back to -1 on leaving Amr::timeStep.
    */
    int level_being_advanced () const noexcept { return which_level_being_advanced; }
    //! Physical time.
    Real cumTime () const noexcept { return cumtime; }
    void setCumTime (Real t) noexcept {cumtime = t;}
    //! Physical time this simulation started
    Real startTime () const noexcept { return start_time; }
    void setStartTime (Real t) noexcept {start_time = t;}
    //! Time step at specified level.
    Real dtLevel (int level) const noexcept { return dt_level[level]; }
    //! Max time step (typically based on physics) at specified level
    Real dtMin (int level) const noexcept { return dt_min[level]; }
    //! Array of time steps at all levels.
    const Vector<Real>& dtLevel () const noexcept { return dt_level; }
    //! Number of subcycled time steps.
    int nCycle (int level) const noexcept { return n_cycle[level]; }
    //! Number of time steps at specified level.
    int levelSteps (int lev) const noexcept { return level_steps[lev]; }
    //! Number of time steps at specified level.
    void setLevelSteps (int lev, int n) noexcept { level_steps[lev] = n; }
    //! Which step are we at for the specified level?
    int levelCount (int lev) const noexcept { return level_count[lev]; }
    //! Which step are we at for the specified level?
    void setLevelCount (int lev, int n) noexcept { level_count[lev] = n; }
    //! Whether to regrid right after restart
    bool RegridOnRestart () const noexcept;
    //! Interval between regridding.
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
    static const std::list<std::string>& stateSmallPlotVars () noexcept { return state_small_plot_vars; }
    //! Is the string the name of a variable in state_plot_vars?
    static bool isStatePlotVar (const std::string& name);
    static bool isStateSmallPlotVar (const std::string& name);
    /**
    * \brief If the string is not the name of a variable in state_plot_vars,
    * add it to state_plot_vars.
    */
    static void addStatePlotVar (const std::string& name);
    static void addStateSmallPlotVar (const std::string& name);
    //! Remove the string from state_plot_vars.
    static void deleteStatePlotVar (const std::string& name);
    //! Clear the list of state_plot_vars.
    static void clearStatePlotVarList ();
    static void clearStateSmallPlotVarList ();
    //!  Fill the list of state_plot_vars with all of the state quantities.
    static void fillStatePlotVarList ();
    static void fillStateSmallPlotVarList ();
    //!  Write out plotfiles (True/False)?
    static bool Plot_Files_Output ();
    /**
    * \brief The names of derived variables to output in the
    * plotfile.  They can be set using the amr.derive_plot_vars
    * variable in a ParmParse inputs file.
    */
    static const std::list<std::string>& derivePlotVars () noexcept { return derive_plot_vars; }
    static const std::list<std::string>& deriveSmallPlotVars () noexcept { return derive_small_plot_vars; }
    //! Is the string the name of a variable in derive_plot_vars?
    static bool isDerivePlotVar (const std::string& name) noexcept;
    static bool isDeriveSmallPlotVar (const std::string& name) noexcept;
    /**
    * \brief If the string is not the name of a variable in
    * derive_plot_vars, add it to derive_plot_vars.
    */
    static void addDerivePlotVar (const std::string& name);
    static void addDeriveSmallPlotVar (const std::string& name);
    //! Remove the string from derive_plot_vars.
    static void deleteDerivePlotVar (const std::string& name);
    static void deleteDeriveSmallPlotVar (const std::string& name);
    //! Clear the list of derive_plot_vars.
    static void clearDerivePlotVarList ();
    static void clearDeriveSmallPlotVarList ();
    //!  Fill the list of derive_plot_vars with all derived quantities.
    static void fillDerivePlotVarList ();
    static void fillDeriveSmallPlotVarList ();

    static void Initialize ();
    static void Finalize ();
    //! AmrLevel lev.
    AmrLevel& getLevel (int lev) noexcept { return *amr_level[lev]; }
    //! Array of AmrLevels.
    Vector<std::unique_ptr<AmrLevel> >& getAmrLevels () noexcept;
    //! Total number of cells.
    long cellCount () noexcept;
    //! Number of cells at given level.
    long cellCount (int lev) noexcept;
    //! Total number of grids.
    int numGrids () noexcept;
    //! Number of grids at given level.
    int numGrids (int lev) noexcept;
    //! More work to be done?
    int okToContinue () noexcept;
    //! Regrid only!
    void RegridOnly (Real time, bool do_io = true);
    //! Should we regrid this level?
    bool okToRegrid (int level) noexcept;
    //! Array of BoxArrays read in to initially define grid hierarchy
    static const BoxArray& initialBa (int level) noexcept
        { BL_ASSERT(level-1 < initial_ba.size()); return initial_ba[level-1]; }
    //! Number of levels at which the grids are initially specified
    static int initialBaLevels () noexcept { return initial_ba.size(); }
    //! Do a complete integration cycle.
    virtual void coarseTimeStep (Real stop_time);

    //! Do a complete integration cycle and return the coarse dt.
    Real coarseTimeStepDt (Real stop_time);
    //! Retrieve derived data. User is responsible for deleting pointer.
    std::unique_ptr<MultiFab> derive (const std::string& name,
				      Real           time,
				      int            lev,
				      int            ngrow);
    //! Name of the restart chkpoint file.
    const std::string& theRestartFile () const noexcept { return restart_chkfile; }
    //! Name of the restart plotfile.
    const std::string& theRestartPlotFile () const noexcept { return restart_pltfile; }
    //! The ith datalog file.  Do with it what you want.
    std::ostream& DataLog (int i);
    //! The filename of the ith datalog file.
    const std::string DataLogName (int i) const noexcept { return datalogname[i]; }
    //! How many datalogs have been opened
    int NumDataLogs () noexcept;
    /**
    * \brief Compute the optimal subcycling pattern.
    * This assumes that anything less than cycle_max[i] is a valid
    * number of subcycles at level[i]. For example
    * if ref_ratio[i] = cycle_max[i] = 4, then 1,2,3,4 are all valid
    * values for n_cycles[i]
    */
    static Real computeOptimalSubcycling (int   n,
                                          int*  best,
                                          Real* dt_max,
                                          Real* est_work,
                                          int*  cycle_max);

    //! Write the plot file to be used for visualization.
    virtual void writePlotFile ();
    int stepOfLastPlotFile () const noexcept {return last_plotfile;}
    //! Write the small plot file to be used for visualization.
    virtual void writeSmallPlotFile ();
    int stepOfLastSmallPlotFile () const noexcept {return last_smallplotfile;}
    //! Write current state into a chk* file.
    virtual void checkPoint ();
    int stepOfLastCheckPoint () const noexcept {return last_checkpoint;}

    const Vector<BoxArray>& getInitialBA() noexcept;

    /**
    * \brief Specialized version:
    * Define BoundaryPointLists that give the intersections
    *    of the external geometry with constant (i,k) and (j,k)
    * These are defined at the coarsest level indexing only.
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
    };

    /**
    * \brief More general version:
    * Define BoundaryPointLists that give the intersections
    *    of the external geometry with constant (i,k),(j,k)
    *    and (i,j).
    * These are defined at the coarsest level indexing only.
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
    };

    BoundaryPointList& getIntersectLoX() noexcept
    {
        return intersect_lox;
    };
    BoundaryPointList& getIntersectHiX() noexcept
    {
        return intersect_hix;
    };
    BoundaryPointList& getIntersectLoY() noexcept
    {
        return intersect_loy;
    };
    BoundaryPointList& getIntersectHiY() noexcept
    {
        return intersect_hiy;
    };
    BoundaryPointList& getIntersectLoZ() noexcept
    {
        return intersect_loz;
    };
    BoundaryPointList& getIntersectHiZ() noexcept
    {
        return intersect_hiz;
    };

#ifdef AMREX_PARTICLES
    //! Redistribute particles
    void RedistributeParticles ();
#endif

    void InstallNewDistributionMap (int lev, const DistributionMapping& newdm);

    bool UsingPrecreateDirectories () noexcept;

protected:

    //! Initialize grid hierarchy -- called by Amr::init.
    void initialInit (Real strt_time, Real stop_time,
                      const BoxArray* lev0_grids = 0, const Vector<int>* pmap = 0);
    //! Read the probin file.
    void readProbinFile (int& init);
    //! Check for valid input.
    void checkInput ();
    //! Restart from a checkpoint file.
    void restart (const std::string& filename);
    //! Define and initialize coarsest level.
    void defBaseLevel (Real start_time, const BoxArray* lev0_grids = 0, const Vector<int>* pmap = 0);
    //! Define and initialize refined levels.
    void bldFineLevels (Real start_time);
    //! Rebuild grid hierarchy finer than lbase.
    virtual void regrid (int  lbase,
                         Real time,
                         bool initial = false) override;
    //! Regrid level 0 on restart.
    virtual void regrid_level_0_on_restart ();
    //! Define new grid locations (called from regrid) and put into new_grids.
    void grid_places (int              lbase,
                      Real             time,
                      int&             new_finest,
                      Vector<BoxArray>& new_grids);

    DistributionMapping makeLoadBalanceDistributionMap (int lev, Real time, const BoxArray& ba) const;
    void LoadBalanceLevel0 (Real time);

    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override;
    virtual BoxArray GetAreaNotToTag (int lev) override;
    virtual void ManualTagsPlacement (int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev) override;

    //! Do a single timestep on level L.
    virtual void timeStep (int  level,
                           Real time,
                           int  iteration,
                           int  niter,
                           Real stop_time);

    // pure virtural function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override
	{ amrex::Abort("How did we get her!"); }
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override
	{ amrex::Abort("How did we get her!"); }
    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) override
	{ amrex::Abort("How did we get her!"); }
    virtual void ClearLevel (int lev) override
	{ amrex::Abort("How did we get her!"); }

    //! Whether to write a plotfile now
    bool writePlotNow () noexcept;
    bool writeSmallPlotNow () noexcept;

    void printGridInfo (std::ostream& os,
                        int           min_lev,
                        int           max_lev);

    void setRecordGridInfo (const std::string&);

    void setRecordRunInfo (const std::string&);

    void setRecordRunInfoTerse (const std::string&);

    void setRecordDataInfo (int i, const std::string&);

    void initSubcycle();
    void initPltAndChk();

    int initInSitu();
    int updateInSitu();
    static int finalizeInSitu();

    //
    // The data ...
    //
    std::string      regrid_grids_file;   //!< Grids file that will bypass regridding.
    std::string      initial_grids_file;  //!< Grids file that will bypass regridding only at initialization.
    Vector<std::unique_ptr<AmrLevel> > amr_level;    //!< Vector of levels
    Real             cumtime;      //!< Physical time variable.
    Real             start_time;   //!< Physical time this simulation started.
    Vector<Real>      dt_level;     //!< Timestep at this level.
    Vector<int>       level_steps;  //!< Number of time steps at this level.
    Vector<int>       level_count;
    Vector<int>       n_cycle;
    std::string      subcycling_mode; //!<Type of subcycling to use.
    Vector<Real>      dt_min;
    bool             isPeriodic[AMREX_SPACEDIM];  //!< Domain periodic?
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
    int              small_plot_int;  //!< How often small plotfile (# of time steps)
    Real             small_plot_per;  //!< How often small plotfile (in units of time)
    Real             small_plot_log_per;  //!< How often small plotfile (in units of log10(time))
    int              write_plotfile_with_checkpoint;  //!< Write out a plotfile whenever we checkpoint
    int              file_name_digits; //!< How many digits to use in the plotfile and checkpoint names
    int              message_int;     //!< How often checking messages touched by user, such as "stop_run"
    std::string      plot_file_root;  //!< Root name of plotfile.
    std::string      small_plot_file_root;  //!< Root name of small plotfile.

    int              which_level_being_advanced; //!< Only >=0 if we are in Amr::timeStep(level,...)

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
    std::string      probin_file;
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

#if defined(BL_USE_SENSEI_INSITU)
    static AmrInSituBridge *insitu_bridge;
#endif

public:
#ifdef USE_PERILLA
    std::vector<std::vector<RegionGraph*> > graphArray;
    std::vector<RegionGraph*> amrGraphArray;
    std::vector<RegionGraph*> &get_graphArray(int level){return graphArray[level];}
#endif

    BoundaryPointList intersect_lox;
    BoundaryPointList intersect_loy;
    BoundaryPointList intersect_loz;
    BoundaryPointList intersect_hix;
    BoundaryPointList intersect_hiy;
    BoundaryPointList intersect_hiz;

    static bool first_smallplotfile;

};

}

#endif /*_Amr_H_*/
