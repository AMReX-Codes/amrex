
#ifndef _Amr_H_
#define _Amr_H_ 

#include <fstream>

#include <Box.H>
#include <Geometry.H>
#include <BoxArray.H>
#include <Array.H>
#include <PArray.H>
#include <BCRec.H>

#ifdef USE_STATIONDATA
#include <StationData.H>
#endif

#ifdef USE_PARTICLES
#include <ParGDB.H>
#endif

class AmrLevel;
class LevelBld;
class BoxDomain;

#ifdef USE_PARTICLES
class AmrParGDB;
#endif

//
//
// Manages hierarchy of levels for time-dependent AMR computations.
//
// The Amr class is designed to manage parts of the computation  which do
// not belong on a single level, like establishing and updating the hierarchy
// of levels, global timestepping, and managing the different AmrLevels
//

class Amr
{
    typedef std::multimap< std::pair<int, int>, double >  BoundaryPointList;

public:
    //
    // The constructor.
    //
    Amr ();
    Amr (const RealBox* rb, int max_level_in, Array<int> n_cell_in, int coord);
    void InitAmr(int max_level_in, Array<int> n_cell_in);
    //
    // The destructor.
    //
    virtual ~Amr ();
    //
    // Init data after construction. Must be called before timestepping.
    //
    virtual void init (Real strt_time, Real stop_time);
    //
    // First part of initialInit
    //
    void InitializeInit (Real strt_time, Real stop_time, 
                         const BoxArray* lev0_grids = 0, const Array<int>* pmap = 0);
    //
    // Second part of initialInit
    //
    void FinalizeInit (Real strt_time, Real stop_time); 
    //
    // Set the timestep on each level.
    //
    void setDtLevel (const Array<Real>& dt_lev);
    //
    // Set the timestep at one level.
    //
    void setDtLevel (Real dt, int lev);
    //
    // Set the dtmin on each level.
    //
    void setDtMin (const Array<Real>& dt_lev);
    //
    // Set the cycle count on each level.
    //
    void setNCycle (const Array<int>& mss);
    //
    // return blocking_factor
    //
    int blockingFactor (int lev) const { return blocking_factor[lev]; }
    //
    // Largest allowable grid.
    //
    int maxGridSize (int lev) const { return max_grid_size[lev]; }
    //
    // Subcycle in time?
    //
    int subCycle () const { return sub_cycle; }
    //
    // How are we subcycling?
    //
    const std::string& subcyclingMode() const { return subcycling_mode; }
    //
    // Maximum allowable level.
    //
    int maxLevel () const { return max_level; }
    //
    // Finest currently existing level.
    //
    int finestLevel () const { return finest_level; }
    //  
    // What is "level" in Amr::timeStep?  This is only relevant if we are still in Amr::timeStep;
    //      it is set back to -1 on leaving Amr::timeStep.
    //  
    int level_being_advanced () const { return which_level_being_advanced; }
    //
    // Physical time.
    //
    Real cumTime () const { return cumtime; }
    void setCumTime (Real t) {cumtime = t;}
    //
    // Physical time this simulation started
    //
    Real startTime () const { return start_time; }
    void setStartTime (Real t) {start_time = t;}
    //
    // Time step at specified level.
    //
    Real dtLevel (int level) const { return dt_level[level]; }
    //
    // Max time step (typically based on physics) at specified level
    //
    Real dtMin (int level) const { return dt_min[level]; }
    //
    // Array of time steps at all levels.
    //
    const Array<Real>& dtLevel () const { return dt_level; }
    //
    // Geometry object.
    //
    const Geometry& Geom (int level) const { return geom[level]; }
    //
    // Number of subcycled time steps.
    //
    int nCycle (int level) const { return n_cycle[level]; }
    //
    // Refinement ratio.
    //
    IntVect refRatio (int _level) const { return ref_ratio[_level]; }
    //
    // Maximum refinement ratio in any direction.
    //
    int MaxRefRatio (int level) const;
    //
    // Refinement ratios between all levels.
    //
    const Array<IntVect>& refRatio () const { return ref_ratio; }
    //
    // Number of time steps at specified level.
    //
    int levelSteps (int lev) const { return level_steps[lev]; }
    //
    // Number of time steps at specified level.
    //
    void setLevelSteps (int lev, int n) { level_steps[lev] = n; }
    //
    // Which step are we at for the specified level?
    //
    int levelCount (int lev) const { return level_count[lev]; }
    //
    // Whether to regrid right after restart
    //
    bool RegridOnRestart () const;
    //
    // Interval between regridding.
    //
    int regridInt (int lev) const { return regrid_int[lev]; }
    //
    // Number of time steps between checkpoint files.
    //
    int checkInt () const { return check_int; }
    //
    // Time between checkpoint files.
    //
    Real checkPer() const { return check_per; }
    //
    // Number of time steps between plot files.
    //
    int plotInt () const { return plot_int; }
    //
    // Time between plot files.
    //
    Real plotPer () const { return plot_per; }
    //
    // Number of time steps between small plot files.
    //
    int smallplotInt () const { return small_plot_int; }
    //
    // Time between plot files.
    //
    Real smallplotPer () const { return small_plot_per; }
    //
    // The names of state variables to output in the
    // plotfile.  They can be set using the amr.plot_vars variable
    // in a ParmParse inputs file. 
    //
    static const std::list<std::string>& statePlotVars () { return state_plot_vars; }
    static const std::list<std::string>& stateSmallPlotVars () { return state_small_plot_vars; }
    //
    // Is the string the name of a variable in state_plot_vars?
    //
    static bool isStatePlotVar (const std::string& name);
    static bool isStateSmallPlotVar (const std::string& name);
    //
    // If the string is not the name of a variable in state_plot_vars,
    // add it to state_plot_vars.  
    //
    static void addStatePlotVar (const std::string& name);
    static void addStateSmallPlotVar (const std::string& name);
    //
    // Remove the string from state_plot_vars.
    //
    static void deleteStatePlotVar (const std::string& name);
    //
    // Clear the list of state_plot_vars.
    //
    static void clearStatePlotVarList ();
    static void clearStateSmallPlotVarList ();
    //
    //  Fill the list of state_plot_vars with all of the state quantities.
    //
    static void fillStatePlotVarList ();
    //
    //  Write out plotfiles (True/False)?
    //
    static bool Plot_Files_Output ();
    //
    // The names of derived variables to output in the
    // plotfile.  They can be set using the amr.derive_plot_vars 
    // variable in a ParmParse inputs file.
    //
    static const std::list<std::string>& derivePlotVars () { return derive_plot_vars; }
    //
    // Is the string the name of a variable in derive_plot_vars?
    //
    static bool isDerivePlotVar (const std::string& name);
    //
    // If the string is not the name of a variable in 
    // derive_plot_vars, add it to derive_plot_vars.  
    //
    static void addDerivePlotVar (const std::string& name);
    //
    // Remove the string from derive_plot_vars.
    //
    static void deleteDerivePlotVar (const std::string& name);
    //
    // Clear the list of derive_plot_vars.
    //
    static void clearDerivePlotVarList ();
    //
    //  Fill the list of derive_plot_vars with all derived quantities.
    //
    static void fillDerivePlotVarList ();

    static void Initialize ();
    static void Finalize ();
    //
    // Number of buffer cells in error estimator.
    //
    int nErrorBuf (int lev) const { return n_error_buf[lev]; }
    //
    // Minimum allowable grid efficiency.
    //
    Real gridEff () const { return grid_eff; }
    //
    // Number of cells to define proper nesting 
    //
    int nProper () const { return n_proper; }
    //
    // List of grids at each level.
    //
    const BoxArray& boxArray (int lev) const;
#ifdef USE_PARTICLES
    //
    // List of particle_grids at each level.
    //
    const BoxArray& ParticleBoxArray (int lev) const;
#endif
    //
    // AmrLevel lev.
    //
    AmrLevel&  getLevel (int lev) { return amr_level[lev]; }
    //
    // PArray of AmrLevels.
    //
    PArray<AmrLevel>& getAmrLevels ();
    //
    // Total number of cells.
    //
    long cellCount ();
    //
    // Number of cells at given level.
    //
    long cellCount (int lev);
    //
    // Total number of grids.
    //
    int numGrids ();
    //
    // Number of grids at given level.
    //
    int numGrids (int lev);
    //
    // More work to be done?
    //
    int okToContinue ();
    //
    // Regrid only!
    //
    void RegridOnly (Real time);
    //
    // Should we regrid this level?
    //
    bool okToRegrid (int level);
    //
    // Array of BoxArrays read in to initially define grid hierarchy
    //
    static const BoxArray& initialBa (int level) 
        { BL_ASSERT(level-1 < initial_ba.size()); return initial_ba[level-1]; }
    //
    // Number of levels at which the grids are initially specified 
    //
    static int initialBaLevels () { return initial_ba.size(); }
    //
    // Do a complete integration cycle.
    //
    virtual void coarseTimeStep (Real stop_time);
    //
    // Do a complete integration cycle and return the coarse dt.
    //
    Real coarseTimeStepDt (Real stop_time);
    //
    // Retrieve derived data. User is responsible for deleting pointer.
    //
    MultiFab* derive (const std::string& name,
                      Real           time,
                      int            lev,
                      int            ngrow);
    //
    // Name of the restart chkpoint file.
    //
    const std::string& theRestartFile () const { return restart_chkfile; }
    //
    // Name of the restart plotfile.
    //
    const std::string& theRestartPlotFile () const { return restart_pltfile; }
    //
    // The ith datalog file.  Do with it what you want.
    //
    std::ostream& DataLog (int i);
    //
    // The filename of the ith datalog file.
    //
    const std::string DataLogName (int i) const { return datalogname[i]; }
    //
    // How many datalogs have been opened
    //
    int NumDataLogs ();
    //
    // Compute the optimal subcycling pattern.
    // This assumes that anything less than cycle_max[i] is a valid
    // number of subcycles at level[i]. For example 
    // if ref_ratio[i] = cycle_max[i] = 4, then 1,2,3,4 are all valid 
    // values for n_cycles[i]
    //
    static Real computeOptimalSubcycling (int   n, 
                                          int*  best,
                                          Real* dt_max, 
                                          Real* est_work, 
                                          int*  cycle_max);
    //
    // Write the plot file to be used for visualization.
    //
    virtual void writePlotFile ();
    int stepOfLastPlotFile () const {return last_plotfile;}
    //
    // Write the small plot file to be used for visualization.
    //
    virtual void writeSmallPlotFile ();
    int stepOfLastSmallPlotFile () const {return last_smallplotfile;}
    //
    // Write current state into a chk* file.
    //
    virtual void checkPoint ();
    int stepOfLastCheckPoint () const {return last_checkpoint;}

    const Array<BoxArray>& getInitialBA();

    //
    // Specialized version:
    // Define BoundaryPointLists that give the intersections
    //    of the external geometry with constant (i,k) and (j,k)
    // These are defined at the coarsest level indexing only.
    //
    void setBoundaryGeometry(BoundaryPointList& IntersectLoX,
                             BoundaryPointList& IntersectHiX,
                             BoundaryPointList& IntersectLoY,
                             BoundaryPointList& IntersectHiY)
    {
        intersect_lox = IntersectLoX;
        intersect_hix = IntersectHiX;
        intersect_loy = IntersectLoY;
        intersect_hiy = IntersectHiY;
    };

    //
    // More general version:
    // Define BoundaryPointLists that give the intersections
    //    of the external geometry with constant (i,k),(j,k)
    //    and (i,j).
    // These are defined at the coarsest level indexing only.
    //
    void setBoundaryGeometry(BoundaryPointList& IntersectLoX,
                             BoundaryPointList& IntersectHiX,
                             BoundaryPointList& IntersectLoY,
                             BoundaryPointList& IntersectHiY,
                             BoundaryPointList& IntersectLoZ,
                             BoundaryPointList& IntersectHiZ)
    {
        intersect_lox = IntersectLoX;
        intersect_hix = IntersectHiX;
        intersect_loy = IntersectLoY;
        intersect_hiy = IntersectHiY;
        intersect_loz = IntersectLoZ;
        intersect_hiz = IntersectHiZ;
    };

    BoundaryPointList& getIntersectLoX() 
    {
        return intersect_lox;
    };
    BoundaryPointList& getIntersectHiX() 
    {
        return intersect_hix;
    };
    BoundaryPointList& getIntersectLoY() 
    {
        return intersect_loy;
    };
    BoundaryPointList& getIntersectHiY() 
    {
        return intersect_hiy;
    };
    BoundaryPointList& getIntersectLoZ() 
    {
        return intersect_loz;
    };
    BoundaryPointList& getIntersectHiZ() 
    {
        return intersect_hiz;
    };

#ifdef USE_PARTICLES
    //
    // Add one particle at a time.
    //
    void addOneParticle (int id_in, int cpu_in, 
                         std::vector<double>& xloc, std::vector<double>& attributes);
    //
    // Redistribute after adding particles
    //
    void RedistributeParticles ();
    //
    // Get particle data out of the Particle Class
    //
    void GetParticleIDs       (Array<int>&  part_ids);
    void GetParticleCPU       (Array<int>&  part_cpu);
    void GetParticleLocations (Array<Real>& part_locs);
    void GetParticleData      (Array<Real>& part_data, int start_comp, int num_comp);

    AmrParGDB* GetParGDB () const { return m_gdb; }
#endif

    void AddProcsToSidecar(int nSidecarProcs, int prevSidecarProcs);
    void AddProcsToComp(int nSidecarProcs, int prevSidecarProcs);
    void RedistributeGrids(int how);
    void BroadcastBoundaryPointList(BoundaryPointList &bpl, int myLocalId, int rootId, MPI_Comm comm);
    void PrintData(std::ostream& os);
    static void BroadcastBCRec(BCRec &bcrec, int myLocalId, int rootId, MPI_Comm localComm);

protected:
    //
    // Used by grid_places().
    //
    static void ProjPeriodic (BoxList&        bd,
                              const Geometry& geom);
    //
    // Initialize grid hierarchy -- called by Amr::init.
    //
    void initialInit (Real strt_time, Real stop_time, 
                      const BoxArray* lev0_grids = 0, const Array<int>* pmap = 0);
    //
    // Read the probin file.
    //
    void readProbinFile (int& init);
    //
    // Check for valid input.
    //
    void checkInput ();
    //
    // Restart from a checkpoint file.
    //
    void restart (const std::string& filename);
    //
    // Define and initialize coarsest level.
    //
    void defBaseLevel (Real start_time, const BoxArray* lev0_grids = 0, const Array<int>* pmap = 0);
    //
    // Define and initialize refined levels.
    //
    void bldFineLevels (Real start_time);
    //
    // Rebuild grid hierarchy finer than lbase.
    //
    virtual void regrid (int  lbase,
                         Real time,
                         bool initial = false);
    //
    // Regrid level 0 on restart. 
    //
    virtual void regrid_level_0_on_restart ();
    //
    // Define new grid locations (called from regrid) and put into new_grids.
    //
    void grid_places (int              lbase,
                      Real             time,
                      int&             new_finest,
                      Array<BoxArray>& new_grids);

    //
    // If ngrids < nprocs, break up the grids to better load balance
    //
    void impose_refine_grid_layout (int              lbase,
                                    int              new_finest,
                                    Array<BoxArray>& new_grids);

    //
    // Do a single timestep on level L.
    //
    virtual void timeStep (int  level,
                           Real time,
                           int  iteration,
                           int  niter,
                           Real stop_time);

    //
    // Whether to write a plotfile now
    //
    bool writePlotNow ();
    bool writeSmallPlotNow ();

    void printGridInfo (std::ostream& os,
                        int           min_lev,
                        int           max_lev);

    void printGridSummary (std::ostream& os,
                           int           min_lev,
                           int           max_lev);

    void setRecordGridInfo (const std::string&);

    void setRecordRunInfo (const std::string&);

    void setRecordRunInfoTerse (const std::string&);

    void setRecordDataInfo (int i, const std::string&);
    
    void initSubcycle();
    void initPltAndChk();

    //
    // The data ...
    //
    std::string      regrid_grids_file;   // Grids file that will bypass regridding.
    std::string      initial_grids_file;  // Grids file that will bypass regridding only at initialization.
    int              max_level;    // Maximum allowed level.
    int              finest_level; // Current finest level.
    Array<IntVect>   ref_ratio;    // Refinement ratios [0:finest_level-1]
    PArray<AmrLevel> amr_level;    // Array of levels
    Array<Geometry>  geom;         // Index extent of domain at all levels.
    Real             cumtime;      // Physical time variable.
    Real             start_time;   // Physical time this simulation started.
    Array<Real>      dt_level;     // Timestep at this level.
    Array<int>       level_steps;  // Number of time steps at this level.
    Array<int>       level_count;
    Array<int>       n_cycle;
    std::string      subcycling_mode; //Type of subcycling to use.
    Array<Real>      dt_min;
    bool             isPeriodic[BL_SPACEDIM];  // Domain periodic?
    Array<int>       regrid_int;      // Interval between regridding.
    Array<int>       n_error_buf;     // Buffer cells around each tagged cell.
    Array<int>       blocking_factor; // Blocking factor in grid generation (by level).
    Array<int>       max_grid_size;   // Maximum allowable grid size (by level).
    Real             grid_eff;        // Grid efficiency.
    int              n_proper;        // # cells required for proper nesting.
    int              last_checkpoint; // Step number of previous checkpoint.
    int              check_int;       // How often checkpoint (# time steps).
    Real             check_per;       // How often checkpoint (units of time).
    std::string      check_file_root; // Root name of checkpoint file.
    int              last_plotfile;   // Step number of previous plotfile.
    int              last_smallplotfile;   // Step number of previous small plotfile.
    int              plot_int;        // How often plotfile (# of time steps)
    Real             plot_per;        // How often plotfile (in units of time)
    int              small_plot_int;  // How often small plotfile (# of time steps)
    Real             small_plot_per;  // How often small plotfile (in units of time)
    int              write_plotfile_with_checkpoint;  // Write out a plotfile whenever we checkpoint
    int              file_name_digits; // How many digits to use in the plotfile and checkpoint names
    int              message_int;     // How often checking messages touched by user, such as "stop_run"
    std::string      plot_file_root;  // Root name of plotfile.
    std::string      small_plot_file_root;  // Root name of small plotfile.

    int              which_level_being_advanced; // Only >=0 if we are in Amr::timeStep(level,...)

#ifdef USE_STATIONDATA
    StationData      station;
#endif

    int              verbose;
    int              record_grid_info;
    int              record_run_info;
    int              record_run_info_terse;
    std::ofstream    gridlog;
    std::ofstream    runlog;
    std::ofstream    runlog_terse;
    PArray<std::fstream> datalog;
    Array<std::string> datalogname;
    int              sub_cycle;
    std::string      restart_chkfile;
    std::string      restart_pltfile;
    std::string      probin_file;
    LevelBld*        levelbld;
    bool             abort_on_stream_retry_failure;
    int              stream_max_tries;
    int              rebalance_grids;

    bool             bUserStopRequest;
    //
    // The static data ...
    //
    static std::list<std::string> state_plot_vars;  // State Vars to dump to plotfile 
    static std::list<std::string> state_small_plot_vars;  // State Vars to dump to small plotfile 
    static std::list<std::string> derive_plot_vars; // Derived Vars to dump to plotfile 
    static bool                   first_plotfile;
    // Array of BoxArrays read in to initially define grid hierarchy
    static Array<BoxArray> initial_ba;
    // Array of BoxArrays read in to externally define grid hierarchy at each regrid
    static Array<BoxArray> regrid_ba;

public:
    // Should we keep the coarser grids fixed (and not regrid those levels) at all?
    static bool useFixedCoarseGrids;
    // Up to what level should we keep the coarser grids fixed (and not regrid those levels)?
    static int useFixedUpToLevel;

    BoundaryPointList intersect_lox;
    BoundaryPointList intersect_loy;
    BoundaryPointList intersect_loz;
    BoundaryPointList intersect_hix;
    BoundaryPointList intersect_hiy;
    BoundaryPointList intersect_hiz;

    static bool                   first_smallplotfile;

#ifdef USE_PARTICLES
    AmrParGDB* m_gdb;
#endif

private:
    //
    // Disallowed
    //
    Amr (const Amr& rhs);
    Amr& operator= (const Amr& rhs);
};

#endif /*_Amr_H_*/
