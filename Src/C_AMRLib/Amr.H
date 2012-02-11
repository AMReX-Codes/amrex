
#ifndef _Amr_H_
#define _Amr_H_ 

#include <fstream>

#include <Box.H>
#include <Geometry.H>
#include <BoxArray.H>
#include <Array.H>
#include <PArray.H>
#ifdef USE_STATIONDATA
#include <StationData.H>
#endif

class AmrLevel;
class LevelBld;
class BoxDomain;

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
public:
    //
    // The constructor.
    //
    Amr ();
    //
    // The destructor.
    //
    ~Amr ();
    //
    // Init data after construction. Must be called before timestepping.
    //
    virtual void init (Real strt_time, Real stop_time);
    //
    // Set the timestep on each level.
    //
    void setDtLevel (const Array<Real>& dt_lev);
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
    int blockingFactor (int lev) const;
    //
    // Largest allowable grid.
    //
    int maxGridSize (int lev) const;
    //
    // Subcycle in time?
    //
    int subCycle () const;
    //
    // Maximum allowable level.
    //
    int maxLevel () const;
    //
    // Finest currently existing level.
    //
    int finestLevel () const;
    //
    // Physical time.
    //
    Real cumTime () const;
    //
    // Time step at specified level.
    //
    Real dtLevel (int level) const;
    //
    // Array of time steps at all levels.
    //
    const Array<Real>& dtLevel () const;
    //
    // Geometry object.
    //
    const Geometry& Geom (int level) const;
    //
    // Number of subcycled time steps.
    //
    int nCycle (int level) const;
    //
    // Refinement ratio.
    //
    IntVect refRatio (int coarse_level) const;
    //
    // Maximum refinement ratio in any direction.
    //
    int MaxRefRatio (int level) const;
    //
    // Refinement ratios between all levels.
    //
    const Array<IntVect>& refRatio () const;
    //
    // Number of time steps at specified level.
    //
    int levelSteps (int lev) const;
    //
    // Which step are we at for the specified level?
    //
    int levelCount (int lev) const;
    //
    // Whether to regrid right after restart
    //
    bool RegridOnRestart () const;
    //
    // Interval between regridding.
    //
    int regridInt (int lev) const;
    //
    // Number of time steps between checkpoint files.
    //
    int checkInt () const;
    //
    // Time between checkpoint files.
    //
    Real checkPer() const;
    //
    // Number of time steps between plot files.
    //
    int plotInt () const;
    //
    // Time between plot files.
    //
    Real plotPer () const;
    //
    // The names of state variables to output in the
    // plotfile.  They can be set using the amr.plot_vars variable
    // in a ParmParse inputs file. 
    //
    static const std::list<std::string>& statePlotVars ();
    //
    // Is the string the name of a variable in state_plot_vars?
    //
    static bool isStatePlotVar (const std::string& name);
    //
    // If the string is not the name of a variable in state_plot_vars,
    // add it to state_plot_vars.  
    //
    static void addStatePlotVar (const std::string& name);
    //
    // Remove the string from state_plot_vars.
    //
    static void deleteStatePlotVar (const std::string& name);
    //
    // Clear the list of state_plot_vars.
    //
    static void clearStatePlotVarList ();
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
    static const std::list<std::string>& derivePlotVars ();
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
    int nErrorBuf (int lev) const;
    //
    // Minimum allowable grid efficiency.
    //
    Real gridEff () const;
    //
    // Number of cells to define proper nesting 
    //
    int nProper () const;
    //
    // List of grids at each level.
    //
    const BoxArray& boxArray (int lev) const;
    //
    // AmrLevel lev.
    //
    AmrLevel&  getLevel (int lev);
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
    // Do a complete integration cycle.
    //
    virtual void coarseTimeStep (Real stop_time);
    //
    // Retrieve derived data. User is responsible for deleting pointer.
    //
    MultiFab* derive (const std::string& name,
                      Real           time,
                      int            lev,
                      int            ngrow);
    //
    // Name of the restart file.
    //
    const std::string& theRestartFile () const;
    //
    // The ith datalog file.  Do with it what you want.
    //
    std::ostream& DataLog (int i);
    //
    // How many datalogs have been opened
    //
    int NumDataLogs ();

protected:
    //
    // Used by grid_places().
    //
    static void ProjPeriodic (BoxList&        bd,
                              const Geometry& geom);
    //
    // Write the plot file to be used for visualization.
    //
    void writePlotFile (const std::string& root, int num);
    //
    // Initialize grid hierarchy.
    //
    void initialInit (Real strt_time, Real stop_time);
    //
    // Read the probin file.
    //
    void readProbinFile (int& init);
    //
    // Check for valid input.
    //
    void checkInput ();
    //
    // Write current state into a chk* file.
    //
    void checkPoint ();
    //
    // Restart from a checkpoint file.
    //
    void restart (const std::string& filename);
    //
    // Define and initialize coarsest level.
    //
    void defBaseLevel (Real start_time);
    //
    // Define and initialize refined levels.
    //
    void bldFineLevels (Real start_time);
    //
    // Rebuild grid hierarchy finer than lbase.
    //
    void regrid (int  lbase,
                 Real time,
                 bool initial = false);
    //
    // Define new grid locations (called from regrid) and put into new_grids.
    //
    void grid_places (int              lbase,
                      Real             time,
                      int&             new_finest,
                      Array<BoxArray>& new_grids);
    //
    // Do a single timestep on level L.
    //
    void timeStep (int  level,
                   Real dt,
                   int  iter,
                   int  ncycle,
                   Real stop_time);

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
    //
    // The data ...
    //
    std::string      grids_file;   // Grids file that will bypass regridding.
    int              max_level;    // Maximum allowed level.
    int              finest_level; // Current finest level.
    Array<IntVect>   ref_ratio;    // Refinement ratios [0:finest_level-1]
    PArray<AmrLevel> amr_level;    // Array of levels
    Array<Geometry>  geom;         // Index extent of domain at all levels.
    Real             cumtime;      // Physical time variable.
    Array<Real>      dt_level;     // Timestep at this level.
    Array<int>       level_steps;  // Number of time steps at this level.
    Array<int>       level_count;
    Array<int>       n_cycle;
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
    std::string          check_file_root; // Root name of checkpoint file.
    int              last_plotfile;   // Step number of previous plotfile.
    int              plot_int;        // How often plotfile (# of time steps)
    Real             plot_per;        // How often plotfile (in units of time)
    int              file_name_digits; // How many digits to use in the plotfile and checkpoint names
    std::string      plot_file_root;  // Root name of plotfile.
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
    PArray<std::ofstream> datalog;
    int              sub_cycle;
    std::string      restart_file;
    std::string      probin_file;
    LevelBld*        levelbld;
    //
    // The static data ...
    //
    static std::list<std::string> state_plot_vars;  // State Vars to dump to plotfile 
    static std::list<std::string> derive_plot_vars; // Derived Vars to dump to plotfile 
    static bool                   first_plotfile;

private:
    //
    // Disallowed
    //
    Amr (const Amr& rhs);
    Amr& operator= (const Amr& rhs);
};

inline
AmrLevel&
Amr::getLevel (int lev)
{
    return amr_level[lev];
}

inline
int
Amr::blockingFactor (int lev) const
{
    return blocking_factor[lev];
}

inline
int
Amr::maxGridSize (int lev) const
{
    return max_grid_size[lev];
}

inline
int
Amr::maxLevel () const
{
    return max_level;
}

inline
int
Amr::finestLevel () const
{
    return finest_level;
}

inline
IntVect
Amr::refRatio (int level) const
{
    return ref_ratio[level];
}

inline
const Array<IntVect>&
Amr::refRatio () const
{
    return ref_ratio;
}

inline
Real
Amr::dtLevel (int level) const
{
    return dt_level[level];
}

inline
const Array<Real>&
Amr::dtLevel () const
{
    return dt_level;
}

inline
const Geometry&
Amr::Geom (int level) const
{
    return geom[level];
}

inline
int
Amr::checkInt () const
{
    return check_int;
}

inline
Real
Amr::checkPer () const
{
    return check_per;
}

inline
int
Amr::plotInt () const
{
    return plot_int;
}

inline
Real
Amr::plotPer () const
{
    return plot_per;
}

inline
const std::list<std::string>&
Amr::statePlotVars ()
{
    return state_plot_vars;
}

inline
const std::list<std::string>&
Amr::derivePlotVars ()
{
    return derive_plot_vars;
}

inline
int
Amr::nCycle (int level) const
{
    return n_cycle[level];
}

inline
int
Amr::levelSteps (int i) const
{
    return level_steps[i];
}

inline
int
Amr::levelCount (int i) const
{
    return level_count[i];
}

inline
Real
Amr::cumTime () const
{
    return cumtime;
}

inline
int
Amr::regridInt (int lev) const
{
    return regrid_int[lev];
}

inline
int
Amr::nErrorBuf (int lev) const
{
    return n_error_buf[lev];
}

inline
Real
Amr::gridEff () const
{
    return grid_eff;
}

inline
int
Amr::subCycle () const
{
    return sub_cycle;
}

inline
int
Amr::nProper () const
{
    return n_proper;
}

inline
const std::string&
Amr::theRestartFile () const
{
    return restart_file;
}

#endif /*_Amr_H_*/
