//BL_COPYRIGHT_NOTICE

#ifndef _Amr_H_
#define _Amr_H_ 

//
// $Id: Amr.H,v 1.22 1999-04-22 17:13:19 marc Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <fstream>
using std::ofstream;
#else
#include <fstream.h>
#endif

#include <Box.H>
#include <BoxList.H>
#include <RealBox.H>
#include <Geometry.H>
#include <BoxArray.H>
#include <Array.H>
#include <PArray.H>
#include <FArrayBox.H>
#include <aString.H>
#include <StationData.H>

class AmrLevel;
class LevelBld;

//
//@Man:
//@Memo: Manages hierarchy of levels for time-dependent AMR computations.
/*@Doc:

  The Amr class is designed to manage parts of the computation  which do
  not belong on a single level, like establishing and updating the hierarchy
  of levels, global timestepping, and managing the different AmrLevels
*/

class Amr
{
public:
    //
    //@ManDoc: The constructor.
    //
    Amr ();
    //
    //@ManDoc: The destructor.
    //
    ~Amr ();
    //
    //@ManDoc: Init data after construction. Must be called before timestepping.
    //
    void init (Real strt_time, Real stop_time);
    //
    //@ManDoc: Set the timestep on each level.
    //
    void setDtLevel (const Array<Real>& dt_lev);
    //
    //@ManDoc: Set the dtmin on each level.
    //
    void setDtMin(const Array<REAL>& dt_lev);
    //
    //@ManDoc: Set the cycle count on each level.
    //
    void setNCycle (const Array<int>& mss);
    //
    //@ManDoc: Largest allowable grid.
    //
    int maxGridSize () const;
    //
    //@ManDoc: Subcycle in time?
    //
    int subCycle () const;
    //
    //@ManDoc: Maximum allowable level.
    //
    int maxLevel () const;
    //
    //@ManDoc: Finest currently existing level.
    //
    int finestLevel () const;
    //
    //@ManDoc: Physical time.
    //
    Real cumTime () const;
    //
    //@ManDoc: Time step at specified level.
    //
    Real dtLevel (int level) const;
    //
    //@ManDoc: Array of time steps at all levels.
    //
    const Array<Real>& dtLevel () const;
    //
    //@ManDoc: Geometry object.
    //
    const Geometry& Geom (int level) const;
    //
    //@ManDoc: Number of subcycled time steps.
    //
    int nCycle (int level) const;
    //
    //@ManDoc: Refinement ratio.
    //
    IntVect refRatio (int coarse_level) const;
    //
    //@ManDoc: Maximum refinement ratio in any direction.
    //
    int MaxRefRatio (int level) const;
    //
    //@ManDoc: Refinement ratios between all levels.
    //
    const Array<IntVect>& refRatio () const;
    //
    //@ManDoc: Number of time steps at specified level.
    //
    int levelSteps (int lev) const;
    //
    //@ManDoc: Interval between regridding.
    //
    int regridInt (int lev) const;
    //
    //@ManDoc: Number of time steps between checkpoint files.
    //
    int checkInt () const;
    //
    //@ManDoc: Time between checkpoint files.
    //
    Real checkPer() const;
    //
    //@ManDoc: Number of time steps between plot files.
    //
    int plotInt () const;
    //
    //@ManDoc: Time between plot files.
    //
    Real plotPer () const;

    /*@ManDoc: The names of state or derived variables to output in the
               plotfile.  They can be set using the amr.plot\_vars variable
               in a ParmParse inputs file.  If no variables are specified
               in the inputs file this defaults to a single-element array
               containing the string "ALL".
    */
    const Array<aString>& plotVars () const;

    /*@ManDoc: Is the string the name of a variable in plot\_vars?
               Also returns true if plot\_vars.length() == 1 and
               plot\_vars[0] == "ALL".
    */
    bool isPlotVar (const aString& name) const;
    //
    //@ManDoc: Number of buffer cells in error estimator.
    //
    int nErrorBuf (int lev) const;
    //
    //@ManDoc: Minimum allowable grid efficiency.
    //
    Real gridEff () const;
    //
    //@ManDoc: Number of cells to define proper nesting 
    //
    int nProper () const;
    //
    //@ManDoc: Blocking factor in grid generation.
    //
    int blockingFactor () const;
    //
    //@ManDoc: List of grids at each level.
    //
    const BoxArray& boxArray (int lev) const;
    //
    //@ManDoc: AmrLevel lev.
    //
    AmrLevel&  getLevel (int lev);
    //
    //@ManDoc: PArray of AmrLevels.
    //
    PArray<AmrLevel>& getAmrLevels ();
    //
    //@ManDoc: Total number of cells.
    //
    long cellCount ();
    //
    //@ManDoc: Number of cells at given level.
    //
    long cellCount (int lev);
    //
    //@ManDoc: Total number of grids.
    //
    int numGrids ();
    //
    //@ManDoc: Number of grids at given level.
    //
    int numGrids (int lev);
    //
    //@ManDoc: More work to be done?
    //
    virtual int okToContinue ();
    //
    //@ManDoc: Do a complete integration cycle.
    //
    void coarseTimeStep (Real stop_time);
    //
    //@ManDoc: Retrieve derived data. User is responsible for deleting pointer.
    //
    MultiFab* derive (const aString& name,
                      Real           time,
                      int            lev,
                      int            ngrow);
    //
    //@ManDoc: Write a plotfile to disk.
    //
    void writePlotFile ();
    //
    //@ManDoc: Name of the restart file.
    //
    const aString& theRestartFile () const;
    //
    //@ManDoc: The datalog file.  Do with it what you want.
    //
    ostream& DataLog() { return datalog; }

protected:
    //
    // Write the plot file to be used for visualization.
    //
    void writePlotFile (const aString& root, int num);
    //
    // Initialize grid hierarchy.
    //
    void initialInit (Real strt_time, Real stop_time);
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
    void restart (const aString& filename);
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
                   int  ncycle);

    void printGridInfo (ostream& os,
                        int      min_lev,
                        int      max_lev);

    void setRecordGridInfo (const aString&);

    void setRecordRunInfo (const aString&);

    void setRecordDataInfo (const aString&);
    //
    // The data ...
    //
    aString          grids_file;   // Grids file that will bypass regridding.
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
    int              blocking_factor; // Blocking factor in grid generation.
    Real             grid_eff;        // Grid efficiency.
    int              max_grid_size;   // Maximum allowable grid size.
    int              n_proper;        // # cells required for proper nesting.
    int              last_checkpoint; // Step number of previous checkpoint.
    int              check_int;       // How often checkpoint (# time steps).
    Real             check_per;       // How often checkpoint (units of time).
    aString          check_file_root; // Root name of checkpoint file.
    int              last_plotfile;   // Step number of previous plotfile.
    int              plot_int;        // How often plotfile (# of time steps)
    Real             plot_per;        // How often plotfile (in units of time)
    aString          plot_file_root;  // Root name of plotfile.
    Array<aString>   plot_vars;       // Vars to dump to plotfile or "ALL".
    StationData      station;
    int              verbose;
    int              record_grid_info;
    int              record_run_info;
    ofstream         gridlog;
    ofstream         runlog;
    ofstream         datalog;
    int              sub_cycle;
    aString          restart_file;
    aString          probin_file;
    LevelBld*        levelbld;

private:
    //
    // Disallowed
    //
    Amr (const Amr& rhs);
    Amr& operator= (const Amr& rhs);
};

//
// Inlines.
//

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
const Array<aString>&
Amr::plotVars () const
{
    return plot_vars;
}

inline
int
Amr::maxGridSize () const
{
    return max_grid_size;
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
int
Amr::nCycle (int level) const
{
    return n_cycle[level];
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
Amr::levelSteps (int i) const
{
    return level_steps[i];
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
Amr::blockingFactor () const
{
    return blocking_factor;
}

inline
int
Amr::nProper () const
{
    return n_proper;
}

inline
void
Amr::writePlotFile ()
{
    writePlotFile(plot_file_root,level_steps[0]);
}

inline
const aString&
Amr::theRestartFile () const
{
    return restart_file;
}

#endif /*_Amr_H_*/
