//BL_COPYRIGHT_NOTICE

#ifndef _Amr_H_
#define _Amr_H_ 

//
// $Id: Amr.H,v 1.9 1998-01-06 23:43:56 lijewski Exp $
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

class AmrLevel;
class LevelBld;

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

    Amr ();

    ~Amr ();
    //
    // Init ata after construction - must be called before timestepping.
    //
    void init ();

    void setDtLevel (const Array<Real>& dt_lev);

    void setNCycle (const Array<int>& mss);
    //
    // Largest allowable grid.
    //
    int maxGridSize () const;
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
    //  Physical time.
    //
    Real cumTime () const;
    //
    // Time step at specified level.
    //
    Real dtLevel (int level) const;
    //
    // Array of time steps at all levels.
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
    // Interval between regridding.
    //
    int  regridInt(int lev) const;
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
    int plotInt() const;
    //
    // Time between plot files.
    //
    Real plotPer() const;
    //
    // Number of buffer cells in error estimator.
    //
    int nErrorBuf (int lev) const;
    //
    // Minimum allowable grid efficiency.
    //
    Real gridEff () const;
    //
    // Blocking factor in grid generation.
    //
    int blockingFactor () const;
    //
    // List of grids at each level.
    //
    const BoxArray& boxArray (int lev) const;
    //
    // AmrLevel lev.
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

    virtual int okToContinue();
    //
    // Do a complete integration cycle.
    //
    void coarseTimeStep (Real stop_time);
    //
    // Retrieve derived data.
    //
    // The user is responsible for deleting this pointer when done with it.
    //
    MultiFab* derive (const aString& name,
                      Real           time,
                      int            lev);
    //
    // Retrieve derived data.
    //
    // The user is responsible for deleting this pointer when done with it.
    //
    FArrayBox* derive (const Box&     b,
                       const aString& name,
                       Real           time,
                       int            lev);

    void writePlotFile ();
    //
    // Name of the restart file.
    //
    const aString& theRestartFile () const;

protected:
    //
    // Write the plot file to be used for visualization.
    //
    void writePlotFile (const aString& root, int num);
    //
    // Initialize grid hierarchy.
    //
    void initialInit ();
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
    void restart(const aString& filename);
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
                 Real time);
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

    void setRecordDataInfo(const aString&);
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
    int              trace;
    int              debug;
    int              verbose;
    int              silent;
    int              record_grid_info;
    int              record_run_info;
    ofstream         gridlog;
    ofstream         runlog;
    ofstream         datalog;
    int              sub_cycle;
    aString          restart_file;
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

#endif
