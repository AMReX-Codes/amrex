#ifndef _Amr_H_
#define _Amr_H_ 

//
// $Id: Amr.H,v 1.2 1997-11-22 00:44:22 lijewski Exp $
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

///  Manages hierarchy of levels for time-dependent AMR computations

class Amr {

/*  The Amr class is designed to manage parts of the computation
    which do not belong on a single level, like establishing and updating
    the hierarchy of levels, global timestepping, and managing the different
    AmrLevels
*/

public:

/// administrative functions

//  default constructor
    Amr();

//  delete all memory
    ~Amr();

//  initialize amr data after construction - must be called before timestepping
    void init( );

//  set initial time variables
    void setInitialTime(REAL strt_time);

/// undocumented hacks
///
    void setDtLevel(const Array<REAL>& dt_lev);
///
    void setNCycle(const Array<int>& mss);

/// access functions

//  largest allowable grid
    int  maxGridSize() const;

//  subcycle in time?
    int  subCycle() const;

//  maximum allowable level
    int  maxLevel() const;

//  finest currently existing level
    int  finestLevel() const;

//  physical time
    REAL cumTime() const;

//  time step at specified level
    REAL dtLevel(int level) const;

//  array of time steps at all levels
    const Array<REAL>& dtLevel() const;

//  geometry object
    const Geometry& Geom(int level) const;

//  number of subcycled time steps
    int  nCycle(int level) const;

//  refinement ratio
    IntVect refRatio( int coarse_level ) const;

//  maximum refinement ratio in any direction
    int MaxRefRatio( int level ) const;

//  refinement ratios between all levels
    const Array<IntVect>& refRatio() const;

//  number of time steps at specified level
    int  levelSteps(int lev) const;

//  interval between regridding
    int  regridInt(int lev) const;

//  number of time steps between checkpoint files
    int  checkInt() const {
       return check_int;
    }

//  time between checkpoint files
    REAL  checkPer() const {
       return check_per;
    }

//  number of time steps between plot files
    int  plotInt() const {
       return plot_int;
    }

//  time between plot files
    REAL  plotPer() const {
       return plot_per;
    }

//  number of buffer cells in error estimator
    int  nErrorBuf(int lev) const;

//  minimum allowable grid efficiency
    REAL gridEff() const;

//  blocking factor in grid generation
    int  blockingFactor() const;

//  list of grids at each level
    const BoxArray& boxArray(int lev) const;

//  AmrLevel lev
    AmrLevel&  getLevel(int lev);

//  PArray of AmrLevels
    PArray<AmrLevel>&  getAmrLevels();

/// info functions

//  total number of cells
    long cellCount();

//  number of cells at given level
    long cellCount(int lev);

//  total number of grids
    int  numGrids();

//  number of grids at given level
    int  numGrids(int lev);

///  query
//
    virtual int okToContinue();
///
//  do a complete integration cycle
    void  coarseTimeStep(REAL stop_time);

//  retrieve derived data
    PArray<FARRAYBOX>* derive(const aString& name, REAL time, int lev);

/*   Alternate form.  Note: the user is responsible for deleting this pointer
     when done with it
*/
    FARRAYBOX* derive(const BOX& b, const aString& name, REAL time, int lev);

    void writePlotFile();
    //
    // Print out some Amr-wide I/O statistics.
    //
    static void PrintIOStats (ostream& os);

protected:

/// these functions are internal worker routines

//  write the plot file to be used for visualization
    void writePlotFile(const aString& root, int num);

//  initialize grid hierarchy
    void initialInit();

//  check for valid input
    void checkInput();

//  write current state into a chk* file
    void checkPoint();

//  restart from a checkpoint file
    void restart(const aString& filename);

//  define and initialize coarsest level
    void defBaseLevel(REAL start_time);

//  define and initialize refined levels
    void bldFineLevels(REAL start_time);

//  rebuild grid hierarchy finer than lbase
    void regrid(int lbase, REAL time);

//  define new grid locations (called from regrid) and put into new_grids
    void grid_places(int lbase, REAL time, int& new_finest,
                     Array<BoxArray>& new_grids);

//  do a single timestep on level L
    void timeStep(int level, REAL dt, int iter, int ncycle);
// 
    void printGridInfo(ostream& os, int min_lev, int max_lev);

/// function to set data members
    void setRecordGridInfo(const aString&);
/// 
    void setRecordRunInfo(const aString&);

protected:

//  name of grids file that will bypass regridding
    aString grids_file;

//  maximum allowed level
    int max_level;

//  current finest level
    int finest_level;

//  array of refinement ratios [0:finest_level-1]
    Array<IntVect> ref_ratio;

//  array of levels
    PArray<AmrLevel> amr_level;

//  index extent of domain at all levels
    Array<Geometry>  geom;

// physical time variable
    REAL cumtime;

//  timestep at this level
    Array<REAL> dt_level;

//  number of time steps at this level
    Array<int>  level_steps;
///
    Array<int>  level_count;
///
    Array<int>  n_cycle;
///
    Array<REAL> dt_min;

//  is domain periodic in any index direction?
    bool isPeriodic[BL_SPACEDIM];

//  regridding variables

//  interval between regridding
    Array<int> regrid_int;

//  number of buffer cells around each tagged cell also tagged for refinement
    Array<int> n_error_buf;

//  blocking factor in grid generation
    int blocking_factor;

//  grid efficiency
    REAL grid_eff;

//  maximum allowable grid size
    int max_grid_size;

//  number of cells required for proper nesting of levels
    int n_proper;

///  checkpoint variables

//   step number of previous checkpoint
     int      last_checkpoint;

//   how often to checkpoint (measured in number of time steps)
     int      check_int;

//   how often to checkpoint (measured in units of time)
     REAL     check_per;

//   root name of checkpoint file
     aString  check_file_root;

///  plotfile variables

//   step number of previous plotfile
     int      last_plotfile  ;

//   how often to write plotfile (measured in number of time steps)
     int      plot_int;

//   how often to write plotfile (measured in units of time)
     REAL      plot_per;

//   root name of plotfile
     aString  plot_file_root;

/// logging information and flags - these are static functions so can be used without object

    int trace;
///
    int debug;
///
    int verbose;
///
    int silent;
///
    int record_grid_info;
///
    int record_run_info;
///
    ofstream gridlog;
///
    ofstream runlog;
///
    int sub_cycle;
///
    aString restart_file;

    LevelBld* levelbld;

private:
    //
    // We accumulate # bytes written to disk in pltfile & chkfile headers.
    // VisMF::TheBytesWrittenToDisk() accumulates the actual MF byte count.
    //
    static double BytesWrittenToDisk;

    Amr(const Amr&);
    Amr& operator=(const Amr&);
};

//
// Inlines.
//

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
REAL
Amr::dtLevel (int level) const
{
    return dt_level[level];
}

inline
const Array<REAL>&
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
REAL
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
REAL
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

#endif
