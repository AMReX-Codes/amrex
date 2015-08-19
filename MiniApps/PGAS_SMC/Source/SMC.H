#ifndef _SMC_H_
#define _SMC_H_

#include <BC_TYPES.H>
#include <AmrLevel.H>

#include <iostream>
#include <iomanip>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using std::istream;
using std::ostream;

enum StateType {State_Type = 0,
                NUM_STATE_TYPE};

class SMC
    :
    public AmrLevel
{

public:
    //
    //Default constructor.  Builds invalid object.
    //
    SMC () {;}
    //
    //The basic constructor.
    //
    SMC (Amr&            papa,
	 int             lev,
	 const Geometry& level_geom,
	 const BoxArray& bl,
	 Real            time);
    //
    //The destructor.
    //
    virtual ~SMC ();

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // begin pure virtual functions
    //
    /*A string written as the first item in writePlotFile() at
      level zero. It is so we can distinguish between different
      types of plot files. For SMC it has the form: SMC-Vnnn
    */
    virtual std::string thePlotFileType () const;
    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
				ostream&       os,
				VisMF::How     how);
    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData ();
    //
    //Initialize data on this level from another SMC (during regrid).
    //
    virtual void init (AmrLevel& old);
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init ();
    //
    //Proceed with next timestep?
    //
    virtual int okToContinue ();
    //
    //Advance grids at this level in time.
    //
    virtual Real advance (Real time,
			  Real dt,
			  int  iteration,
			  int  ncycle);
    //
    //Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
				   int                   sub_cycle,
				   Array<int>&           n_cycle,
				   const Array<IntVect>& ref_ratio,
				   Array<Real>&          dt_level,
				   Real                  stop_time);
    //
    //Compute new `dt'.
    //
    virtual void computeNewDt (int                   finest_level,
			       int                   sub_cycle,
			       Array<int>&           n_cycle,
			       const Array<IntVect>& ref_ratio,
			       Array<Real>&          dt_min,
			       Array<Real>&          dt_level,
			       Real                  stop_time,
			       int                   post_regrid_flag);
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration);
    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase,
			      int new_finest);   
    //
    //Do work after a restart().
    //
    virtual void post_restart ();
    //
    //Do work after init().
    //
    virtual void post_init (Real stop_time);
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (TagBoxArray& tb,
			   int          clearval,
			   int          tagval,
			   Real         time,
			   int          n_error_buf = 0, int ngrow = 0)
	{;}
    //
    // end pure virtual functions
    // $$$$$$$$$$$$$$$$$$$$$$$$$$

    //
    //Restart from a checkpoint file.
    //
    virtual void restart (Amr&     papa,
			  istream& is,
			  bool bReadSpecial = false);
    //
    //Call AmrLevel::checkPoint
    //
    virtual void checkPoint(const std::string& dir,
			    std::ostream&      os,
			    VisMF::How         how,
			    bool               dump_old);        
    //
    // Modify list of variables to be plotted
    // 
    virtual void setPlotVariables ();
    //
    //Define data descriptors.
    //
    static void variableSetUp ();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();    
    //
    //Estimate time step.
    //
    Real estTimeStep (Real dt_old);
    //
    //Compute initial time step.
    //
    Real initialTimeStep ();
    
    static int ncons, nprim;

    MultiFab prim;

protected:

    static void read_params ();

//    static void init_network();
    
    SMC& getLevel (int lev);

    // static data members
    static int  verbose;    
    static int  cfl_int;
    static Real cfl;
    static Real init_shrink;
    static Real change_max;
};    

//
// Inlines.
//

inline
SMC&
SMC::getLevel (int lev)
{
    return *(SMC *) &parent->getLevel(lev);
}

#endif
