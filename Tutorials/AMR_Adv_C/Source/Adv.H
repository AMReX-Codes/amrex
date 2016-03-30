#ifndef _Adv_H_
#define _Adv_H_

#include <AmrLevel.H>
#include <FluxRegister.H>
#include <BL_CXX11.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

enum StateType { State_Type = 0,
                 NUM_STATE_TYPE };

//
// AmrLevel-derived class for hyperbolic conservation equations
//

class Adv
    :
    public AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    Adv ();
    //
    //The basic constructor.
    //
    Adv (Amr&            papa,
	 int             lev,
	 const Geometry& level_geom,
	 const BoxArray& bl,
	 Real            time);
    //
    //The destructor.
    //
    virtual ~Adv () BL_OVERRIDE;
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (Amr&          papa,
                          std::istream& is,
			  bool          bReadSpecial = false) BL_OVERRIDE;

    /*A string written as the first item in writePlotFile() at
               level zero. It is so we can distinguish between different
               types of plot files. For Adv it has the form: Adv-Vnnn
    */
    virtual std::string thePlotFileType () const BL_OVERRIDE;

    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how) BL_OVERRIDE;
    //
    //Define data descriptors.
    //
    static void variableSetUp ();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData () BL_OVERRIDE;
    //
    //Initialize data on this level from another Adv (during regrid).
    //
    virtual void init (AmrLevel& old) BL_OVERRIDE;
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () BL_OVERRIDE;
    //
    //Proceed with next timestep?
    //
    virtual int okToContinue () BL_OVERRIDE { return 1; }
    //
    //Advance grids at this level in time.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) BL_OVERRIDE;
    //
    //Estimate time step.
    //
    Real estTimeStep (Real dt_old);
    //
    //Compute initial time step.
    //
    Real initialTimeStep ();
    //
    //Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Array<int>&           n_cycle,
                                   const Array<IntVect>& ref_ratio,
                                   Array<Real>&          dt_level,
                                   Real                  stop_time) BL_OVERRIDE;
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
                               int                   post_regrid_flag) BL_OVERRIDE;
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration) BL_OVERRIDE;
    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase, int new_finest) BL_OVERRIDE { ; }
    //
    //Do work after a restart().
    //
    virtual void post_restart () BL_OVERRIDE { ; }
    //
    //Do work after init().
    //
    virtual void post_init (Real stop_time) BL_OVERRIDE;
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
			   int          n_error_buf = 0, int ngrow = 0) BL_OVERRIDE;

    static int  NUM_STATE;
    static int  NUM_GROW;

protected:

    static void read_params ();

    Adv& getLevel (int lev);

    FluxRegister& getFluxReg ();
    FluxRegister& getFluxReg (int lev);

    void reflux ();

    void avgDown ();

    void avgDown (int state_indx);

    //
    // The data.
    //
    FluxRegister*        flux_reg;
    //
    // Static data members.
    //
    static int       verbose;
    static Real      cfl;
    static int       do_reflux;
};    

//
// Inlines.
//

inline
Adv&
Adv::getLevel (int lev)
{
    return *(Adv *) &parent->getLevel(lev);
}

inline
FluxRegister&
Adv::getFluxReg ()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
FluxRegister&
Adv::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}
#endif /*_Adv_H_*/
