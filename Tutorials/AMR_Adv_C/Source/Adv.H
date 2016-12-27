#ifndef _Adv_H_
#define _Adv_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>

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
    public amrex::AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    Adv ();
    //
    //The basic constructor.
    //
    Adv (amrex::Amr&     papa,
	 int             lev,
	 const amrex::Geometry& level_geom,
	 const amrex::BoxArray& bl,
	 const amrex::DistributionMapping& dm,
	 Real            time);
    //
    //The destructor.
    //
    virtual ~Adv () override;
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (amrex::Amr&   papa,
                          std::istream& is,
			  bool          bReadSpecial = false) override;

    /*A string written as the first item in writePlotFile() at
               level zero. It is so we can distinguish between different
               types of plot files. For Adv it has the form: Adv-Vnnn
    */
    virtual std::string thePlotFileType () const override;

    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                amrex::VisMF::How  how) override;
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
    virtual void initData () override;
    //
    //Initialize data on this level from another Adv (during regrid).
    //
    virtual void init (amrex::AmrLevel& old) override;
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () override;
    //
    //Proceed with next timestep?
    //
    virtual int okToContinue () override { return 1; }
    //
    //Advance grids at this level in time.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) override;
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
                                   amrex::Array<int>&           n_cycle,
                                   const amrex::Array<amrex::IntVect>& ref_ratio,
                                   amrex::Array<Real>&          dt_level,
                                   Real                  stop_time) override;
    //
    //Compute new `dt'.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               amrex::Array<int>&           n_cycle,
                               const amrex::Array<amrex::IntVect>& ref_ratio,
                               amrex::Array<Real>&          dt_min,
                               amrex::Array<Real>&          dt_level,
                               Real                  stop_time,
                               int                   post_regrid_flag) override;
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration) override;
    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase, int new_finest) override { ; }
    //
    //Do work after a restart().
    //
    virtual void post_restart () override { ; }
    //
    //Do work after init().
    //
    virtual void post_init (Real stop_time) override;
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
			   int          n_error_buf = 0, int ngrow = 0) override;

    static int  NUM_STATE;
    static int  NUM_GROW;

protected:

    static void read_params ();

    Adv& getLevel (int lev);

    amrex::FluxRegister& getFluxReg ();
    amrex::FluxRegister& getFluxReg (int lev);

    void reflux ();

    void avgDown ();

    void avgDown (int state_indx);

    //
    // The data.
    //
    amrex::FluxRegister*        flux_reg;
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
amrex::FluxRegister&
Adv::getFluxReg ()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
amrex::FluxRegister&
Adv::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}
#endif /*_Adv_H_*/
