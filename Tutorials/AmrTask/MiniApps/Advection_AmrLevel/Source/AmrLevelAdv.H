#ifndef _AmrLevelAdv_H_
#define _AmrLevelAdv_H_

#ifdef AMREX_PARTICLES
#include <AMReX_AmrParticles.H>
#endif

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>

#include <memory>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

enum StateType { Phi_Type = 0,
                 NUM_STATE_TYPE };

//
// AmrLevel-derived class for hyperbolic conservation equations
//

class AmrLevelAdv
    :
    public amrex::AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    AmrLevelAdv ();
    //
    //The basic constructor.
    //
    AmrLevelAdv (amrex::Amr&     papa,
	         int             lev,
                 const amrex::Geometry& level_geom,
                 const amrex::BoxArray& bl,
                 const amrex::DistributionMapping& dm,
                 amrex::Real            time);
    //
    //The destructor.
    //
    virtual ~AmrLevelAdv () override;
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (amrex::Amr&   papa,
                          std::istream& is,
			  bool          bReadSpecial = false) override;

    virtual void checkPoint (const std::string& dir,
			     std::ostream&      os,
			     amrex::VisMF::How  how = amrex::VisMF::NFiles,
			     bool               dump_old = true) override;

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
    //Initialize data on this level from another AmrLevelAdv (during regrid).
    //
    virtual void init (amrex::AmrLevel& old) override;
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () override;
    //
    //Advance grids at this level in time.
    //
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle) override;
    //
    //Estimate time step.
    //
    amrex::Real estTimeStep (amrex::Real dt_old);
    //
    //Compute initial time step.
    //
    amrex::Real initialTimeStep ();
    //
    //Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   amrex::Vector<int>&           n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                  stop_time) override;
    //
    //Compute new `dt'.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               amrex::Vector<int>&           n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                  stop_time,
                               int                   post_regrid_flag) override;
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration) override;

    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase, int new_finest) override;
    //
    //Do work after a restart().
    //
    virtual void post_restart () override;
    //
    //Do work after init().
    //
    virtual void post_init (amrex::Real stop_time) override;
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           amrex::Real         time,
			   int          n_error_buf = 0, int ngrow = 0) override;

#ifdef AMREX_PARTICLES
    static amrex::AmrTracerParticleContainer* theTracerPC () { return TracerPC.get(); }
#endif

    static int  NUM_STATE;
    static int  NUM_GROW;

protected:

    static void read_params ();

    // inline functioned implemented below
    AmrLevelAdv& getLevel (int lev);
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
    static int          verbose;
    static amrex::Real  cfl;
    static int          do_reflux;

#ifdef AMREX_PARTICLES
    void init_particles ();
    static int       do_tracers;
    static std::unique_ptr<amrex::AmrTracerParticleContainer> TracerPC;
#endif
};    

//
// Inlines.
//

inline
AmrLevelAdv&
AmrLevelAdv::getLevel (int lev)
{
    return *(AmrLevelAdv *) &parent->getLevel(lev);
}

inline
amrex::FluxRegister&
AmrLevelAdv::getFluxReg ()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
amrex::FluxRegister&
AmrLevelAdv::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}
#endif /*_AmrLevelAdv_H_*/
