#ifndef _AmrLevelAdv_H_
#define _AmrLevelAdv_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_PhysBCFunct.H>

#ifdef AMREX_PARTICLES
#include <AMReX_AmrParticles.H>
#endif

#include <memory>
#include <iostream>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

// Problem-dependent struct. Defined under Exec directory.
#include "Prob_Parm.H"

enum StateType { Phi_Type = 0,
                 NUM_STATE_TYPE };

/**
 * AmrLevel-derived class for hyperbolic conservation equations
 */
class AmrLevelAdv
    :
    public amrex::AmrLevel
{
public:

    /**
     * Default constructor.  Builds invalid object.
     */
    AmrLevelAdv ();

    /**
     * The basic constructor.
     */
    AmrLevelAdv (amrex::Amr&     papa,
                 int             lev,
                 const amrex::Geometry& level_geom,
                 const amrex::BoxArray& bl,
                 const amrex::DistributionMapping& dm,
                 amrex::Real            time);

    /**
     * The destructor.
     */
    virtual ~AmrLevelAdv () override;

    /**
     * Restart from a checkpoint file.
     */
    virtual void restart (amrex::Amr&   papa,
                          std::istream& is,
                          bool          bReadSpecial = false) override;

    /**
     * Write a checkpoint file.
     */
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             amrex::VisMF::How  how = amrex::VisMF::NFiles,
                             bool               dump_old = true) override;

    /**
     * Write a plotfile to specified directory.
     */
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                amrex::VisMF::How  how) override;

    /**
     * Define data descriptors.
     */
    static void variableSetUp ();

    /**
     * Cleanup data descriptors at end of run.
     */
    static void variableCleanUp ();

    /**
     * Initialize grid data at problem start-up.
     */
    virtual void initData () override;

    /**
     * Initialize data on this level from another AmrLevelAdv (during regrid).
     */
    virtual void init (amrex::AmrLevel& old) override;

    /**
     * Initialize data on this level after regridding if old level did not previously exist
     */
    virtual void init () override;

    /**
     * Advance grids at this level in time.
     */
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle) override;

    /**
     * Estimate time step.
     */
    amrex::Real estTimeStep (amrex::Real dt_old);

    /**
     * Compute initial time step.
     */
    amrex::Real initialTimeStep ();

    /**
     * Compute initial `dt'.
     */
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   amrex::Vector<int>&           n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                  stop_time) override;

    /**
     * Compute new `dt'.
     */
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               amrex::Vector<int>&           n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                  stop_time,
                               int                   post_regrid_flag) override;

    /**
     * Do work after timestep().
     */
    virtual void post_timestep (int iteration) override;

    /**
     * Do work after regrid().
     */
    virtual void post_regrid (int lbase, int new_finest) override;

    /**
     * Do work after a restart().
     */
    virtual void post_restart () override;

    /**
     * Do work after init().
     */
    virtual void post_init (amrex::Real stop_time) override;

    /**
     * Error estimation for regridding.
     */
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int                 clearval,
                           int                 tagval,
                           amrex::Real         time,
                           int          n_error_buf = 0, int ngrow = 0) override;

#ifdef AMREX_PARTICLES
    static amrex::AmrTracerParticleContainer* theTracerPC () { return TracerPC.get(); }
#endif

    static int  NUM_STATE;
    static int  NUM_GROW;

    // ProbParm class instances
    static ProbParm* h_prob_parm;
    static ProbParm* d_prob_parm;

    // Variables for error tagging
    static int max_phierr_lev;
    static int max_phigrad_lev;
    static amrex::Vector<amrex::Real> phierr;
    static amrex::Vector<amrex::Real> phigrad;

#if !defined(AMREX_USE_CUDA)
protected:
#endif

    // Function to read parameters from input file.
    static void read_params ();

    /**
     * Function to read tagging parameters from input file.
     * See Tagging_params.cpp for implementation.
     */
    static void get_tagging_params ();

    /**
     * Function to perform advection for a single level.
     * See Adv.cpp for implementation.
     */
    AMREX_GPU_HOST
    static void advect (const amrex::Real /*time*/,
                        const amrex::Box& bx,
                        amrex::GpuArray<amrex::Box,BL_SPACEDIM> nbx,
                        const amrex::FArrayBox& statein,
                        amrex::FArrayBox& stateout,
                        AMREX_D_DECL(const amrex::FArrayBox& vx,
                                     const amrex::FArrayBox& vy,
                                     const amrex::FArrayBox& vz),
                        AMREX_D_DECL(amrex::FArrayBox& fx,
                                     amrex::FArrayBox& fy,
                                     amrex::FArrayBox& fz),
                        amrex::GpuArray<amrex::Real,BL_SPACEDIM> dx,
                        const amrex::Real dt);

    // inline functions implemented below
    AmrLevelAdv& getLevel (int lev);
    amrex::FluxRegister& getFluxReg ();
    amrex::FluxRegister& getFluxReg (int lev);

    void reflux ();

    void avgDown ();

    void avgDown (int state_indx);

    /*
     * The data.
     */
    amrex::FluxRegister*        flux_reg;

    /*
     * Static data members.
     */
    static int          verbose;
    static amrex::Real  cfl;
    static int          do_reflux;

#ifdef AMREX_PARTICLES
    void init_particles ();
    static int       do_tracers;
    static std::unique_ptr<amrex::AmrTracerParticleContainer> TracerPC;
#endif

};

/**
 * Boundary condition function that does nothing.
 * See bc_nullfill.cpp for implementation.
 */
void nullfill (amrex::Box const& bx, amrex::FArrayBox& data,
               const int dcomp, const int numcomp,
               amrex::Geometry const& geom, const amrex::Real time,
               const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
               const int scomp);

/**
 * Problem-dependent function that computes Godunov velocities for each face.
 * See face_velocity_?D_K.H under Exec directory for implementation.
 */
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void get_face_velocity(const amrex::Real time,
                       amrex::FArrayBox& vx,
                       amrex::FArrayBox& vy,
                       amrex::GpuArray<amrex::Real,BL_SPACEDIM> dx,
                       amrex::GpuArray<amrex::Real,BL_SPACEDIM> prob_lo);

/*
 * Inlines.
 */

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
