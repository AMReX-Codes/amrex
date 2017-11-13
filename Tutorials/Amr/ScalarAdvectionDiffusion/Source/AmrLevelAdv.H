#ifndef _AmrLevelAdv_H_
#define _AmrLevelAdv_H_


#include <AMReX_AmrLevel.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_YAFluxRegister.H>

#include <memory>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

enum StateType { Phi_Type = 0,
                 NUM_STATE_TYPE };

namespace amrex
{
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
    amrex::Real advanceGodunov (amrex::Real time,
                                amrex::Real dt,
                                int  iteration,
                                int  ncycle) ;
    amrex::Real advanceMOLRK2  (amrex::Real time,
                                amrex::Real dt,
                                int  iteration,
                                int  ncycle) ;

    amrex::Real advanceMOLRK3  (Real time,
                                Real dt,
                                int  iteration,
                                int  ncycle);

    amrex::Real  advanceMOLRK4 (Real time,
                                Real dt,
                                int  iteration,
                                int  ncycle);
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
    virtual void computeInitialDt (int                                  finest_level,
                                   int                                  sub_cycle,
                                   amrex::Vector<int>&                  n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                          stop_time) override;
    //
    //Compute new `dt'.
    //
    virtual void computeNewDt (int                                  finest_level,
                               int                                  sub_cycle,
                               amrex::Vector<int>&                  n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                          stop_time,
                               int                                  post_regrid_flag) override;
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

    static int  NUM_STATE;
    static int  NUM_GROW;

    amrex::Real
    compute_dPhiDt_godunov  (const MultiFab& Sborder, MultiFab& dPhiDt, Real time, Real dt,
                             YAFluxRegister* fr_as_crse, YAFluxRegister* fr_as_fine, int iteration);
    amrex::Real
    compute_dPhiDt_MOL2ndOrd (const MultiFab& Sborder, MultiFab& dPhiDt, Real time, Real dt,
                              YAFluxRegister* fr_as_crse, YAFluxRegister* fr_as_fine, int iteration);


    amrex::Real
    compute_dPhiDt_MOL4thOrd (const MultiFab& Sborder, MultiFab& dPhiDt, Real time, Real dt,
                              YAFluxRegister* fr_as_crse, YAFluxRegister* fr_as_fine, int iteration);
                               

  protected:

    static void read_params ();

    // inline functioned implemented below
    AmrLevelAdv& getLevel (int lev);
    amrex::YAFluxRegister& getFluxReg ();
    amrex::YAFluxRegister& getFluxReg (int lev);

    void reflux ();

    void avgDown ();

    void avgDown (int state_indx);

    //
    // The data.
    //
    amrex::YAFluxRegister*        flux_reg;
    //
    // Static data members.
    //
    static int          verbose;
    static amrex::Real  cfl;
    static amrex::Real  diffco;
    static int          do_reflux;

    //only persistent if not on finest level
    MultiFab m_k1;
    MultiFab m_k2;
    MultiFab m_k3;
    MultiFab m_k4;
 
    void
    fillGhostCellsRK4 (MultiFab & a_phiC, //phi on the coarser level
                       const int& a_stage); //rk4 stage to fill ghost for;
    void
    fillGhostCellsRK3 (MultiFab & a_phiC, //phi on the coarser level
                       const int& a_stage); //rk3 stage to fill ghost for;

    bool m_use_fixed_dt, m_use_limiting;
    std::string m_algorithm;
    Real m_fixed_dt;
    void initSwitches();
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
  amrex::YAFluxRegister&
  AmrLevelAdv::getFluxReg ()
  {
    BL_ASSERT(flux_reg);
    return *flux_reg;
  }

  inline
  amrex::YAFluxRegister&
  AmrLevelAdv::getFluxReg (int lev)
  {
    return getLevel(lev).getFluxReg();
  }
}
#endif /*_AmrLevelAdv_H_*/
