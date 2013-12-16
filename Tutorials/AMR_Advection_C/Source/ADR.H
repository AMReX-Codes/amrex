#ifndef _ADR_H_
#define _ADR_H_

#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>

#include <iostream>

using std::istream;
using std::ostream;

enum StateType { State_Type = 0,
                 NUM_STATE_TYPE };

//
// AmrLevel-derived class for hyperbolic conservation equations for stellar media
//

class ADR
    :
    public AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    ADR ();
    //
    //The basic constructor.
    //
    ADR (Amr&            papa,
            int             lev,
            const Geometry& level_geom,
            const BoxArray& bl,
            Real            time);
    //
    //The destructor.
    //
    virtual ~ADR ();
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (Amr&     papa,
                          istream& is,
			  bool bReadSpecial = false);
    //
    //Call AmrLevel::checkPoint and then add radiation info
    //
    virtual void checkPoint(const std::string& dir,
                            std::ostream&      os,
                            VisMF::How         how,
                            bool               dump_old);

    /*A string written as the first item in writePlotFile() at
               level zero. It is so we can distinguish between different
               types of plot files. For ADR it has the form: ADR-Vnnn
    */
    virtual std::string thePlotFileType () const;

    virtual void setPlotVariables ();
    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                ostream&       os,
                                VisMF::How     how);
    //
    //Define data descriptors.
    //
    static void variableSetUp ();
    //
    //Define tagging functions.
    //
    static void ErrorSetUp ();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData ();

    //
    //Get problem metrics.
    //
    MultiFab* Area ();
    MultiFab& Area (int dir);
    //
    //The volume of the multifab.
    //
    MultiFab& Volume ();
    //
    //Set time levels of state data.
    //
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new);
    //
    //Initialize data on this level from another ADR (during regrid).
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
    //Allocate data at old time.
    //
    virtual void allocOldData ();
    //
    //Remove data at old time.
    //
    virtual void removeOldData ();
    //
    //Print information about energy budget.
    //
    void do_energy_diagnostics();
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration);
    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep (Real cumtime);
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
			   int          n_error_buf = 0, int ngrow = 0);

    // Returns a MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the MultiFab is built on the appropriately
    // grown BoxArray.
    MultiFab* derive (const std::string& name,
                      Real               time,
                      int                ngrow);

    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    void derive (const std::string& name,
                 Real               time,
                 MultiFab&          mf,
                 int                dcomp);

    static int numGrow();

#ifdef REACTIONS
    void react_first_half_dt(MultiFab& S_new, 
                             Real time, 
                             Real dt);

    void react_second_half_dt(MultiFab& S_new, 
                              Real time, 
                              Real dt);

    void strang_chem (MultiFab&  state,
                      Real       time,
                      Real       dt);

    void reactState (FArrayBox&  Snew,
                     FArrayBox&  Sold,
                     const Box&  bx,
                     Real        time,
                     Real        dt_react);
#endif

    void getOldSource (Real old_time,                Real dt, MultiFab& Rhs);
    void getNewSource (Real old_time, Real new_time, Real dt, MultiFab& Rhs);

#ifdef DIFFUSION
    void getDiffusionTerm (Real time, MultiFab& DiffTerm, int comp);
    void add_diffusion_to_old_source (MultiFab& S_new, MultiFab& OldDiffTerm, 
                                      Real prev_time, int comp);
    void time_center_diffusion       (MultiFab& S_new, MultiFab& OldDiffTerm, 
                                      Real  cur_time, Real dt, int comp);
#endif

    Real volWgtSum (const std::string& name, Real time);

    Real volWgtSquaredSum (const std::string& name, Real time);

    Real locWgtSum (const std::string& name, Real time, int idir);

    static int       NUM_STATE;
    static int       Density, Xvel, Yvel, Zvel;
    static int       FirstAdv,  LastAdv,  NumAdv;
    static int       FirstSpec, LastSpec, NumSpec;

protected:

    //
    // Need to pass ugdnv out of cmpflx to do level sets
    //
    MultiFab* u_gdnv;

    static void adr_network_init ();
    static void read_params ();

    ADR& getLevel (int lev);

    FluxRegister& getFluxReg ();
    FluxRegister& getFluxReg (int lev);

    void reflux ();

    void enforce_nonnegative_species (MultiFab& S_new);

    void enforce_consistent_e (MultiFab& S);

    void avgDown ();

    void avgDown (int state_indx);

    void buildMetrics ();

    // integrate derived quantities over domain

    Real sumDerive (const std::string& name, Real time);

    void sum_integrated_quantities ();

    //
    // Interpolate cell-centered cync correction from coarse to fine.
    //
    enum SyncInterpType
    {
        PC_T,
        CellCons_T,
        CellConsLin_T,
        CellConsProt_T
    };

    void SyncInterp (MultiFab&      CrseSync,
                     int            c_lev,
                     MultiFab&      FineSync,
                     int            f_lev,
                     IntVect&       ratio,
                     int            src_comp,
                     int            dest_comp,
                     int            num_comp,
                     int            increment,
                     Real           dt_clev,
                     int**          bc_orig_qty,
                     SyncInterpType which_interp = CellCons_T,
                     int            state_comp   = -1);

    //
    // The data.
    //
    MultiFab             volume;
    MultiFab             area[BL_SPACEDIM];
    MultiFab             dLogArea[1];
    Array< Array<Real> > radius;
    FluxRegister*        flux_reg;
    //
    // Static data members.
    //
    static bool      dump_old;
    static int       radius_grow;
    static int       verbose;
    static Real      cfl;
    static Real      init_shrink;
    static Real      change_max;
    static int       do_reflux;
    static ErrorList err_list;
    static BCRec     phys_bc;
    static int       NUM_GROW;

    static int       allow_untagging;
    static int       normalize_species;

#ifdef DIFFUSION
    static int       diffuse_spec;
#endif

    bool             FillPatchedOldState_ok;

    // permits reactions to be turned on and off -- mostly for efficiency's sake
    static int do_react;

    // if true, define an additional source term 
    static int add_ext_src;

    // There can be only one Diffusion object, it covers all levels:
    static class Diffusion *diffusion;
};    

//
// Inlines.
//

inline
int
ADR::numGrow()
{
    return NUM_GROW;
}

inline
MultiFab*
ADR::Area ()
{
    return area;
}

inline
MultiFab&
ADR::Area (int dir)
{
    return area[dir];
}

inline
MultiFab&
ADR::Volume ()
{
    return volume;
}

inline
ADR&
ADR::getLevel (int lev)
{
    return *(ADR *) &parent->getLevel(lev);
}

inline
FluxRegister&
ADR::getFluxReg ()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
FluxRegister&
ADR::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}
#endif /*_ADR_H_*/
