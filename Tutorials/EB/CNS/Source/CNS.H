#ifndef CNS_H_
#define CNS_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_ErrorList.H>

class CNS
    :
    public amrex::AmrLevel
{
public:

    CNS ();
    CNS (amrex::Amr&            papa,
         int                    lev,
         const amrex::Geometry& level_geom,
         const amrex::BoxArray& bl,
         const amrex::DistributionMapping& dm,
         amrex::Real            time);
    virtual ~CNS ();

    CNS (const CNS& rhs) = delete;
    CNS& operator= (const CNS& rhs) = delete;
    
    // Restart from a checkpoint file.
    virtual void restart (amrex::Amr&     papa,
                          std::istream&   is,
                          bool            bReadSpecial = false) override;

    // Write checkpoint
    virtual void checkPoint(const std::string& dir,
                            std::ostream&      os,
                            amrex::VisMF::How  how = amrex::VisMF::NFiles,
                            bool               dump_old = true);

    // Write a plotfile to specified directory.
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                amrex::VisMF::How  how);

    // Initialize data on this level from another CNS (during regrid).
    virtual void init (amrex::AmrLevel& old) override;

    // Initialize data on this level after regridding if old level did not previously exist
    virtual void init () override;
    
    // Initialize grid data at problem start-up.
    virtual void initData () ;

    // Advance grids at this level in time.
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle) override;

    virtual void computeInitialDt (int                                 finest_level,
                                   int                                 sub_cycle,
                                   amrex::Array<int>&                  n_cycle,
                                   const amrex::Array<amrex::IntVect>& ref_ratio,
                                   amrex::Array<amrex::Real>&          dt_level,
                                   amrex::Real                         stop_time) override;

    virtual void computeNewDt (int                                 finest_level,
                               int                                 sub_cycle,
                               amrex::Array<int>&                  n_cycle,
                               const amrex::Array<amrex::IntVect>& ref_ratio,
                               amrex::Array<Real>&                 dt_min,
                               amrex::Array<Real>&                 dt_level,
                               amrex::Real                         stop_time,
                               int                                 post_regrid_flag) override;

    virtual void post_regrid (int lbase, int new_finest) override;

    // Do work after timestep().
    virtual void post_timestep (int iteration) override;

    // Do work after init().
    virtual void post_init (amrex::Real stop_time) override;

    // Error estimation for regridding.
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int                 clearval,
                           int                 tagval,
                           amrex::Real         time,
                           int                 n_error_buf = 0,
                           int                 ngrow = 0) override;

    // Define data descriptors.
    static void variableSetUp ();

    // Define tagging functions.
    static void ErrorSetUp ();

    // Cleanup data descriptors at end of run.
    static void variableCleanUp ();

protected:

    static void read_params ();
    static int numGrow() { return NUM_GROW; };

    const amrex::MultiFab& volFrac () const { return volfrac; }

    CNS& getLevel (int lev) { return dynamic_cast<CNS&>(parent->getLevel(lev)); }

    void avgDown ();

    void buildMetrics ();

#if 0
    void init_eb (const amrex::Geometry& level_geom,
                  const amrex::BoxArray& ba,
                  const amrex::DistributionMapping& dm);
#endif

    amrex::Real estTimeStep ();

    // Compute initial time step.
    amrex::Real initialTimeStep ();

    void computeTemp (amrex::MultiFab& State, int ng);

    void compute_dSdt (amrex::MultiFab& S, amrex::MultiFab& dSdt, amrex::Real dt);

    amrex::MultiFab volfrac;

    static constexpr int NUM_GROW = 4;
    enum StateVariable {
        Density = 0, Xmom, Ymom, Zmom, Eden, Eint, Temp,
        NUM_STATE
    };

    enum StateDataType {
        State_Type = 0,
        NUM_STATEDATA_TYPE
    };

    static amrex::ErrorList err_list;
    static amrex::BCRec     phys_bc;

    // Parameters
    static int verbose;
    static amrex::IntVect hydro_tile_size;
    static amrex_real cfl;
};

#endif
