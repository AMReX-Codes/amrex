#ifndef _Darcy_H_
#define _Darcy_H_

#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <LevelBld.H>
#include <DarcySNES.H>

#include <iostream>

using std::istream;
using std::ostream;

enum StateType { State_Type = 0,
                 NUM_STATE_TYPE };


class DarcyBld
  :
  public LevelBld
{
  virtual void variableSetUp ();
  virtual void variableCleanUp ();
  virtual AmrLevel *operator() ();
  virtual AmrLevel *operator() (Amr&            papa,
                                int             lev,
                                const Geometry& level_geom,
                                const BoxArray& ba,
                                Real            time);
};


//
// AmrLevel-derived class for Darcy flow
//

class Darcy
  :
  public AmrLevel
{
public:
  typedef BaseFab<int> IntFab;
  typedef FabArray<IntFab> MultiIntFab; 


  //
  //Default constructor.  Builds invalid object.
  //
  Darcy ();
  //
  //The basic constructor.
  //
  Darcy (Amr&            papa,
         int             lev,
         const Geometry& level_geom,
         const BoxArray& bl,
         Real            time);

  static void CleanupStatics ();
  static void InitializeStaticVariables ();
  //
  //The destructor.
  //
  virtual ~Darcy ();
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
    types of plot files. For Darcy it has the form: Darcy-Vnnn
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
  //Cleanup data descriptors at end of run.
  //
  static void variableCleanUp ();

  virtual void setPhysBoundaryValues (FArrayBox& dest,
                                      int        state_indx,
                                      Real       time,
                                      int        dest_comp,
                                      int        src_comp,
                                      int        num_comp);
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
  //Initialize data on this level from another Darcy (during regrid).
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

  bool multilevel_advance(Real  t, 
                          Real  dt, 
                          Real& dt_suggest);

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


  void build_layout ();
  void build_snes ();

  static int       NUM_STATE;
  static int       RhoSat, Pressure;

protected:

  static void read_params ();

  void PostStaticInitialize ();

  void BuildMaterialDataServer ();

  Darcy& getLevel (int lev);

  void avgDown ();

  void avgDown (int state_indx);

  void buildMetrics ();

  //
  // The data.
  //
  MultiFab             volume;
  MultiFab             area[BL_SPACEDIM];
  MultiFab             dLogArea[1];
  Array< Array<Real> > radius;
  bool                 post_restart_flag; // true during first timestep after restart
  //
  // Static data members.
  //
  static bool      dump_old;
  static int       radius_grow;
  static int       verbose;
  static Real      init_shrink;
  static Real      change_max;
  static ErrorList err_list;
  static BCRec     phys_bc;
  static Real      dt_suggest_for_next;
  static Real      dt_max;

  static int NUM_PHASES;
  static int NUM_COMPS_PER_PHASE;

  static Layout* layout;
  static MLBoundary* mlb;
  static DarcySNES* snes;
  static int verbose_SNES;
};    

//
// Inlines.
//

inline
MultiFab*
Darcy::Area ()
{
  return area;
}

inline
MultiFab&
Darcy::Area (int dir)
{
  return area[dir];
}

inline
MultiFab&
Darcy::Volume ()
{
  return volume;
}

inline
Darcy&
Darcy::getLevel (int lev)
{
  return *(Darcy *) &parent->getLevel(lev);
}

#endif /*_Darcy_H_*/
