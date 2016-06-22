
#ifndef _AmrLevel_H_
#define _AmrLevel_H_ 

#include <REAL.H>
#include <Geometry.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <Derive.H>
#include <BCRec.H>
#include <Interpolater.H>
#include <Amr.H>
#include <SlabStat.H>
#include <DistributionMapping.H>
#include <StateDescriptor.H>
#include <StateData.H>
#include <VisMF.H>

#include <map>

class TagBox;
class TagBoxArray;

//
// Virtual base class for managing individual levels.
//
// AmrLevel functions both as a container for state data on a level
// and also manages the advancement of data in time.
//

class AmrLevel
{
  friend class FillPatchIterator;
  friend class FillPatchIteratorHelper;

public:
    //
    // What time are we at?
    //
    enum TimeLevel { AmrOldTime,
                     AmrHalfTime,
                     AmrNewTime,
                     Amr1QtrTime,
                     Amr3QtrTime,
                     AmrOtherTime };
    //
    // The destructor.
    //
    virtual ~AmrLevel ();
    //
    // A string written as the first item in writePlotFile() at
    // level zero. This MUST be defined by each derived class.
    // It is so we can distinguish between different types of
    // plot files.  This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual std::string thePlotFileType () const = 0;
    //
    // Write plot file stuff to specified directory.  This is a
    // pure virtual function and hence MUST be implemented by
    // derived classes.
    //
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how = VisMF::OneFilePerCPU) = 0;
    //
    // Write small plot file stuff to specified directory.  
    // Unlike writePlotFile, this is NOT a pure virtual function 
    // so implementation by derived classes is optional.
    //
    virtual void writeSmallPlotFile (const std::string& dir,
                                     std::ostream&      os,
				     VisMF::How         how = VisMF::OneFilePerCPU) {};
    //
    // Write current state to checkpoint file.
    //
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             VisMF::How         how = VisMF::OneFilePerCPU,
                             bool               dump_old = true);
    //
    // Restart from a checkpoint file.
    //
    virtual void restart (Amr&          papa,
                          std::istream& is,
			  bool          bReadSpecial = false);

    //
    // Old checkpoint may have different number of states than the new source code.
    //
    virtual void set_state_in_checkpoint (Array<int>& state_in_checkpoint);

#ifdef USE_PARTICLES
    //
    // Make the particle distribution map
    //
    void make_particle_dmap ();
#endif

    //
    // Is name a state variable?
    //
    static bool isStateVariable (const std::string& name,
                                int&               state_indx,
                                int&               ncomp);

    static void FlushFPICache ();
    //
    // Compute the initial time step.  This is a pure virtual function
    // and hence MUST be implemented by derived classes.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Array<int>&           n_cycle,
                                   const Array<IntVect>& ref_ratio,
                                   Array<Real>&          dt_level,
                                   Real                  stop_time) = 0;
    //
    // Compute the next time step.  This is a pure virtual function
    // and hence MUST be implemented by derived classes.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               Array<int>&           n_cycle,
                               const Array<IntVect>& ref_ratio,
                               Array<Real>&          dt_min,
                               Array<Real>&          dt_level,
                               Real                  stop_time,
                               int                   post_regrid_flag) = 0;
    //
    // Do an integration step on this level.  Returns maximum safe
    // time step.  This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) = 0;
    //
    // Contains operations to be done after a timestep.  This is a
    // pure virtual function and hence MUST be implemented by derived
    // classes.
    //
    virtual  void post_timestep (int iteration) = 0;
    //
    // Contains operations to be done only after a full coarse
    // timestep.  The default implementation does nothing.
    //
    virtual void postCoarseTimeStep (Real time);
    //
    // Operations to be done after restart.  This is a pure virtual
    // function and hence MUST be implemented by derived classes.
    //
    virtual  void post_restart () = 0;
    //
    // Operations to be done after regridding (like avgDown).
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual  void post_regrid (int lbase,
                               int new_finest) = 0;
    //
    // Operations to be done after initialization.
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual  void post_init (Real stop_time) = 0;
    //
    // Is it ok to continue the calculation?
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual  int okToContinue () = 0;
    //
    // Should I regrid with this level as base level?
    // This test is only evaluated if regrid_int > 0 and 
    // level_count >= regrid_int as well. Defaults to true.
    //
    virtual  int okToRegrid ();
    //
    // Init grid data at problem start-up.
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual void initData () = 0;
    //
    // Set the time levels of state data.
    //
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new);
    //
    // Alloc space for old time data.
    //
    virtual void allocOldData ();
    //
    // Delete old-time data.
    //
    virtual void removeOldData ();
    //
    // Init data on this level from another AmrLevel (during regrid).
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual void init (AmrLevel &old) = 0;
    //
    // Init data on this level after regridding if old AmrLevel
    // did not previously exist. This is a pure virtual function
    // and hence MUST be implemented by derived classes.
    //
    virtual void init () = 0;
    //
    // Reset data to initial time by swapping new and old time data.
    //
    void reset ();
    //
    // Returns this AmrLevel.
    //
    int Level () const { return level; }
    //
    // List of grids at this level.
    //
    const BoxArray& boxArray () const { return grids; }
    const BoxArray& getEdgeBoxArray (int dir) const;
    const BoxArray& getNodalBoxArray () const;
    //
    // Number of grids at this level.
    //
    int numGrids () const { return grids.size(); }
#ifdef USE_PARTICLES
    //
    // List of particle_grids at this level.
    //
    const BoxArray& ParticleBoxArray () const { return particle_grids; }
    //
    // DistributionMap for particle_grids at this level.
    //
    const DistributionMapping& ParticleDistributionMap () const { return particle_dmap; }
    //
    // Number of particle_grids at this level.
    //
    int numParticleGrids () const { return particle_grids.size(); }
#endif
    //
    // Returns the indices defining physical domain.
    //
    const Box& Domain () const { return geom.Domain(); }
    //
    // Timestep n at this level.
    //
    int nStep () const { return parent->levelSteps(level); }
    //
    // Returns the geometry object.
    //
    const Geometry& Geom () const { return geom; }
    //
    const IntVect& fineRatio () const { return fine_ratio; }
    //
    // Returns number of cells on level.
    //
    long countCells () const;

    //
    // Get the area not to tag.
    //
    const BoxArray& getAreaNotToTag();
    const Box& getAreaToTag();
    //
    // Constuct the area not to tag.
    //
    void constructAreaNotToTag();
    //
    // Set the area not to tag.
    //
    void setAreaNotToTag(BoxArray& ba);

    //
    // Error estimation for regridding. This is a pure virtual
    // function and hence MUST be implemented by derived classes.
    //
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
			   int          n_error_buf = 0,
                           int          ngrow = 0) = 0;
    //
    // Interpolate from coarse level to the valid area in dest.
    //
    void FillCoarsePatch (MultiFab& dest,
                          int       dcomp,
                          Real      time,
                          int       state_idx,
                          int       scomp,
                          int       ncomp,
			  int       nghost = 0);
    //
    // Function to set physical boundary conditions.
    //
    virtual void setPhysBoundaryValues (FArrayBox& dest,
                                        int        state_indx,
                                        Real       time,
                                        int        dest_comp,
                                        int        src_comp,
                                        int        num_comp);
    //
    // Returns a MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the MultiFab is built on the appropriately
    // grown BoxArray.
    //
    virtual MultiFab* derive (const std::string& name,
                              Real               time,
                              int                ngrow);
    //
    // This version of derive() fills the dcomp'th component of mf
    // with the derived quantity.
    //
    virtual void derive (const std::string& name,
                         Real               time,
                         MultiFab&          mf,
                         int                dcomp);
    //
    // State data object.
    //
    StateData& get_state_data (int state_indx) { return state[state_indx]; }
    //
    // State data at old time.
    //
    MultiFab& get_old_data (int state_indx) { return state[state_indx].oldData(); }
    //
    // State data at old time.
    //
    const MultiFab& get_old_data (int state_indx) const { return state[state_indx].oldData(); }
    //
    // State data at new time.
    //
    MultiFab& get_new_data (int state_indx) { return state[state_indx].newData(); }
    //
    // State data at new time.
    //
    const MultiFab& get_new_data (int state_indx) const { return state[state_indx].newData(); }
    //
    // Returns list of Descriptors.
    //
    static const DescriptorList& get_desc_lst () { return desc_lst; }
    //
    // Returns list of derived variables.
    //
    static DeriveList& get_derive_lst ();

#ifdef USE_SLABSTAT
    //
    // Returns list of slab stats.
    //
    static SlabStatList& get_slabstat_lst ();
#endif
    //
    // Boundary condition access function.
    //
    Array<int> getBCArray (int State_Type,
			   int gridno,
			   int scomp,
			   int ncomp);
    //
    // Get state data at specified index and time.
    //
    MultiFab& get_data (int  state_indx,
                        Real time);
    //
    // Hack to allow override of (non-fine-fine) fillpatched boundary data
    //
    virtual void set_preferred_boundary_values (MultiFab& S,
                                                int       state_index,
                                                int       scomp,
                                                int       dcomp,
                                                int       ncomp,
                                                Real      time) const;
    // 
    // Called in grid_places after other tagging routines to modify
    // the list of tagged points.  Default implementation does nothing.
    //
    virtual void manual_tags_placement (TagBoxArray&    tags,
                                        Array<IntVect>& bf_lev);
    //
    // Modify list of variables to be plotted
    // 
    virtual void setPlotVariables ();
    //
    // Modify list of variables to be plotted
    // 
    virtual void setSmallPlotVariables ();
    //
    // Estimate the amount of work required to advance Just this level
    // based on the number of cells.
    // This estimate can be overwritten with different methods
    //
    virtual Real estimateWork();
    //
    // Returns one the TimeLevel enums.
    // Asserts that time is between AmrOldTime and AmrNewTime.
    // 
    TimeLevel which_time (int  state_indx,
                          Real time) const;

    //
    // Does the AmrLevel want Amr to write a plotfile now?
    //
    virtual bool writePlotNow ();

    //
    // Does the AmrLevel want Amr to write a small plotfile now?
    //
    virtual bool writeSmallPlotNow ();

#ifdef USE_PARTICLES
    //
    // Add one particle a time to a ParticleContainer owned by this AmrLevel.
    //
    virtual void addOneParticle (int id_in, int cpu_in,
                                 std::vector<double>& xloc, 
                                 std::vector<double>& attributes) {;}

    //
    // Add this function which can be called from the parent 
    //
    virtual void particle_redistribute (int lbase = 0, bool init = false) {;}

    //
    // Get particle data out of the Particle Class and send it back out through 
    //     an Amr* amrptr.
    //
    virtual void GetParticleIDs       (Array<int>&  part_ids ) {;}
    virtual void GetParticleCPU       (Array<int>&  part_cpu ) {;}
    virtual void GetParticleLocations (Array<Real>& part_locs) {;}
    virtual void GetParticleData      (Array<Real>& part_data, int start_comp, int num_comp) {;}

#endif

    static void FillPatch(AmrLevel& amrlevel,
			  MultiFab& leveldata,
			  int       boxGrow,
			  Real      time,
			  int       index,
			  int       scomp,
			  int       ncomp,
	                  int       dcomp=0);
    
    virtual void AddProcsToComp(Amr *aptr, int nSidecarProcs, int prevSidecarProcs,
                                int ioProcNumSCS, int ioProcNumAll, int scsMyId,
			        MPI_Comm scsComm);
    void Check() const;

protected:
    //
    // The constructors -- for derived classes.
    //
    AmrLevel ();

    AmrLevel (Amr&            papa,
              int             lev,
              const Geometry& level_geom,
              const BoxArray& bl,
              Real            time);

    AmrLevel (Amr&            papa,
              int             lev,
              const Geometry& level_geom,
              const BoxArray& bl,
	      const DistributionMapping& dm,
              Real            time);
    //
    // Common code used by all constructors.
    //
    void finishConstructor (); 

    //
    // The Data.
    //
    int                   level;        // AMR level (0 is coarsest).
    Geometry              geom;         // Geom at this level.
    BoxArray              grids;        // Cell-centered locations of grids.
#ifdef USE_PARTICLES
    BoxArray              particle_grids; // Cell-centered locations of grids on which particles live.
    DistributionMapping   particle_dmap; // Distribution map associated with particle_grids.
#endif
    Amr*                  parent;       // Pointer to parent AMR structure.
    IntVect               crse_ratio;   // Refinement ratio to coarser level.
    IntVect               fine_ratio;   // Refinement ratio to finer level.
    static DeriveList     derive_lst;   // List of derived quantities.
    static DescriptorList desc_lst;     // List of state variables.
#ifdef USE_SLABSTAT
    static SlabStatList   slabstat_lst; // List of SlabStats.
#endif
    Array<StateData>      state;        // Array of state data.

    BoxArray              m_AreaNotToTag; //Area which shouldn't be tagged on this level.
    Box                   m_AreaToTag;    //Area which is allowed to be tagged on this level.

private:

    mutable BoxArray      edge_grids[BL_SPACEDIM];  // face-centered grids
    mutable BoxArray      nodal_grids;              // all nodal grids

    //
    // Disallowed.
    //
    AmrLevel (const AmrLevel&);
    AmrLevel& operator = (const AmrLevel&);
};

//
// Forward declaration.
//
class FillPatchIteratorHelper;

class FillPatchIterator
    :
    public MFIter
{
  public:

    friend class AmrLevel;

    FillPatchIterator (AmrLevel& amrlevel,
                       MultiFab& leveldata);

    FillPatchIterator (AmrLevel& amrlevel,
                       MultiFab& leveldata,
                       int       boxGrow,
                       Real      time,
                       int       state_indx,
                       int       scomp,
                       int       ncomp);

    void Initialize (int  boxGrow,
                     Real time,
                     int  state_indx,
                     int  scomp,
                     int  ncomp);

    ~FillPatchIterator ();

    FArrayBox& operator() () { return m_fabs[MFIter::index()]; }

    Box UngrownBox () const { return MFIter::validbox(); }

    MultiFab& get_mf() { return m_fabs; }
    
  private:
    //
    // Disallowed.
    //
    FillPatchIterator ();
    FillPatchIterator (const FillPatchIterator& rhs);
    FillPatchIterator& operator= (const FillPatchIterator& rhs);

    void FillFromLevel0 (Real time, int index, int scomp, int dcomp, int ncomp);
    void FillFromTwoLevels (Real time, int index, int scomp, int dcomp, int ncomp);

    //
    // The data.
    //
    AmrLevel&                         m_amrlevel;
    MultiFab&                         m_leveldata;
    std::vector< std::pair<int,int> > m_range;
    MultiFab                          m_fabs;
    int                               m_ncomp;
};

class FillPatchIteratorHelper
{
public:

    friend class FillPatchIterator;

    FillPatchIteratorHelper (AmrLevel& amrlevel,
                             MultiFab& leveldata);

    FillPatchIteratorHelper (AmrLevel&     amrlevel,
                             MultiFab&     leveldata,
                             int           boxGrow,
                             Real          time,
                             int           state_indx,
                             int           scomp,
                             int           ncomp,
                             Interpolater* mapper);

    void Initialize (int           boxGrow,
                     Real          time,
                     int           state_indx,
                     int           scomp,
                     int           ncomp,
                     Interpolater* mapper);

    ~FillPatchIteratorHelper ();

    void fill (FArrayBox& fab, int dcomp, int idx);

private:
    //
    // Disallowed.
    //
    FillPatchIteratorHelper ();
    FillPatchIteratorHelper (const FillPatchIteratorHelper& rhs);
    FillPatchIteratorHelper& operator= (const FillPatchIteratorHelper& rhs);
    //
    // The data.
    //
    AmrLevel&                  m_amrlevel;
    MultiFab&                  m_leveldata;
    MultiFabCopyDescriptor     m_mfcd;
    Array< Array<MultiFabId> > m_mfid;     // [level][oldnew]
    Interpolater*              m_map;
    std::map<int,Box>          m_ba;
    Real                       m_time;
    int                        m_growsize;
    int                        m_index;
    int                        m_scomp;
    int                        m_ncomp;
    bool                       m_FixUpCorners;

    std::map< int,Array< Array<Box> > >                m_fbox; // [grid][level][validregion]
    std::map< int,Array< Array<Box> > >                m_cbox; // [grid][level][fillablesubbox]
    std::map< int,Array< Array< Array<FillBoxId> > > > m_fbid; // [grid][level][fillablesubbox][oldnew]
};

#endif /*_AmrLevel_H_*/
