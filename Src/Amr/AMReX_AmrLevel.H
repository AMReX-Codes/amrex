
#ifndef AMREX_AmrLevel_H_
#define AMREX_AmrLevel_H_

#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_LayoutData.H>
#include <AMReX_Derive.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Amr.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_StateData.H>
#include <AMReX_VisMF.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBSupport.H>
#include <AMReX_EBInterpolater.H>
#endif

#include <memory>
#include <map>

namespace amrex {

class TagBox;
class TagBoxArray;
template <class T>
class MFGraph;
class RGIter;

/**
* \brief Virtual base class for managing individual levels.
* AmrLevel functions both as a container for state data on a level
* and also manages the advancement of data in time.
*/

class AmrLevel
{
  friend class FillPatchIterator;
  friend class FillPatchIteratorHelper;
  template <class T>
  friend class MFGraph;
  friend class RGIter;
  friend class AsyncFillPatchIterator;

public:
    //! What time are we at?
    enum TimeLevel { AmrOldTime,
                     AmrHalfTime,
                     AmrNewTime,
                     Amr1QtrTime,
                     Amr3QtrTime,
                     AmrOtherTime };
    //! The destructor.
    virtual ~AmrLevel ();
    //! Get the level directory names
    void LevelDirectoryNames (const std::string &dir,
                              std::string &LevelDir,
                              std::string &FullPath);
    //! Create the Level_ directory for checkpoint and plot files
    virtual void CreateLevelDirectory (const std::string &dir);
    /**
    * \brief Set if the Level_ directory was created or to clear the value.
    * CreateLevelDirectory sets levelDirectoryCreated = true
    */
    void SetLevelDirectoryCreated(bool ldc) noexcept { levelDirectoryCreated = ldc; }
    /**
    * \brief A string written as the first item in writePlotFile() at
    * level zero. 
    * It is so we can distinguish between different types of
    * plot files.
    * This default "HyperCLaw-V1.1" is for VisIt software and some of our
    * internal postprocessing routines
    */
    virtual std::string thePlotFileType () const
    {
        static const std::string the_plot_file_type("HyperCLaw-V1.1");
        return the_plot_file_type;
    }
    /**
    * \brief Write plot file stuff to specified directory.  This is a
    * pure virtual function and hence MUST be implemented by
    * derived classes.
    */
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how = VisMF::NFiles);

    //! Do pre-plotfile work to avoid synchronizations while writing the amr hierarchy
    virtual void writePlotFilePre (const std::string& dir,
                                   std::ostream&      os);

    //! Do post-plotfile work to avoid synchronizations while writing the amr hierarchy
    virtual void writePlotFilePost (const std::string& dir,
                                    std::ostream&      os);

    /**
    * \brief Write small plot file stuff to specified directory.  
    * Unlike writePlotFile, this is NOT a pure virtual function 
    * so implementation by derived classes is optional.
    */
    virtual void writeSmallPlotFile (const std::string& dir,
                                     std::ostream&      os,
				     VisMF::How         how = VisMF::NFiles) {};
    //! Write current state to checkpoint file.
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os,
                             VisMF::How         how = VisMF::NFiles,
                             bool               dump_old = true);
    //! Do pre-checkpoint work to avoid synchronizations while writing the amr hierarchy
    virtual void checkPointPre (const std::string& dir,
                                std::ostream&      os);
    //! Do post-checkpoint work to avoid synchronizations while writing the amr hierarchy
    virtual void checkPointPost (const std::string& dir,
                                 std::ostream&      os);
    //! Restart from a checkpoint file.
    virtual void restart (Amr&          papa,
                          std::istream& is,
			  bool          bReadSpecial = false);

    //! Old checkpoint may have different number of states than the new source code.
    virtual void set_state_in_checkpoint (Vector<int>& state_in_checkpoint);

    //! Is name a state variable?
    static bool isStateVariable (const std::string& name,
                                int&               state_indx,
                                int&               ncomp);

    static void FlushFPICache ();
    /**
    * \brief Compute the initial time step.  This is a pure virtual function
    * and hence MUST be implemented by derived classes.
    */
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Vector<int>&           n_cycle,
                                   const Vector<IntVect>& ref_ratio,
                                   Vector<Real>&          dt_level,
                                   Real                  stop_time) = 0;
    /**
    * \brief Compute the next time step.  This is a pure virtual function
    * and hence MUST be implemented by derived classes.
    */
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_min,
                               Vector<Real>&          dt_level,
                               Real                  stop_time,
                               int                   post_regrid_flag) = 0;
    /**
    * \brief Do an integration step on this level.  Returns maximum safe
    * time step.  This is a pure virtual function and hence MUST
    * be implemented by derived classes.
    */
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) = 0;

#ifdef USE_PERILLA
    // For Perilla initialization
    virtual void initPerilla (Real time)=0;
    virtual void finalizePerilla (Real time)=0;
#endif


    /**
    * \brief Contains operations to be done after a timestep.  This is a
    * pure virtual function and hence MUST be implemented by derived
    * classes.
    */
    virtual  void post_timestep (int iteration) = 0;
    /**
    * \brief Contains operations to be done only after a full coarse
    * timestep.  The default implementation does nothing.
    */
    virtual void postCoarseTimeStep (Real time);
    /**
    * \brief Operations to be done after restart.
    */
    virtual void post_restart () {};
    /**
    * \brief Operations to be done after regridding
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual  void post_regrid (int lbase,
                               int new_finest) = 0;
    /**
    * \brief Operations to be done after initialization.
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual  void post_init (Real stop_time) = 0;
    /**
    * \brief Is it ok to continue the calculation?
    */
    virtual  int okToContinue () { return 1; }
    /**
    * \brief Should I regrid with this level as base level?
    * This test is only evaluated if regrid_int > 0 and 
    * level_count >= regrid_int as well. Defaults to true.
    */
    virtual  int okToRegrid ();
    /**
    * \brief Init grid data at problem start-up.
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual void initData () = 0;
    //! Set the time levels of state data.
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new);
    //! Alloc space for old time data.
    virtual void allocOldData ();
    //! Delete old-time data.
    virtual void removeOldData ();
    /**
    * \brief Init data on this level from another AmrLevel (during regrid).
    * This is a pure virtual function and hence MUST be
    * implemented by derived classes.
    */
    virtual void init (AmrLevel &old) = 0;
    /**
    * Init data on this level after regridding if old AmrLevel
    * did not previously exist. This is a pure virtual function
    * and hence MUST be implemented by derived classes.
    */
    virtual void init () = 0;
    //! Reset data to initial time by swapping new and old time data.
    void reset ();
    //! Returns this AmrLevel.
    int Level () const noexcept { return level; }
    //! List of grids at this level.
    const BoxArray& boxArray () const noexcept { return grids; }
    const BoxArray& getEdgeBoxArray (int dir) const noexcept;
    const BoxArray& getNodalBoxArray () const noexcept;
    //
    const DistributionMapping& DistributionMap () const noexcept { return dmap; }
    //
    const FabFactory<FArrayBox>& Factory () const noexcept { return *m_factory; }
    //! Number of grids at this level.
    int numGrids () const noexcept { return grids.size(); }
    //! Number of states at this level.
    int numStates () const noexcept { return state.size(); }
    //! Returns the indices defining physical domain.
    const Box& Domain () const noexcept { return geom.Domain(); }
    //! Timestep n at this level.
    int nStep () const noexcept { return parent->levelSteps(level); }
    //! Returns the geometry object.
    const Geometry& Geom () const noexcept { return geom; }
    //
    const IntVect& fineRatio () const noexcept { return fine_ratio; }
    //! Returns number of cells on level.
    long countCells () const noexcept;

    //! Get the area not to tag.
    const BoxArray& getAreaNotToTag() noexcept;
    const Box& getAreaToTag() noexcept;
    //! Constuct the area not to tag.
    void constructAreaNotToTag();
    //! Set the area not to tag.
    void setAreaNotToTag(BoxArray& ba) noexcept;

    /**
    * \brief Error estimation for regridding. This is a pure virtual
    * function and hence MUST be implemented by derived classes.
    */
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
			   int          n_error_buf = 0,
                           int          ngrow = 0) = 0;
    //! Interpolate from coarse level to the valid area in dest.
    void FillCoarsePatch (MultiFab& dest,
                          int       dcomp,
                          Real      time,
                          int       state_idx,
                          int       scomp,
                          int       ncomp,
			  int       nghost = 0);
    //! Function to set physical boundary conditions.
    virtual void setPhysBoundaryValues (FArrayBox& dest,
                                        int        state_indx,
                                        Real       time,
                                        int        dest_comp,
                                        int        src_comp,
                                        int        num_comp);
    /**
    * \brief Returns a MultiFab containing the derived data for this level.
    * The user is responsible for deleting this pointer when done
    * with it.  If ngrow>0 the MultiFab is built on the appropriately
    * grown BoxArray.
    */
    virtual std::unique_ptr<MultiFab> derive (const std::string& name,
					      Real               time,
					      int                ngrow);
    /**
    * \brief This version of derive() fills the dcomp'th component of mf
    * with the derived quantity.
    */
    virtual void derive (const std::string& name,
                         Real               time,
                         MultiFab&          mf,
                         int                dcomp);
    //! State data object.
    StateData& get_state_data (int state_indx) noexcept { return state[state_indx]; }
    //! State data at old time.
    MultiFab& get_old_data (int state_indx) noexcept { return state[state_indx].oldData(); }
    //! State data at old time.
    const MultiFab& get_old_data (int state_indx) const noexcept { return state[state_indx].oldData(); }
    //! State data at new time.
    MultiFab& get_new_data (int state_indx) noexcept { return state[state_indx].newData(); }
    //! State data at new time.
    const MultiFab& get_new_data (int state_indx) const noexcept { return state[state_indx].newData(); }
    //! Returns list of Descriptors.
    static const DescriptorList& get_desc_lst () noexcept { return desc_lst; }
    //! Returns list of derived variables.
    static DeriveList& get_derive_lst () noexcept;
    //! Returns whether or not we want a post-timestep regrid.
    int postStepRegrid () noexcept { return post_step_regrid; }
    //! Sets a new value for the post-timestep regrid trigger.
    void setPostStepRegrid (int new_val) noexcept { post_step_regrid = new_val; }

    //! Update the distribution maps in StateData based on the size of the map
    void UpdateDistributionMaps ( DistributionMapping& dmap );

    //! Boundary condition access function.
    Vector<int> getBCArray (int State_Type,
			   int gridno,
			   int scomp,
			   int ncomp);
    //! Get state data at specified index and time.
    MultiFab& get_data (int  state_indx, Real time) noexcept;
    //! Hack to allow override of (non-fine-fine) fillpatched boundary data
    virtual void set_preferred_boundary_values (MultiFab& S,
                                                int       state_index,
                                                int       scomp,
                                                int       dcomp,
                                                int       ncomp,
                                                Real      time) const;
    /** 
    * \brief Called in grid_places after other tagging routines to modify
    * the list of tagged points.  Default implementation does nothing.
    */
    virtual void manual_tags_placement (TagBoxArray&    tags,
                                        const Vector<IntVect>& bf_lev);
    //! Modify list of variables to be plotted
    virtual void setPlotVariables ();
    //! Modify list of variables to be plotted
    virtual void setSmallPlotVariables ();
    /**
    * \brief Estimate the amount of work required to advance Just this level
    * based on the number of cells.
    * This estimate can be overwritten with different methods
    */
    virtual Real estimateWork();

    //! Which state data type is for work estimates? -1 means none
    virtual int WorkEstType () { return -1; }

    /**
    * \brief Returns one the TimeLevel enums.
    * Asserts that time is between AmrOldTime and AmrNewTime.
    */ 
    TimeLevel which_time (int  state_indx, Real time) const noexcept;

    //! Does the AmrLevel want Amr to write a plotfile now?
    virtual bool writePlotNow ();

    //! Does the AmrLevel want Amr to write a small plotfile now?
    virtual bool writeSmallPlotNow ();

#ifdef AMREX_PARTICLES
    //! This function can be called from the parent 
    virtual void particle_redistribute (int lbase = 0, bool a_init = false) {;}
#endif

    static void FillPatch (AmrLevel& amrlevel,
                           MultiFab& leveldata,
                           int       boxGrow,
                           Real      time,
                           int       index,
                           int       scomp,
                           int       ncomp,
                           int       dcomp=0);

    static void FillPatchAdd (AmrLevel& amrlevel,
                              MultiFab& leveldata,
                              int       boxGrow,
                              Real      time,
                              int       index,
                              int       scomp,
                             int       ncomp,
                             int       dcomp=0);
    
#ifdef AMREX_USE_EB
    static void SetEBMaxGrowCells (int nbasic, int nvolume, int nfull) noexcept {
        m_eb_basic_grow_cells = nbasic;
        m_eb_volume_grow_cells = nvolume;
        m_eb_full_grow_cells = nfull;
    }
    static int            m_eb_basic_grow_cells;
    static int            m_eb_volume_grow_cells;
    static int            m_eb_full_grow_cells;
    static void SetEBSupportLevel (EBSupport ebs) { m_eb_support_level = ebs; }
    static EBSupport      m_eb_support_level;
#endif

protected:
    //! The constructors -- for derived classes.
    AmrLevel () noexcept;

    AmrLevel (Amr&            papa,
              int             lev,
              const Geometry& level_geom,
              const BoxArray& bl,
	      const DistributionMapping& dm,
              Real            time);

    AmrLevel (const AmrLevel&) = delete;
    AmrLevel& operator = (const AmrLevel&) = delete;

    //! Common code used by all constructors.
    void finishConstructor (); 

    //
    // The Data.
    //
    int                   level;        // AMR level (0 is coarsest).
    Geometry              geom;         // Geom at this level.
    BoxArray              grids;        // Cell-centered locations of grids.
    DistributionMapping   dmap;         // Distribution of grids among processes
    Amr*                  parent;       // Pointer to parent AMR structure.
    IntVect               crse_ratio;   // Refinement ratio to coarser level.
    IntVect               fine_ratio;   // Refinement ratio to finer level.
    static DeriveList     derive_lst;   // List of derived quantities.
    static DescriptorList desc_lst;     // List of state variables.
    Vector<StateData>      state;        // Array of state data.

    BoxArray              m_AreaNotToTag; //Area which shouldn't be tagged on this level.
    Box                   m_AreaToTag;    //Area which is allowed to be tagged on this level.

    int                   post_step_regrid; // Whether or not to do a regrid after the timestep.

    bool                  levelDirectoryCreated;    // for checkpoints and plotfiles

    std::unique_ptr<FabFactory<FArrayBox> > m_factory;

private:

    mutable BoxArray      edge_grids[AMREX_SPACEDIM];  // face-centered grids
    mutable BoxArray      nodal_grids;              // all nodal grids
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

    FArrayBox& operator() () noexcept { return m_fabs[MFIter::index()]; }

    Box UngrownBox () const noexcept { return MFIter::validbox(); }

    MultiFab& get_mf() noexcept { return m_fabs; }

#ifdef USE_PERILLA
    FillPatchIterator (AmrLevel& amrlevel,
                       MultiFab& leveldata,
                       int       boxGrow,
                       Real      time,
                       int       state_indx,
                       int       scomp,
                       int       ncomp,
                       int       f);

    void initFillPatch(int boxGrow, int time, int index, int scomp, int ncomp, int iter);

    void InitializePush (int  boxGrow,
                         Real time,
                         int  state_indx,
                         int  scomp,
                         int  ncomp,
                         int f);

    void InitializePull (int  boxGrow,
                         Real time,
                         int  state_indx,
                         int  scomp,
                         int  ncomp,
                         int f);

    void FillPatchPush (int  boxGrow,
                         Real time,
                         int  state_indx,
                         int  scomp,
                         int  ncomp,
                         int f,
                         unsigned char pushLevel,
                         bool singleT=false);

    void FillPatchPull (int  boxGrow,
                         Real time,
                         int  state_indx,
                         int  scomp,
                         int  ncomp,
                         int f,
                         bool singleT=false);

    void finalizeGraphs()
    {
      //std::cout << "Completing RGs ";

      if(destGraph != NULL)
        {
          //std::cout << destGraph->graphID << " ";
           destGraph->finalizeGraph();
        }
      if(csrcGraph != NULL)
        {
          //std::cout << csrcGraph->graphID << " ";
          csrcGraph->finalizeGraph();
        }
      if(fsrcGraph != NULL)
        {
          //std::cout << fsrcGraph->graphID << " ";
          fsrcGraph->finalizeGraph();
        }
      if(m_rg_crse_patch != NULL)
        {
          //std::cout << m_rg_crse_patch->graphID << " ";
          m_rg_crse_patch->finalizeGraph();
        }

      //std::cout <<" by tg " << tg << std::endl;
    }

    void Reset()
    {
      int tg= perilla::wid();
      //std::cout << "Resetting RGs ";
      if(destGraph != NULL)
        {
          //std::cout << destGraph->graphID << " ";
           destGraph->Reset();
        }
      if(csrcGraph != NULL)
        {
          //std::cout << csrcGraph->graphID << " ";
          csrcGraph->Reset();
        }
      if(fsrcGraph != NULL)
        {
          //std::cout << fsrcGraph->graphID << " ";
          fsrcGraph->Reset();
        }
      if(m_rg_crse_patch != NULL)
        {
          //std::cout << m_rg_crse_patch->graphID << " ";
          m_rg_crse_patch->Reset();
        }
      //std::cout <<" by tg " << tg << std::endl;
    }

    RegionGraph* get_destGraph(){return destGraph;}
    RegionGraph* get_crscGraph(){return csrcGraph;}
    RegionGraph* get_fsrcGraph(){return fsrcGraph;}
#endif

    
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

public:
#ifdef USE_PERILLA
    RegionGraph*                      destGraph;
    RegionGraph*                      csrcGraph;
    RegionGraph*                      fsrcGraph;
    RegionGraph*                      m_rg_crse_patch;
    std::list<RegionGraph*>	      regionList;
    std::list<MultiFab*>	      mfList;
    std::list<StateDataPhysBCFunct*>  stateDataList;



    MultiFab*                         m_mf_crse_patch;
    const FabArrayBase::FPinfo*       m_fpc;
    MultiFab*                         dmf;
    MultiFab*                         dmff;
    Vector<MultiFab*>                 smf;
    Geometry*                         geom;
    StateDataPhysBCFunct*             physbcf;
    bool                              isProperlyNested;
    Vector<MultiFab*>                 smf_crse;
    Vector<Real>                      stime_crse;
    StateDataPhysBCFunct*             physbcf_crse;
    Geometry*                         geom_crse;
    Vector<MultiFab*>                 smf_fine;
    Vector<Real>                      stime_fine;
    StateDataPhysBCFunct*             physbcf_fine;
    Geometry*                         geom_fine;

    Vector<Real>                 stime;
    void FillFromLevel0Push (Real time, int index, int scomp, int dcomp, int ncomp, int f);
    void FillFromLevel0PushOnly (Real time, int index, int scomp, int dcomp, int ncomp, int f, bool singleT);
    void FillFromLevel0Pull (Real time, int index, int scomp, int dcomp, int ncomp, int f, bool singleT);
    void FillFromTwoLevelsPushOnly (Real time, int index, int scomp, int dcomp, int ncomp, int f, unsigned char pushLevel, bool singleT);
    void FillFromTwoLevelsPush (Real time, int index, int scomp, int dcomp, int ncomp, int f, unsigned char pushLevel, bool singleT);
    void FillFromTwoLevelsPull (Real time, int index, int scomp, int dcomp, int ncomp, int f, bool singleT);
    void FillPatchTwoLevelsPush (Amr& amr, MultiFab& mf, Real time,
                                 Vector<MultiFab*>& cmf, Vector<Real>& ct,
                                 Vector<MultiFab*>& fmf, Vector<Real>& ft,
                                 RegionGraph* destGraph, RegionGraph* csrcGraph, RegionGraph* fsrcGraph, int f,
                                 FillPatchIterator* fpIter,
                                 MultiFab *dmf,
                                 MultiFab *dmff,
                                 int scomp, int dcomp, int ncomp,
                                 const Geometry& cgeom, const Geometry& fgeom,
                                 StateDataPhysBCFunct& cbc, StateDataPhysBCFunct& fbc,
                                 const IntVect& ratio,
                                 Interpolater* mapper, const Vector<BCRec>& bcs, unsigned char pushLevel, bool singleT);

    void FillPatchTwoLevelsPull (MultiFab& mf, Real time,
                                 Vector<MultiFab*>& cmf, Vector<Real>& ct,
                                 Vector<MultiFab*>& fmf, Vector<Real>& ft,
                                 RegionGraph* destGraph, RegionGraph* csrcGraph, RegionGraph* fsrcGraph, int f,
                                 FillPatchIterator* fpIter,
                                 int scomp, int dcomp, int ncomp,
                                 const Geometry& cgeom, const Geometry& fgeom,
                                 StateDataPhysBCFunct& cbc, StateDataPhysBCFunct& fbc,
                                 const IntVect& ratio,
                                 Interpolater* mapper, const Vector<BCRec>& bcs, bool singleT);

    void FillPatchSingleLevelPush (Amr& amr, MultiFab& mf, Real time,
                                   Vector<MultiFab*>& smf, Vector<Real>& stime,
                                   RegionGraph* destGraph, RegionGraph* srcGraph, int f,
                                   MultiFab *dmf,
                                   int scomp, int dcomp, int ncomp,
                                   const Geometry& geom, StateDataPhysBCFunct& physbcf, bool singleT);
    void FillPatchSingleLevelPull (MultiFab& mf, Real time,
                                   Vector<MultiFab*>& smf, Vector<Real>& stime,
                                   RegionGraph* destGraph, RegionGraph* srcGraph, int f,
                                   int scomp, int dcomp, int ncomp,
                                   const Geometry& geom, StateDataPhysBCFunct& physbcf, bool singleT);

#endif
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
    Vector< Vector<MultiFabId> > m_mfid;     // [level][oldnew]
    Interpolater*              m_map;
    std::map<int,Box>          m_ba;
    Real                       m_time;
    int                        m_growsize;
    int                        m_index;
    int                        m_scomp;
    int                        m_ncomp;
    bool                       m_FixUpCorners;

    std::map< int,Vector< Vector<Box> > >                m_fbox; // [grid][level][validregion]
    std::map< int,Vector< Vector<Box> > >                m_cbox; // [grid][level][fillablesubbox]
    std::map< int,Vector< Vector< Vector<FillBoxId> > > > m_fbid; // [grid][level][fillablesubbox][oldnew]
};


  //////////////////////////////////////Perilla/////////////////////
#ifdef USE_PERILLA
class AsyncFillPatchIterator
    :
    public MFIter
{
  public:

  friend class AmrLevel;
  friend class RGIter;

    AsyncFillPatchIterator (AmrLevel& amrlevel,
                            MultiFab& leveldata,
                            int       boxGrow,
                            Real      time,
                            int       state_indx,
                            int       scomp,
                            int       ncomp,
                            int iter);

    void initFillPatch(int boxGrow,
                       Real time,
                       int index,
                       int scomp,
                       int ncomp,
                       int iter);

    static void  initialSend(amrex::Vector<amrex::AsyncFillPatchIterator*> afpi,
                             amrex::Vector<amrex::AsyncFillPatchIterator*> upper_afpi,
                             int  boxGrow,
                             Real time,
                             int  state_indx,
                             int  scomp,
                             int  ncomp,
                             int  iter);

    void PushOnly (int  boxGrow,
                   Real time,
                   int  state_indx,
                   int  scomp,
                   int  ncomp,
                   int f,
                   unsigned char pushLevel,
                   bool singleT=false);

    void SendIntraLevel (RGIter& rgi,
                         int  boxGrow,
                         Real time,
                         int  state_indx,
                         int  scomp,
                         int  ncomp,
                         int  iter,
                         int f,
                         bool singleT=false);

    void SendIntraLevel (RGIter* rgi,
                         int  boxGrow,
                         Real time,
                         int  state_indx,
                         int  scomp,
                         int  ncomp,
                         int  iter,
                         int f,
                         bool singleT=false);

    void SendInterLevel (RGIter& rgi,
                         int  boxGrow,
                         Real time,
                         int  state_indx,
                         int  scomp,
                         int  ncomp,
                         int  iter,
                         int f,
                         bool singleT=false);

    void SendInterLevel (RGIter* rgi,
                         int  boxGrow,
                         Real time,
                         int  state_indx,
                         int  scomp,
                         int  ncomp,
                         int  iter,
                         int f,
                         bool singleT=false);

    void Receive (RGIter& rgi,
                   int  boxGrow,
                   Real time,
                   int  state_indx,
                   int  scomp,
                   int  ncomp,
                   int f,
                   bool singleT=false);

    void Receive (RGIter* rgi,
                   int  boxGrow,
                   Real time,
                   int  state_indx,
                   int  scomp,
                   int  ncomp,
                   int f,
                   bool singleT=false);

    void Receive (RGIter& rgi,
                  MultiFab& dest,
                  int  boxGrow,
                  Real time,
                  int  state_indx,
                  int  scomp,
                  int  ncomp,
                  int f,
                  bool singleT=false);

    void Receive (RGIter* rgi,
                  MultiFab& dest,
                  int  boxGrow,
                  const Real time,
                  int  state_indx,
                  int  scomp,
                  int  ncomp,
                  int f,
                  bool singleT=false);

    void PullOnly (int  boxGrow,
                   Real time,
                   int  state_indx,
                   int  scomp,
                   int  ncomp,
                   int f,
                   bool singleT=false);

    void PullOnly (MultiFab& dest,
                   int  boxGrow,
                   Real time,
                   int  state_indx,
                   int  scomp,
                   int  ncomp,
                   int f,
                   bool singleT=false);

    void FillPatchTwoLevelsPush (Amr& amr, MultiFab& mf, Real time,
                                 Vector<MultiFab*>& cmf, Vector<Real>& ct,
                                 Vector<MultiFab*>& fmf, Vector<Real>& ft,
                                 RegionGraph* destGraph, RegionGraph* csrcGraph, RegionGraph* fsrcGraph, int f,
                                 AsyncFillPatchIterator* fpIter,
                                 MultiFab *dmf,
                                 MultiFab *dmff,
                                 int scomp, int dcomp, int ncomp,
                                 const Geometry& cgeom, const Geometry& fgeom,
                                 StateDataPhysBCFunct& cbc, StateDataPhysBCFunct& fbc,
                                 const IntVect& ratio,
                                 Interpolater* mapper, const Vector<BCRec>& bcs, unsigned char pushLevel, bool singleT);

    void FillPatchTwoLevelsPull (MultiFab& mf, Real time,
                                 Vector<MultiFab*>& cmf, Vector<Real>& ct,
                                 Vector<MultiFab*>& fmf, Vector<Real>& ft,
                                 RegionGraph* destGraph, RegionGraph* csrcGraph, RegionGraph* fsrcGraph, int f,
                                 AsyncFillPatchIterator* fpIter,
                                 int scomp, int dcomp, int ncomp,
                                 const Geometry& cgeom, const Geometry& fgeom,
                                 StateDataPhysBCFunct& cbc, StateDataPhysBCFunct& fbc,
                                 const IntVect& ratio,
                                 Interpolater* mapper, const Vector<BCRec>& bcs, bool singleT);

    void FillPatchSingleLevelPush (Amr& amr, MultiFab& mf, Real time,
                                   Vector<MultiFab*>& smf, const Vector<Real>& stime,
                                   RegionGraph* destGraph, RegionGraph* srcGraph, int f,
                                   MultiFab *dmf,
                                   int scomp, int dcomp, int ncomp,
                                   const Geometry& geom, StateDataPhysBCFunct& physbcf, bool singleT);
    void FillPatchSingleLevelPull (MultiFab& mf, Real time,
                                   Vector<MultiFab*>& smf, const Vector<Real>& stime,
                                   RegionGraph* destGraph, RegionGraph* srcGraph, int f,
                                   int scomp, int dcomp, int ncomp,
                                   const Geometry& geom, StateDataPhysBCFunct& physbcf, bool singleT);

    void FillFromTwoLevelsPush (Real time,
                                int index,
                                int scomp,
                                int dcomp,
                                int ncomp,
                                int f,
                                unsigned char pushLevel,
                                bool singleT);
    void FillFromTwoLevelsPull (Real time,
                                int index,
                                int scomp,
                                int dcomp,
                                int ncomp,
                                int f,
                                bool singleT);

    void FillFromTwoLevelsPull (MultiFab& dest,
                                Real time,
                                int index,
                                int scomp,
                                int dcomp,
                                int ncomp,
                                int f,
                                bool singleT);

    ~AsyncFillPatchIterator ();

    FArrayBox& operator() () { return m_fabs[MFIter::index()]; }

    Box UngrownBox () const { return MFIter::validbox(); }

    MultiFab& get_mf() { return m_fabs; }

  //  protected:
    //
    // Disallowed.
    //
    AsyncFillPatchIterator ();
    AsyncFillPatchIterator (const AsyncFillPatchIterator& rhs);
    AsyncFillPatchIterator& operator= (const AsyncFillPatchIterator& rhs);

    //
    // The data.
    //
    AmrLevel&                         m_amrlevel;
    MultiFab&                         m_leveldata;
    std::vector< std::pair<int,int> > m_range;
    MultiFab                          m_fabs;
    int                               m_ncomp;

public:
    bool                              isProperlyNested;

    amrex::Vector<MultiFab*>                  smf;
    amrex::Vector<Real>                       stime;
    StateDataPhysBCFunct*             physbcf;
    Geometry*                         geom;


    amrex::Vector<MultiFab*>                  smf_crse;
    amrex::Vector<Real>                       stime_crse;
    StateDataPhysBCFunct*             physbcf_crse;
    Geometry*                         geom_crse;

    amrex::Vector<MultiFab*>                  smf_fine;
    amrex::Vector<Real>                       stime_fine;
    StateDataPhysBCFunct*             physbcf_fine;
    Geometry*                         geom_fine;


    RegionGraph*                      destGraph;
    RegionGraph*                      csrcGraph;
    RegionGraph*                      fsrcGraph;

    MultiFab*                         m_mf_crse_patch;
    RegionGraph*                      m_rg_crse_patch;
    const FabArrayBase::FPinfo*       m_fpc;

  //PArray<MultiFab>                  raii;
    MultiFab*                         dmf;
    MultiFab*                         dmff;
    std::list<RegionGraph*>           regionList;
    std::list<MultiFab*>              mfList;
    std::list<StateDataPhysBCFunct*>  stateDataList;



    void completeRegionGraphs()
    { 
      //std::cout << "Completing RGs ";
      
      if(destGraph != NULL)
        { 
          //std::cout << destGraph->graphID << " ";
           destGraph->finalizeRegionGraph();
        }
      if(csrcGraph != NULL)
        { 
          //std::cout << csrcGraph->graphID << " ";
          csrcGraph->finalizeRegionGraph();
        }
      if(fsrcGraph != NULL)
        { 
          //std::cout << fsrcGraph->graphID << " ";
          fsrcGraph->finalizeRegionGraph();
        }
      if(m_rg_crse_patch != NULL)
        { 
          //std::cout << m_rg_crse_patch->graphID << " ";
          m_rg_crse_patch->finalizeRegionGraph();
        }
      //std::cout <<" by tg " << tg << std::endl;
    }

    void Reset()
    {
      //std::cout << "Resetting RGs ";
      if(destGraph != NULL)
        {
          //std::cout << destGraph->graphID << " ";
           destGraph->Reset();
        }
      if(csrcGraph != NULL)
        {
          //std::cout << csrcGraph->graphID << " ";
          csrcGraph->Reset();
        }
      if(fsrcGraph != NULL)
        {
          //std::cout << fsrcGraph->graphID << " ";
          fsrcGraph->Reset();
        }
      if(m_rg_crse_patch != NULL)
        {
          //std::cout << m_rg_crse_patch->graphID << " ";
          m_rg_crse_patch->Reset();
        }
      //std::cout <<" by tg " << tg << std::endl;
    }

  // Variables for optimization calls of two level push/pulll

};
#endif


}

#endif /*_AmrLevel_H_*/
