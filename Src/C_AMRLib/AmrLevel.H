#ifndef _AmrLevel_H_
#define _AmrLevel_H_ 
//
// $Id: AmrLevel.H,v 1.56 2005-01-10 22:05:20 car Exp $
//
#include <REAL.H>
#include <Box.H>
#include <RealBox.H>
#include <Geometry.H>
#include <FArrayBox.H>
#include <BoxArray.H>
#include <Array.H>
#include <FabArray.H>
#include <MultiFab.H>
#include <Derive.H>
#include <BCRec.H>
#include <Interpolater.H>
#include <Amr.H>
#include <SlabStat.H>
#include <StateDescriptor.H>
#include <StateData.H>
#include <VisMF.H>

class TagBox;
class TagBoxArray;

//
//@Man:
//@Memo: Virtual base class for managing individual levels.
/*@Doc:

  AmrLevel functions both as a container for state data on a level
  and also manages the advancement of data in time.
*/

class AmrLevel
{
  friend class FillPatchIterator;
  friend class FillPatchIteratorHelper;

public:
    //
    //@ManDoc: What time are we at?
    //
    enum TimeLevel { AmrOldTime,
                     AmrHalfTime,
                     AmrNewTime,
                     Amr1QtrTime,
                     Amr3QtrTime,
                     AmrOtherTime };
    //
    //@ManDoc: The destructor.
    //
    virtual ~AmrLevel ();

    /*@ManDoc: A string written as the first item in writePlotFile() at
               level zero. This MUST be defined by each derived class.
               It is so we can distinguish between different types of
               plot files.  This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual std::string thePlotFileType () const = 0;

    /*@ManDoc: Write plot file stuff to specified directory.  This is a
               pure virtual function and hence MUST be implemented by
               derived classes.
    */
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how = VisMF::OneFilePerCPU) = 0;
    //
    //@ManDoc: Write current state to checkpoint file.
    //
    virtual void checkPoint (const std::string& dir,
                             std::ostream&  os,
                             VisMF::How     how = VisMF::OneFilePerCPU);
    //
    //@ManDoc: Restart from a checkpoint file.
    //
    virtual void restart (Amr&          papa,
                          std::istream& is,
			  bool          bReadSpecial = false);
    //
    //@ManDoc: Is name a state variable?
    //
    static bool isStateVariable (const std::string& name,
                                int&               state_indx,
                                int&               ncomp);

    /*@ManDoc: Compute the initial time step.  This is a pure virtual function
               and hence MUST be implemented by derived classes.
    */
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Array<int>&           n_cycle,
                                   const Array<IntVect>& ref_ratio,
                                   Array<Real>&          dt_level,
                                   Real                  stop_time) = 0;

    /*@ManDoc: Compute the next time step.  This is a pure virtual function
               and hence MUST be implemented by derived classes.
    */
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               Array<int>&           n_cycle,
                               const Array<IntVect>& ref_ratio,
                               Array<Real>&          dt_min,
                               Array<Real>&          dt_level,
                               Real                  stop_time) = 0;

    /*@ManDoc: Do an integration step on this level.  Returns maximum safe
               time step.  This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) = 0;

    /*@ManDoc: Contains operations to be done after a timestep.  This is a
               pure virtual function and hence MUST be implemented by derived
               classes.
    */
    virtual  void post_timestep (int iteration) = 0;

    /*@ManDoc: Contains operations to be done only after a full coarse
               timestep.  The default implementation does nothing.
    */
    virtual void postCoarseTimeStep (Real time);

    /*@ManDoc: Operations to be done after restart.  This is a pure virtual
               function and hence MUST be implemented by derived classes.
    */
    virtual  void post_restart () = 0;

    /*@ManDoc: Operations to be done after regridding (like avgDown).
               This is a pure virtual function and hence MUST be
               implemented by derived classes.
    */
    virtual  void post_regrid (int lbase,
                               int new_finest) = 0;

    /*@ManDoc: Operations to be done after initialization.
               This is a pure virtual function and hence MUST be
               implemented by derived classes.
    */
    virtual  void post_init (Real stop_time) = 0;

    /*@ManDoc: Is it ok to continue the calculation?
               This is a pure virtual function and hence MUST be
               implemented by derived classes.
    */
    virtual  int okToContinue () = 0;

    /*@ManDoc: Should I regrid with this level as base level?
               This test is only evaluated when level\_count >= regrid\_int
               as well. Defaults to true.
    */
    virtual  int okToRegrid ();

    /*@ManDoc: Init grid data at problem start-up.
               This is a pure virtual function and hence MUST be
               implemented by derived classes.
    */
    virtual void initData () = 0;
    //
    //@ManDoc: Set the time levels of state data.
    //
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new);
    //
    //@ManDoc: Alloc space for old time data.
    //
    virtual void allocOldData ();
    //
    //@ManDoc: Delete old-time data.
    //
    virtual void removeOldData ();

    /*@ManDoc: Init data on this level from another AmrLevel (during regrid).
               This is a pure virtual function and hence MUST be
               implemented by derived classes.
    */
    virtual void init (AmrLevel &old) = 0;

    /*@ManDoc: Init data on this level after regridding if old AmrLevel
               did not previously exist. This is a pure virtual function
               and hence MUST be implemented by derived classes.
    */
    virtual void init () = 0;
    //
    //@ManDoc: Reset data to initial time by swapping new and old time data.
    //
    void reset ();
    //
    //@ManDoc: Returns this AmrLevel.
    //
    int Level () const;
    //
    //@ManDoc: List of grids at this level.
    //
    const BoxArray& boxArray () const;
    //
    //@ManDoc: Number of grids at this level.
    //
    int numGrids () const;
    //
    //@ManDoc: Physical locations of grids at this level.
    //
    const Array<RealBox>& gridLocations () const;
    //
    //@ManDoc: Returns the indices defining physical domain.
    //
    const Box& Domain () const;
    //
    //@ManDoc: Timestep n at this level.
    //
    int nStep () const;
    //
    //@ManDoc: Returns the geometry object.
    //
    const Geometry& Geom () const;
    //
    //@ManDoc: Returns number of cells on level.
    //
    long countCells () const;

    /*@ManDoc: Error estimation for regridding. This is a pure virtual
               function and hence MUST be implemented by derived classes.
    */
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
			   int          n_error_buf = 0,
                           int          ngrow = 0) = 0;
    //
    //@ManDoc: Interpolate from coarse level to the valid area in dest.
    //
    void FillCoarsePatch (MultiFab& dest,
                          int       dcomp,
                          Real      time,
                          int       state_idx,
                          int       scomp,
                          int       ncomp);
    //
    //@ManDoc: Function to set physical boundary conditions.
    //
    void setPhysBoundaryValues (int state_indx,
                                int comp,
                                int ncomp,
                                int do_new = 1);
    //
    //@ManDoc: Another function to set physical boundary conditions.
    //
    void setPhysBoundaryValues (int  state_indx,
                                int  comp,
                                int  ncomp,
                                Real time);

    /*@ManDoc: Returns a MultiFab containing the derived data for this level.
               The user is responsible for deleting this pointer when done
               with it.  If ngrow>0 the MultiFab is built on the appropriately
               grown BoxArray.
    */
    MultiFab* derive (const std::string& name,
                      Real               time,
                      int                ngrow);

    /*@ManDoc: This version of derive() fills the dcomp'th component of mf
               with the derived quantity.
    */
    void derive (const std::string& name,
                 Real               time,
                 MultiFab&          mf,
                 int                dcomp);
    //
    //@ManDoc: State data object.
    //
    StateData& get_state_data (int state_indx);
    //
    //@ManDoc: State data at old time.
    //
    MultiFab& get_old_data (int state_indx);
    //
    //@ManDoc: State data at old time.
    //
    const MultiFab& get_old_data (int state_indx) const;
    //
    //@ManDoc: State data at new time.
    //
    MultiFab& get_new_data (int state_indx);
    //
    //@ManDoc: State data at new time.
    //
    const MultiFab& get_new_data (int state_indx) const;
    //
    //@ManDoc: Returns list of Descriptors.
    //
    static const DescriptorList& get_desc_lst ();
    //
    //@ManDoc: Returns list of derived variables.
    //
    static DeriveList& get_derive_lst ();
    //
    //@ManDoc: Returns list of slab stats.
    //
    static SlabStatList& get_slabstat_lst ();
    //
    //@ManDoc: Boundary condition access function.
    //
    Array<int> getBCArray (int State_Type,
                           int gridno,
                           int scomp,
                           int ncomp);
    //
    //@ManDoc: Get state data at specified index and time.
    //
    MultiFab& get_data (int  state_indx,
                        Real time);
    //
    //@ManDoc: Hack to allow override of (non-fine-fine) fillpatched boundary data
    //
    virtual void set_preferred_boundary_values (MultiFab& S,
                                                int       state_index,
                                                int       scomp,
                                                int       dcomp,
                                                int       ncomp,
                                                Real      time) const;
    // 
    //@ManDoc: Called in grid\_places after other tagging routines to modify
    //         the list of tagged points.  Default implementation does nothing.
    //
    virtual void manual_tags_placement (TagBoxArray&    tags,
                                        Array<IntVect>& bf_lev);
    //
    //@ManDoc: Modify list of variables to be plotted
    // 
    virtual void setPlotVariables ();
    //
    //@ManDoc: Returns one the TimeLevel enums.
    //         Asserts that time is between AmrOldTime and AmrNewTime.
    // 
    TimeLevel which_time (int  state_indx,
                          Real time) const;
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
    //
    // Common code used by all constructors.
    //
    void finishConstructor (); 
    //
    // The Data.
    //
    int level;                        // AMR level (0 is coarsest).
    Geometry geom;                    // Geom at this level.
    BoxArray grids;                   // Cell-centered locations of grids.
    Array<RealBox> grid_loc;          // Physical locations of each grid.
    Amr* parent;                      // Pointer to parent AMR structure.
    IntVect crse_ratio;               // Refinement ratio to coarser level.
    IntVect fine_ratio;               // Refinement ratio to finer level.
    static DeriveList derive_lst;     // List of derived quantities.
    static DescriptorList desc_lst;   // List of state variables.
    static SlabStatList slabstat_lst; // List of SlabStats.
    Array<StateData> state;           // Array of state data.

private:
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

    FArrayBox& operator() () { return m_fabs[index()]; }

    bool isValid ();

    void operator++ ();

    const Box& UngrownBox () const { return MFIter::validbox(); }

  private:
    //
    // Disallowed.
    //
    FillPatchIterator ();
    FillPatchIterator (const FillPatchIterator& rhs);
    FillPatchIterator& operator= (const FillPatchIterator& rhs);
    //
    // The data.
    //
    AmrLevel&                         m_amrlevel;
    MultiFab&                         m_leveldata;
    std::vector< std::pair<int,int> > m_range;
    PArray<FillPatchIteratorHelper>   m_fph;
    MultiFab                          m_fabs;
    int                               m_ncomp;
};

class FillPatchIteratorHelper
    :
    public MFIter
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

    bool isValid ();

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
    AmrLevel&                    m_amrlevel;
    MultiFab&                    m_leveldata;
    MultiFabCopyDescriptor       m_mfcd;
    Array<Array<MultiFabId> >    m_mfid;     // [level][oldnew]
    Interpolater*                m_map;
    Array<Array<Array<Box> > >   m_finebox;  // [grid][level][validregion]
    Array<Array<Array<Box> > >   m_crsebox;  // [grid][level][fillablesubbox]
    Array<Array<Array<Array<FillBoxId> > > > m_fbid;// [grid][level][fillablesubbox][oldnew]
    BoxArray                     m_ba;
    Real                         m_time;
    int                          m_growsize;
    int                          m_index;
    int                          m_scomp;
    int                          m_ncomp;
    bool                         m_init;
    bool                         m_FixUpCorners;
};

#endif /*_AmrLevel_H_*/
