//BL_COPYRIGHT_NOTICE

#ifndef _AmrLevel_H_
#define _AmrLevel_H_ 

//
// $Id: AmrLevel.H,v 1.17 1998-06-27 21:34:42 lijewski Exp $
//

#include <REAL.H>
#include <Box.H>
#include <RealBox.H>
#include <Geometry.H>
#include <FArrayBox.H>
#include <BoxList.H>
#include <TagBox.H>
#include <BoxArray.H>
#include <Array.H>
#include <FabArray.H>
#include <MultiFab.H>
#include <Derive.H>
#include <BCRec.H>
#include <Interpolater.H>
#include <Amr.H>
#include <StateDescriptor.H>
#include <StateData.H>
#include <VisMF.H>

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

public:
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
    virtual aString thePlotFileType () const = 0;

    /*@ManDoc: Write plot file stuff to specified directory.  This is a
               pure virtual function and hence MUST be implemented by
               derived classes.
    */
    virtual void writePlotFile (const aString& dir,
                                ostream&       os,
                                VisMF::How     how = VisMF::OneFilePerCPU) = 0;
    //
    //@ManDoc: Write current state to checkpoint file.
    //
    virtual void checkPoint (const aString& dir,
                             ostream&       os,
                             VisMF::How     how = VisMF::OneFilePerCPU);
    //
    //@ManDoc: Restart from a checkpoint file.
    //
    virtual void restart (Amr&     papa,
                          istream& is);
    //
    //@ManDoc: Is name a state variable?
    //
    static int isStateVariable (const aString& name,
                                int&           state_indx,
                                int&           n_comp);

    /*@ManDoc: Compute the initial time step.  This is a pure virtual function
               and hence MUST be implemented by derived classes.
    */
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Array<int>&           n_cycle,
                                   const Array<IntVect>& ref_ratio,
                                   Array<Real>&          dt_level,
                                   REAL                  stop_time) = 0;

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
    virtual  void post_timestep () = 0;

    /*@ManDoc: Contains operations to be done only after a full coarse
               timestep.  The default implementation does nothing.
    */
    virtual void postCoarseTimeStep (Real time) {}

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
    virtual  void post_init (REAL stop_time) = 0;

    /*@ManDoc: Is it ok to continue the calculation?
               This is a pure virtual function and hence MUST be
               implemented by derived classes.
    */
    virtual  int okToContinue () = 0;

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
    int Level () { return level; }
    //
    //@ManDoc: List of grids at this level.
    //
    const BoxArray& boxArray () const { return grids; }
    //
    //@ManDoc: Number of grids at this level.
    //
    int numGrids () { return grids.length(); }
    //
    //@ManDoc: Physical locations of grids at this level.
    //
    const Array<RealBox>& gridLocations () const { return grid_loc; }
    //
    //@ManDoc: Returns the indices defining physical domain.
    //
    const Box& Domain () const { return geom.Domain(); }
    //
    //@ManDoc: Timestep n at this level.
    //
    int nStep () { return parent->levelSteps(level); }
    //
    //@ManDoc: Returns the geometry object.
    //
    const Geometry& Geom () const { return geom; }
    //
    //@ManDoc: Returns number of cells on level.
    //
    long countCells ();

    /*@ManDoc: Error estimation for regridding. This is a pure virtual
               function and hence MUST be implemented by derived classes.
    */
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time) = 0;
    //
    //@ManDoc: Interpolate up from coarse level.
    //
    void FillCoarsePatch (MultiFab&     dest,
                          int           dest_comp,
                          Real          time,
                          int           state_indx,
                          int           src_comp,
                          int           ncomp,
                          Interpolater* mapper = 0);
    //
    //@ManDoc: Function to set physical boundary conditions.
    //
    void setPhysBoundaryValues (int state_indx,
                                int comp,
                                int ncomp,
                                int do_new = 1)
    {
        state[state_indx].FillBoundary(geom.CellSize(),geom.ProbDomain(),
                                       comp,ncomp,do_new);
    }
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
    MultiFab* derive (const aString& name,
                      Real           time,
                      int            ngrow);
    //
    //@ManDoc: State data object.
    //
    StateData& get_state_data (int state_indx) { return state[state_indx]; }
    //
    //@ManDoc: State data at old time.
    //
    MultiFab& get_old_data (int state_indx)
    {
        return state[state_indx].oldData();
    }
    //
    //@ManDoc: State data at new time.
    //
    MultiFab& get_new_data (int state_indx)
    {
        return state[state_indx].newData();
    }
    //
    //@ManDoc: Returns list of Descriptors.
    //
    const DescriptorList& get_desc_lst () const
    {
        return desc_lst;
    }
    //
    //@ManDoc: Boundary condition access function.
    //
    Array<int> getBCArray (int State_Type,
                           int gridno,
                           int strt_comp,
                           int num_comp);
    //
    //@ManDoc: Get state data at specified index and time.
    //
    MultiFab& get_data (int  state_indx,
                        Real time);

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
    int level;                      // AMR level (0 is coarsest).
    Geometry geom;                  // Geom at this level.
    BoxArray grids;                 // Cell-centered locations of grids.
    Array<RealBox> grid_loc;        // Physical locations of each grid.
    Amr* parent;                    // Pointer to parent AMR structure.
    IntVect crse_ratio;             // Refinement ratio to coarser level.
    IntVect fine_ratio;             // Refinement ratio to finer level.
    static DeriveList derive_lst;   // List of derived quantities.
    static DescriptorList desc_lst; // List of state variables.
    Array<StateData> state;         // Array of state data.

private:
    //
    // Disallowed.
    //
    AmrLevel(const AmrLevel&);
    AmrLevel& operator = (const AmrLevel&);
};

class FillPatchIterator
    :
    public MultiFabIterator
{
  public:

    enum FillType { FillLocally, FillRemotely, FillFromBoundary, Unfillable };
    //
    // Constructors.
    //
    FillPatchIterator (AmrLevel& amrlevel,
                       MultiFab& leveldata);

    FillPatchIterator (AmrLevel&     amrlevel,
                       MultiFab&     leveldata,
                       int           boxGrow,
                       int           dest_comp,
                       Real          time,
                       int           state_indx,
                       int           src_comp,
                       int           ncomp,
                       Interpolater* mapper = 0);

    void Initialize (int           boxGrow,
                     int           dest_comp,
                     Real          time,
                     int           state_indx,
                     int           src_comp,
                     int           ncomp,
                     Interpolater* mapper = 0);

    ~FillPatchIterator ();

    const FArrayBox& operator() () const { return m_fab; }

    FArrayBox& operator() ()  { return m_fab; }

    bool isValid (bool bDoSync = true);

    const Box& UngrownBox() const { return validbox(); }

  private:
    //
    // Disallowed.
    //
    FillPatchIterator ();
    FillPatchIterator (const FillPatchIterator& rhs);
    FillPatchIterator& operator= (const FillPatchIterator& rhs);

  private:
    //
    // The data.
    //
    AmrLevel&                    m_amrlevel;
    MultiFab&                    m_leveldata;
    MultiFabCopyDescriptor       m_mfcd;
    Array<Array<MultiFabId> >    m_mfid;     // [level][oldnew]
    Array<IntVect>               m_ratio;    // [level]
    Array<Interpolater*>         m_map;      // [level]
    Array<Array<Array<Box> > >   m_finebox;  // [grid][level][fillablesubbox]
    Array<Array<Array<Array<FillBoxId> > > > m_fbid;// [grid][level][fillablesubbox][oldnew]
    BoxArray                     m_ba;
    FArrayBox                    m_fab;
    Real                         m_time;
    int                          m_growsize;
    int                          m_stateindex;
    int                          m_scomp;
    int                          m_dcomp;
    int                          m_ncomp;
    bool                         m_init;
};

#endif /*_AmrLevel_H_*/
