//BL_COPYRIGHT_NOTICE

#ifndef _AmrLevel_H_
#define _AmrLevel_H_ 

//
// $Id: AmrLevel.H,v 1.9 1997-12-11 23:27:47 lijewski Exp $
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
// Virtual base class for managing individual levels.
//
// AmrLevel functions both as a container for state data on a level
// and also manages the advancement of data in time.
//
class AmrLevel
{
  friend class FillPatchIterator;

public:

    virtual ~AmrLevel();

#ifdef BL_PARALLEL_IO
    //
    // A string written as the first item in writePlotFile() at level zero.
    // This MUST be defined by each derived class.
    // It is so we can distinguish between different types of plot files.
    //
    virtual aString thePlotFileType () const = 0;
    //
    // Write plot file stuff to specified directory.
    //
    virtual void writePlotFile (const aString& dir,
                                ostream&       os,
                                VisMF::How     how = VisMF::OneFilePerCPU) = 0;
    //
    // Write current state to checkpoint file.
    //
    virtual void checkPoint (const aString& dir,
                             ostream&       os,
                             VisMF::How     how = VisMF::OneFilePerCPU);
#else
    //
    // Write plot file.
    //
    virtual void writePlotFile (ostream &os) = 0;
    //
    // Write current state to checkpoint file.
    //
    virtual void checkPoint (ostream& os);
#endif
    //
    // Restart from a checkpoint file.
    //
    virtual void restart (Amr&     papa,
                          istream& is);
    //
    // Is name a state variable?
    //
    int isStateVariable (const aString& name,
                         int&           state_indx,
                         int&           n_comp);
    //
    // Compute the initial time step.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Array<int>&           n_cycle,
                                   const Array<IntVect>& ref_ratio,
                                   Array<Real>&          dt_level) = 0;
    //
    // Compute the next time step.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               Array<int>&           n_cycle,
                               const Array<IntVect>& ref_ratio,
                               Array<Real>&          dt_min,
                               Array<Real>&          dt_level,
                               Real                  stop_time) = 0;
    //
    // Do an integration step on this level - returns maximum safe time step.
    //
    virtual  Real advance (Real time,
                           Real dt,
                           int  iteration,
                           int  ncycle) = 0;
    //
    // Contains operations to be done after a timestep.
    //
    virtual  void post_timestep () = 0;
    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep (Real time) {}
    //
    // Operations to be done after restart.
    //
    virtual  void post_restart () = 0;
    //
    // Operations to be done after regridding (like avgDown).
    //
    virtual  void post_regrid (int lbase,
                               int new_finest) = 0;
    //
    // Operations to be done after initialization.
    //
    virtual  void post_init () = 0;
    //
    // Is it ok to continue the calculation?
    //
    virtual  int okToContinue () = 0;
    //
    // Init grid data at problem start-up.
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
    //
    virtual void init (AmrLevel &old) = 0;
    //
    // Init data on this level after regridding if old AmrLevel
    // did not previously exist.
    //
    virtual void init () = 0;
    //
    // Reset data to initial time by swapping new and old time data.
    //
    void reset ();
    //
    // Returns this AmrLevel.
    //
    int Level () { return level; }
    //
    // List of grids at this level.
    //
    const BoxArray& boxArray () const { return grids; }
    //
    // Number of grids at this level.
    //
    int numGrids () { return grids.length(); }
    //
    // Physical locations of grids at this level.
    //
    const Array<RealBox>& gridLocations () { return grid_loc; }
    //
    // Returns the indices defining physical domain.
    //
    const Box& Domain () { return geom.Domain(); }
    //
    // timestep n at this level.
    //
    int nStep () { return parent->levelSteps(level); }
    //
    // Returns a geometry object.
    //
    const Geometry& Geom () { return geom; }
    //
    // Returns number of cells on level.
    //
    long countCells ();
    //
    // Error estimation for regridding.
    //
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time) = 0;
    //
    // General fill function.
    //
    void FillPatch (FArrayBox&   dest,
                    int          dest_comp,
                    Real         time,
                    int          state_indx,
                    int          src_comp,
                    int          ncomp,
                   Interpolater* mapper = 0);
    //
    // Interpolate up from coarse level.
    //
    void FillCoarsePatch (FArrayBox&    dest,
                          int           dest_comp,
                          Real          time,
                          int           state_indx,
                          int           src_comp,
                          int           ncomp,
                          Interpolater* mapper = 0);

    void FillCoarsePatch (MultiFab&     dest,
                          int           dest_comp,
                          Real          time,
                          int           state_indx,
                          int           src_comp,
                          int           ncomp,
                          Interpolater* mapper = 0);
    //
    // Function to set physical boundary conditions.
    //
    void setPhysBoundaryValues (int state_indx,
                                int comp,
                                int ncomp,
                                int do_new = 1)
    {
        state[state_indx].FillBoundary(geom.CellSize(),geom.ProbDomain(),
                                       comp,ncomp,do_new);
    }

    void setPhysBoundaryValues (int  state_indx,
                                int  comp,
                                int  ncomp,
                                Real time);
    //
    // Returns a pointer to a fab containing the derived data.
    //
    // The user is responsible for deleting this pointer when done with it.
    //
    FArrayBox* derive (const Box&     bx,
                       const aString& name,
                       Real           time);
    //
    // Returns a PArray to fabs containing the derived data.  
    //
    // The user is responsible for deleting this pointer when done with it.
    //
    MultiFab* derive (const aString& name,
                      Real           time);
    //
    // State data object.
    //
    StateData& get_state_data (int state_indx) { return state[state_indx]; }
    //
    // State data at old time.
    //
    MultiFab& get_old_data (int state_indx)
    {
        return state[state_indx].oldData();
    }
    //
    // State data at new time.
    //
    MultiFab& get_new_data(int state_indx)
    {
        return state[state_indx].newData();
    }

    const DescriptorList & get_desc_lst()
    {
        return desc_lst;
    }
    //
    // Boundary condition access function.
    //
    int* getBCArray (int State_Type,
                     int gridno,
                     int strt_comp,
                     int num_comp);

    MultiFab& get_data(int state_indx, Real time);

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
    // FillPatch for derived quantities.
    //
    void FillDerive (FArrayBox&     dest,
                     const Box&     subbox,
                     const aString& name,
                     Real           time);
    //
    // Common code used by all constructors.
    //
    void finishConstructor();    
    //
    // Pointwise state data.
    //
    void probe (ostream&       os,
                const IntVect& iv,
                int            radius,
                Real           time,
                int            state_indx,
                int            src_comp,
                int            num_comp);
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
    FillPatchIterator(AmrLevel &amrlevel, MultiFab &leveldata);

    FillPatchIterator(AmrLevel &amrlevel,
                      MultiFab &leveldata,
                      int boxGrow,
                      int dest_comp, Real time,
                      int state_indx, int src_comp, int ncomp,
                      Interpolater* mapper = 0);

    void Initialize (int           boxGrow,
                     int           dest_comp,
                     Real          time,
                     int           state_indx,
                     int           src_comp,
                     int           ncomp,
                     Interpolater* mapper = 0);

    ~FillPatchIterator();

    const FArrayBox& operator()() const { return currentFillPatchedFab; }

    FArrayBox& operator()()  { return currentFillPatchedFab; }

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
    AmrLevel& amrLevel;
    MultiFab& levelData;
    MultiFabCopyDescriptor multiFabCopyDesc;
    Array< Array<MultiFabId> > stateDataMFId;   // [level][oldnew]
    Array<IntVect> cumulativeRefRatios;             // [level]
    Array<Interpolater *> map;                  // [level]
    Array< Array< Array< Array<FillBoxId> > > > fillBoxId;
                                // [grid][level][fillablesubbox][oldnew]
    Array< Array< Array<Box> > > savedFineBox;  // [grid][level][fillablesubbox]

    BoxArray localMFBoxes;
    int  growSize;
    int  stateIndex;
    int  srcComp;
    int  destComp;
    int  nComp;
    Real interpTime;
    FArrayBox currentFillPatchedFab;
    bool bIsInitialized;
};

#endif
