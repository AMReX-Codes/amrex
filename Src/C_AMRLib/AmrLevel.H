//BL_COPYRIGHT_NOTICE

#ifndef _AmrLevel_H_
#define _AmrLevel_H_ 

//
// $Id: AmrLevel.H,v 1.5 1997-12-02 00:08:55 lijewski Exp $
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

// ----------------------------------------------------------------
// -------------  AmrLevel Class  ---------------------------------
// ----------------------------------------------------------------

/// Virtual base class for managing individual levels
class AmrLevel
{

  friend class FillPatchIterator;


/* AmrLevel functions both as a container for state data on a level
  and also manages the advancement of data in time.
*/

protected:

//  AMR level (0 is coarsest)
    int level;

//  geom at this level
    Geometry geom;

//  cell-centered boxes giving locations of grids
    BoxArray grids;

//  physical locations of each grid
    Array<RealBox> grid_loc;

//  pointer to parent AMR structure
    Amr * parent;

//  refinement ratio to coarser level
    IntVect crse_ratio;

//  refinement ratio to finer level
    IntVect fine_ratio;

//  list of derived quantities
    static DeriveList       derive_lst;

//  list of state variables
    static DescriptorList   desc_lst;

//  array of state data
    Array<StateData> state;

protected:

//  default constructor
    AmrLevel();

//  full constructor
    AmrLevel(Amr &papa, int lev, const Geometry &level_geom,
	     const BoxArray& bl, Real time);

public:
//  destructor
    virtual ~AmrLevel();

#ifdef BL_PARALLEL_IO
//  write current state to checkpoint file
    virtual void checkPoint (const aString& dir,
                             ostream&       os,
                             VisMF::How     how = VisMF::OneFilePerCPU);
#else
//  write current state to checkpoint file
    virtual void checkPoint (ostream& os);
#endif

//  restart from a checkpoint file
    virtual void restart (Amr&     papa,
                          istream& is);

#ifdef BL_PARALLEL_IO
    //
    // A string written as the first item in writePlotFile() at level zero.
    // This MUST be defined by each derived class.
    // It is so we can distinguish between different types of plot files.
    //
    virtual aString thePlotFileType () const = 0;

//  write plot file stuff to specified directory
    virtual void writePlotFile (const aString& dir,
                                ostream&       os,
                                VisMF::How     how = VisMF::OneFilePerCPU) = 0;
#else
//  write plot file
    virtual void writePlotFile (ostream &os) = 0;
#endif

//  is <aString> a state variable?
    int isStateVariable( const aString &name,
			 int& state_indx, int& n_comp);

//  compute the initial time step
    virtual void computeInitialDt(int finest_level, int sub_cycle,
				  Array<int> &n_cycle,
				  const Array<IntVect> &ref_ratio,
				  Array<Real> &dt_level) = 0;

//  compute the next time step
    virtual void computeNewDt(int finest_level, int sub_cycle,
			      Array<int> &n_cycle,
			      const Array<IntVect> &ref_ratio,
			      Array<Real> &dt_min,
			      Array<Real> &dt_level,
			      Real stop_time) = 0;


//  do an integration step on this level - returns maximum safe time step
    virtual  Real advance(Real time, Real dt, int iteration, int ncycle) = 0;

//  contains operations to be done after a timestep (e.g. refluxing)  -- must be    defined in derived class
    virtual  void post_timestep() = 0;

//  contains operations to be done only after a full coarse timestep 
    virtual void postCoarseTimeStep(Real time) {}

//  operations to be done after restart -- must be defined in derived class
    virtual  void post_restart() = 0;

 // operations to be done after regridding (like avgDown) -- must be defined        in derived class
    virtual  void post_regrid(int lbase, int new_finest) = 0;

//  operations to be done after initialization -- must be defined in derived        class
    virtual  void post_init() = 0;

//  is it ok to continue the calculation?
    virtual  int okToContinue() = 0;

//  init grid data at problem start-up
    virtual void initData() = 0;

//  set the time levels of state data
    virtual void setTimeLevel(Real time, Real dt_old, Real dt_new);

//  alloc space for old time data
    virtual void allocOldData();

//  delete old-time data
    virtual void removeOldData();

//  init data on this level from another AmrLevel (during regrid)
    virtual void init(AmrLevel &old) = 0;

//  init data on this level after regridding if old AmrLevel did not previously     exist
    virtual void init() = 0;

//  reset data to initial time by swapping new and old time data
    void reset();

///
/// simple access functions
//  returns this AmrLevel
    int Level() {return level; }

//  list of grids at this level
    const BoxArray& boxArray() const {return grids;}

//  number of grids at this level
    int numGrids() { return grids.length(); }

//  physical locations of grids at this level
    const Array<RealBox>& gridLocations() { return grid_loc; }

//  returns the indices defining physical domain
    const Box& Domain() { return geom.Domain(); }

//  timestep n at this level
    int nStep() { return parent->levelSteps(level); }

//  returns a geometry object
    const Geometry& Geom() {return geom;}

//  returns number of cells on level
    long countCells();

//  error estimation for regridding
    virtual void errorEst(TagBoxArray &tb, int clearval, int tagval,
                          Real time) = 0;
///
 //  general fill function
    void FillPatch(FArrayBox &dest, int dest_comp, Real time,
		   int state_indx, int src_comp, int ncomp,
		   Interpolater *mapper = 0);
///
//  interpolate up from coarse level
    void FillCoarsePatch(FArrayBox &dest,
		         int dest_comp, Real time,
		         int state_indx, int src_comp, int ncomp,
		         Interpolater *mapper = 0);

///
//  function to set physical boundary conditions
    void setPhysBoundaryValues(int state_indx, int comp, int ncomp,
			       int do_new = 1)
    { state[state_indx].FillBoundary(geom.CellSize(),geom.ProbDomain(),
				 comp,ncomp,do_new); }

    void setPhysBoundaryValues(int state_indx, int comp, int ncomp,
			       Real time);

/// functions to derive specific quantities

// returns a pointer to a fab containing the derived data.  
// Note: if you use this, you are responsible for proper memory management
    FArrayBox* derive(const Box& bx, const aString &name, Real time);

//  returns a PArray to fabs containing the derived data.  
//  Note: if you use this, you are responsible for proper memory management
    PArray<FArrayBox>* derive(const aString &name, Real time);

//  state data object
    StateData& get_state_data(int state_indx)
    { return state[state_indx]; }

//  state data at old time
    MultiFab& get_old_data(int state_indx)
    { return state[state_indx].oldData(); }

//  state data at new time
    MultiFab& get_new_data(int state_indx)
    { return state[state_indx].newData(); }

    const DescriptorList & get_desc_lst()
    { return desc_lst; }

//  boundary condition access function
    int* getBCArray(int State_Type, int gridno, int strt_comp, int num_comp);

    MultiFab& get_data(int state_indx, Real time);

protected:

//  FillPatch for derived quantities
    void FillDerive(FArrayBox &dest, const Box& subbox,
		   const aString &name, Real time);

//  common code used by all constructors
    void finishConstructor();    

//  pointwise state data
    void probe(ostream &os, IntVect iv, int radius, Real time,
	       int state_indx, int src_comp, int num_comp);

private:
    AmrLevel(const AmrLevel&);
    AmrLevel& operator = (const AmrLevel&);
};




// ----------------------------------------------------------------------------
class FillPatchIterator : public MultiFabIterator {

  public:

    enum FillType { FillLocally, FillRemotely, FillFromBoundary, Unfillable };

    FillPatchIterator(AmrLevel &amrlevel,
                      MultiFab &leveldata,
                      const int boxGrow,
                      int dest_comp, Real time,
                      int state_indx, int src_comp, int ncomp,
                      Interpolater *mapper = 0);
    ~FillPatchIterator();
    const FArrayBox &operator()() const { return currentFillPatchedFab; }
    FArrayBox &operator()()  { return currentFillPatchedFab; }
    //inline FabArrayIterator &operator++();
    bool isValid(bool bDoSync = true);

    const Box &UngrownBox() const;

  private:
    // no, no, no
    FillPatchIterator();
    FillPatchIterator(const FillPatchIterator &);
    FillPatchIterator &operator=(const FillPatchIterator &);

  private:
    struct FabCopyMap {
      int fabIndex;
      Box commBox;
    };

    AmrLevel &amrLevel;
    MultiFabCopyDescriptor multiFabCopyDesc;
    Array< Array<MultiFabId> > stateDataMFId;   // [level][oldnew]
    //Array<int> cumulativeRefRatios;             // [level]
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

};

// ----------------------------------------------------------------------------


#endif
