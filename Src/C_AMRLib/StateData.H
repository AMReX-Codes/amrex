
#ifndef _StateData_H_
#define _StateData_H_ 

#include <Box.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <BCRec.H>
#include <Array.H>
#include <VisMF.H>
#include <DistributionMapping.H>
#include <PhysBCFunct.H>
#include <Geometry.H>

class RealBox;
class StateDescriptor;

//
// Current and previous level-time data.
//
// StateData holds state data on a level for the current and previous time step.
//

class StateData
{
public:
    //
    // The default constructor.
    //
    StateData ();
    //
    // Constructor that properly initializes data members.
    //
    StateData (const Box&             p_domain,
               const BoxArray&        grds,
               const StateDescriptor* d,
               Real                   cur_time,
               Real                   dt);
    //
    // The destructor.
    //
    ~StateData ();
    //
    // Takes an empty StateData object dest and defines it with the same properties
    // as the input StateData src. Does a deep copy of the data from src into dest.
    //
    static void Initialize (StateData& dest, const StateData& src);
    //
    // Initializes data members if you used default constructor.
    //
    void define (const Box&             p_domain,
                 const BoxArray&        grds,
                 const StateDescriptor& d, 
                 Real                   cur_time,
                 Real                   dt);

    void define (const Box&             p_domain,
                 const BoxArray&        grds,
		 const DistributionMapping& dm,
                 const StateDescriptor& d, 
                 Real                   cur_time,
                 Real                   dt);
    //
    // Copies old data from another StateData object and sets the same time level.
    // If old data is uninitialized, allocates it with same properties as the input data.
    //
    void copyOld (const StateData& state);
    //
    // Copies new data from another StateData object and sets the same time level.
    // If new data is uninitialized, allocates it with the same properties as the input data.
    //
    void copyNew (const StateData& state);
    //
    // Allocates space for old timestep data.
    //
    void allocOldData ();
    //
    // Deletes the space used by the old timestep data.
    //
    void removeOldData () { delete old_data; old_data = 0; }
    //
    // Reverts back to initial state.
    //
    void reset ();
    //
    // Old data becomes new data and new time is incremented by dt.
    //
    void swapTimeLevels (Real dt);
    //
    // Swaps old data with a new MultiFab.
    //
    void replaceOldData ( MultiFab* mf );
    //
    // Swaps new data with a new MultiFab.
    //
    void replaceNewData ( MultiFab* mf );

    //
    // Sets time of old and new data.
    //
    void setTimeLevel (Real t_new,
                       Real dt_old,
                       Real dt_new);
    //
    // Sets time of old data.
    //
    void setOldTimeLevel (Real t_old);
    //
    // Sets time of new data.
    //
    void setNewTimeLevel (Real t_new);

    void RegisterData (MultiFabCopyDescriptor& multiFabCopyDesc,
                       Array<MultiFabId>&      mfid);

    void InterpAddBox (MultiFabCopyDescriptor& multiFabCopyDesc,
		       Array<MultiFabId>&      mfid,
		       BoxList*                returnedUnfillableBoxes,
		       Array<FillBoxId>&       returnedFillBoxIds,
		       const Box&              subbox,
		       Real                    time,
		       int                     src_comp,
		       int                     dest_comp,
		       int                     num_comp,
		       bool                    extrap = false);

    void InterpFillFab (MultiFabCopyDescriptor&  fabCopyDesc,
			const Array<MultiFabId>& mfid,
			const Array<FillBoxId>&  fillBoxIds,
			FArrayBox&               dest,
			Real                     time,
			int                      src_comp,
			int                      dest_comp,
			int                      num_comp,
			bool                     extrap = false);

    //
    // Set physical bndry values
    //
    void FillBoundary (FArrayBox&     dest,
                       Real           time,
                       const Real*    dx,
                       const RealBox& prob_domain,
                       int            dest_comp,
                       int            src_comp,
                       int            num_comp = 1);
    //
    // Write the state data to a checkpoint file.
    //
    void checkPoint (const std::string& name,
                     const std::string& fullpathname,
                     std::ostream&      os,
                     VisMF::How         how,
                     bool               dump_old = true);
    //
    // Initializes state data from a checkpoint file.
    //
    void restart (std::istream&          is,
                  const StateDescriptor& d,
                  const std::string&     restart_file,
		  bool                   bReadSpecial = false);
    //
    // or from another similar state
    //
    void restart (const StateDescriptor& d,
		  const StateData& rhs);
    //
    // Returns the StateDescriptor.
    //
    const StateDescriptor* descriptor () const { return desc; }
    //
    // Returns the valid domain.
    //
    const Box& getDomain () const { return domain; }
    //
    // Returns the BoxArray.
    //
    const BoxArray& boxArray () const { return grids; }
    //
    // Returns the current time.
    //
    Real curTime () const { return 0.5*(new_time.start + new_time.stop); }
    //
    // Returns the previous time.
    //
    Real prevTime () const { return 0.5*(old_time.start + old_time.stop); }
    //
    // Returns the new data.
    //
    MultiFab& newData () { BL_ASSERT(new_data != 0); return *new_data; }
    //
    // Returns the new data.
    //
    const MultiFab& newData () const { BL_ASSERT(new_data != 0); return *new_data; }
    //
    // Returns the old data.
    //
    MultiFab& oldData () { BL_ASSERT(old_data != 0); return *old_data; }
    //
    // Returns the old data.
    //
    const MultiFab& oldData () const { BL_ASSERT(old_data != 0); return *old_data; }
    //
    // Returns the FAB of new data at grid index `i'.
    //
    FArrayBox& newGrid (int i) { BL_ASSERT(new_data != 0); return (*new_data)[i]; }
    //
    // Returns the FAB of old data at grid index `i'.
    //
    FArrayBox& oldGrid (int i) { BL_ASSERT(old_data != 0); return (*old_data)[i]; }
    //
    // Returns boundary conditions of specified component on the specified grid.
    //
    BCRec getBC (int comp, int i) const;
    //
    // Prints out the time interval.
    //
    void printTimeInterval (std::ostream& os) const;
    //
    // True if there is any old data available.
    //
    bool hasOldData () const { return old_data != 0; }
    //
    // True if there is any new data available.
    //
    bool hasNewData () const { return new_data != 0; }

    void getData (PArray<MultiFab>& data,
		  std::vector<Real>& datatime,
		  Real time) const;

    void AddProcsToComp(const StateDescriptor &sdPtr,
                        int ioProcNumSCS, int ioProcNumAll,
                        int scsMyId, MPI_Comm scsComm);

    void Check() const;


private:

    struct TimeInterval
    {
        Real start;
        Real stop;
    };
    //
    // Pointer to data descriptor.
    //
    const StateDescriptor* desc;
    //
    // Problem domain.
    //
    Box domain;
    //
    // Grids defined at this level.
    //
    BoxArray grids;
    //
    // Time variable assoc with new data.
    //
    TimeInterval new_time;
    //
    // Time variable assoc with old data.
    //
    TimeInterval old_time;
    //
    // Pointer to new-time data.
    //
    MultiFab* new_data;
    //
    // Pointer to previous time data.
    //
    MultiFab* old_data;
};

class StateDataPhysBCFunct
    : public PhysBCFunct
{
public:
    StateDataPhysBCFunct (StateData& sd, int sc, const Geometry& geom_);
    virtual void doit (MultiFab& mf, int dcomp, int ncomp, Real time);
private:
    StateData* statedata;
    int src_comp;
    const Geometry& geom;
    
};

#endif /*_StateData_H_*/
