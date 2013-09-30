
#ifndef _StateData_H_
#define _StateData_H_ 

#include <Box.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <BCRec.H>
#include <Array.H>
#include <PArray.H>
#include <VisMF.H>

class RealBox;
class StateDescriptor;

namespace BoxLib
{
    void readBoxArray (BoxArray&, std::istream&, bool);
}
//
// Current and previous level-time data.
//
// StateData holds state data on a level for the current and previous time step.
// It may also hold additional node data for SDC or RK.
//

const Array<Real> NullArrayReal;

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
    // Constructor that properly initializes data members including node data.
    //
    StateData (const Box&             p_domain,
               const BoxArray&        grds,
               const StateDescriptor* d,
               Real                   cur_time,
               Real                   dt,
	       const Array<Real>&     t_nodes);
    //
    // The destructor.
    //
    ~StateData ();
    //
    // Initializes data members if you used default constructor.
    //
    void define (const Box&             p_domain,
                 const BoxArray&        grds,
                 const StateDescriptor& d, 
                 Real                   cur_time,
                 Real                   dt,
		 const Array<Real>&     t_nodes = NullArrayReal);

    //
    // Allocates space for old timestep data.
    //
    void allocOldData ();
    //
    // Deletes the space used by the old timestep data.
    //
    void removeOldData () { delete old_data; old_data = 0; }
    //
    // Allocates space for mid timestep data.
    //
    void allocMidData ();
    void setMidData (MultiFab* md, Real time);
    //
    // Deletes the space used by the mid timestep data.
    //
    void removeMidData () { mid_data.clear(); }
    //
    // Reverts back to initial state.
    //
    void reset ();
    //
    // Old data becomes new data and new time is incremented by dt.
    //
    void swapTimeLevels (Real dt);
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
    // Returns i-th mid time.
    //
    Real midTime (int i) const { return 0.5*(mid_time[i].start + mid_time[i].stop); }  
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
    // Returns the i-th mid data.
    //
    MultiFab& midData (int i) { BL_ASSERT(mid_data[i].defined(i)); return mid_data[i]; }
    //
    // Returns the i-th mid data.
    //
    const MultiFab& midData (int i) const { BL_ASSERT(mid_data[i].defined(i)); return mid_data[i]; }
    //
    // Return the size of the mid data
    //
    int sizeMidData () const { return mid_data.size(); }
    //
    // Returns the FAB of new data at grid index `i'.
    //
    FArrayBox& newGrid (int i) { BL_ASSERT(new_data != 0); return (*new_data)[i]; }
    //
    // Returns the FAB of old data at grid index `i'.
    //
    FArrayBox& oldGrid (int i) { BL_ASSERT(old_data != 0); return (*old_data)[i]; }
    //
    // Returns the FAB of 'imid'-th mid data at grid index `i'.
    //
    FArrayBox& midGrid (int imid, int i) { BL_ASSERT(mid_data.defined(imid)); return mid_data[imid][i]; }
    //
    // Returns boundary conditions of specified component on the specified grid.
    //
    const BCRec getBC (int comp, int i) const;
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
    //
    // True if there is mid data available. 
    //
    bool hasMidData (int i=0) const { return mid_data.defined(i); }

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
    // Number of middle nodes
    //
    int n_mid_nodes;
    //
    // middle node times normilized to (0,1)
    //
    Array<Real> t_mid_nodes;
    //
    // "Physical" Time variable assoc with middle nodes
    //
    Array<TimeInterval> mid_time;
    //
    // Pointer to new-time data.
    //
    MultiFab* new_data;
    //
    // Pointer to previous time data.
    //
    MultiFab* old_data;
    //
    // Array of pointers to data at middle times
    //
    PArray<MultiFab> mid_data;
};

#endif /*_StateData_H_*/
