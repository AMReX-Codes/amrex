
#ifndef _StateData_H_
#define _StateData_H_ 

#include <Box.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <BCRec.H>
#include <Array.H>
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
    // Initializes data members if you used default constructor.
    //
    void define (const Box&             p_domain,
                 const BoxArray&        grds,
                 const StateDescriptor& d, 
                 Real                   cur_time,
                 Real                   dt);

    //
    // Allocates space for old timestep data.
    //
    void allocOldData ();
    //
    // Deletes the space used by the old timestep data.
    //
    void removeOldData ();
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

    void linInterpAddBox (MultiFabCopyDescriptor& multiFabCopyDesc,
                          Array<MultiFabId>&      mfid,
                          BoxList*                returnedUnfillableBoxes,
                          Array<FillBoxId>&       returnedFillBoxIds,
                          const Box&              subbox,
                          Real                    time,
                          int                     src_comp,
                          int                     dest_comp,
                          int                     num_comp,
                          bool                    extrap = false);

    void linInterpFillFab (MultiFabCopyDescriptor&  fabCopyDesc,
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
    const StateDescriptor* descriptor () const;
    //
    // Returns the valid domain.
    //
    const Box& getDomain () const;
    //
    // Returns the BoxArray.
    //
    const BoxArray& boxArray () const;
    //
    // Returns the current time.
    //
    Real curTime () const;
    //
    // Returns the previous time.
    //
    Real prevTime () const;
    //
    // Returns the new data.
    //
    MultiFab& newData ();
    //
    // Returns the new data.
    //
    const MultiFab& newData () const;
    //
    // Returns the old data.
    //
    MultiFab& oldData ();
    //
    // Returns the old data.
    //
    const MultiFab& oldData () const;
    //
    // Returns the FAB of new data at grid index `i'.
    //
    FArrayBox& newGrid (int i);
    //
    // Returns the FAB of old data at grid index `i'.
    //
    FArrayBox& oldGrid (int i);
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
    bool hasOldData () const;
    //
    // True if there is any new data available.
    //
    bool hasNewData () const;

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

inline
void
StateData::removeOldData ()
{
    delete old_data;
    old_data = 0;
}

inline
const StateDescriptor*
StateData::descriptor () const
{
    return desc;
}

inline
const Box&
StateData::getDomain () const
{
    return domain;
}

inline
const BoxArray&
StateData::boxArray () const
{
    return grids;
}

inline
Real
StateData::curTime () const
{
    return 0.5*(new_time.start + new_time.stop);
}

inline
Real
StateData::prevTime () const
{
    return 0.5*(old_time.start + old_time.stop);
}

inline
MultiFab&
StateData::newData ()
{
    BL_ASSERT(new_data != 0);
    return *new_data;
}

inline
const MultiFab&
StateData::newData () const
{
    BL_ASSERT(new_data != 0);
    return *new_data;
}

inline
MultiFab&
StateData::oldData ()
{
    BL_ASSERT(old_data != 0);
    return *old_data;
}

inline
const MultiFab&
StateData::oldData () const
{
    BL_ASSERT(old_data != 0);
    return *old_data;
}

inline
FArrayBox&
StateData::newGrid (int i)
{
    BL_ASSERT(new_data != 0);
    return (*new_data)[i];
}

inline
FArrayBox&
StateData::oldGrid (int i)
{
    BL_ASSERT(old_data != 0);
    return (*old_data)[i];
}

inline
bool
StateData::hasOldData () const
{
    return old_data != 0;
}

inline
bool
StateData::hasNewData () const
{
    return new_data != 0;
}

#endif /*_StateData_H_*/
