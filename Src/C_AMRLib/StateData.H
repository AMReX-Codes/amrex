#ifndef _StateData_H_
#define _StateData_H_ 
//
// $Id: StateData.H,v 1.19 2001-08-01 21:50:46 lijewski Exp $
//
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
//@Man:
//@Memo: Current and previous level-time data.
/*@Doc:

  StateData holds state data on a level for the current and previous time step.
*/

class StateData
{
public:
    //
    //@ManDoc: The default constructor.
    //
    StateData ();
    //
    //@ManDoc: Constructor that properly initializes data members.
    //
    StateData (const Box&             p_domain,
               const BoxArray&        grds,
               const StateDescriptor* d,
               Real                   cur_time,
               Real                   dt);
    //
    //@ManDoc: The destructor.
    //
    ~StateData ();
    //
    //@ManDoc: Initializes data members if you used default constructor.
    //
    void define (const Box&             p_domain,
                 const BoxArray&        grds,
                 const StateDescriptor& d, 
                 Real                   cur_time,
                 Real                   dt);
    //
    //@ManDoc: Sets up the boundary conditions.
    //
    void buildBC ();
    //
    //@ManDoc: Allocates space for old timestep data.
    //
    void allocOldData ();
    //
    //@ManDoc: Deletes the space used by the old timestep data.
    //
    void removeOldData ();
    //
    //@ManDoc: Reverts back to initial state.
    //
    void reset ();
    //
    //@ManDoc: Old data becomes new data and new time is incremented by dt.
    //
    void swapTimeLevels (Real dt);
    //
    //@ManDoc: Sets time of old and new data.
    //
    void setTimeLevel (Real t_new,
                       Real dt_old,
                       Real dt_new);
    //
    //@ManDoc: Sets time of old data.
    //
    void setOldTimeLevel (Real t_old);
    //
    //@ManDoc: Sets time of new data.
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
    //@ManDoc: Set physical bndry values for all state variables in range.
    //
    void FillBoundary (FArrayBox&     dest,
                       Real           time,
                       const Real*    dx,
                       const RealBox& prob_domain,
                       int            dest_comp,
                       int            src_comp,
                       int            num_comp = 1);
    //
    //@ManDoc: Set physical bndry values for all state variables in range.
    //
    void FillBoundary (const Real*    dx,
                       const RealBox& prob_domain,
                       int            src_comp,
                       int            num_comp = 1,
                       int            do_new = 1);
    //
    //@ManDoc: Write the state data to a checkpoint file.
    //
    void checkPoint (const std::string& name,
                     const std::string& fullpathname,
                     std::ostream&      os,
                     VisMF::How         how,
                     bool               dump_old = true);
    //
    //@ManDoc: Initializes state data from a checkpoint file.
    //
    void restart (std::istream&          is,
                  const StateDescriptor& d,
                  const std::string&     restart_file,
		  bool                   bReadSpecial = false);
    //
    //@ManDoc: Returns the StateDescriptor.
    //
    const StateDescriptor* descriptor () const;
    //
    //@ManDoc: Returns the valid domain.
    //
    const Box& getDomain () const;
    //
    //@ManDoc: Returns the BoxArray.
    //
    const BoxArray& boxArray () const;
    //
    //@ManDoc: Returns the current time.
    //
    Real curTime () const;
    //
    //@ManDoc: Returns the previous time.
    //
    Real prevTime () const;
    //
    //@ManDoc: Returns the new data.
    //
    MultiFab& newData ();
    //
    //@ManDoc: Returns the new data.
    //
    const MultiFab& newData () const;
    //
    //@ManDoc: Returns the old data.
    //
    MultiFab& oldData ();
    //
    //@ManDoc: Returns the old data.
    //
    const MultiFab& oldData () const;
    //
    //@ManDoc: Returns the FAB of new data at grid index `i'.
    //
    FArrayBox& newGrid (int i);
    //
    //@ManDoc: Returns the FAB of old data at grid index `i'.
    //
    FArrayBox& oldGrid (int i);
    //
    //@ManDoc: Returns an array of boundary condition records for the 
    //         specified component.  Each element in the array corresponds
    //         to a grid in the state.
    //
    Array<BCRec>& getBCs (int comp);
    //
    //@ManDoc: Returns boundary conditions of specified component on the
    //         specified grid.
    //
    const BCRec& getBC (int comp, int i) const;
    //
    //@ManDoc: Prints out the time interval.
    //
    void printTimeInterval (std::ostream& os) const;
    //
    //@ManDoc: True if there is any old data available.
    //
    bool hasOldData () const;
    //
    //@ManDoc: True if there is any new data available.
    //
    bool hasNewData () const;

protected:

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
    //
    // Boundary conditions for the grids at this level.
    //
    Array< Array<BCRec> > bc;
};

#endif /*_StateData_H_*/
