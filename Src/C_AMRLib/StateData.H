//BL_COPYRIGHT_NOTICE

#ifndef _StateData_H_
#define _StateData_H_ 

//
// $Id: StateData.H,v 1.11 1998-05-06 17:50:50 lijewski Exp $
//

#include <Box.H>
#include <RealBox.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <BCRec.H>
#include <Array.H>
#include <StateDescriptor.H>
#include <VisMF.H>

//
//@Man:
//@Memo: Current and previous level-time data.
/*@Doc:

  StateData holds state data on a level for the current and previous time step.
*/

class StateData
{

public:

#ifdef __GNUG__
    bool operator== (const StateData& rhs) const { return true; }
    StateData(const StateData& rhs);
    StateData& operator= (const StateData& rhs);
#endif
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
               Real                   dt)
    {
        define(p_domain, grds, *d, cur_time, dt);
    }
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
    void allocOldData ()
    {
        if (old_data == 0)
        {
            old_data = new MultiFab(grids,desc->nComp(),desc->nExtra(),Fab_allocate);
        }
    }
    //
    //@ManDoc: Deletes the space used by the old timestep data.
    //
    void removeOldData ()
    {
        delete old_data;
        old_data = 0;
    }
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
    void checkPoint (const aString& name,
                     const aString& fullpathname,
                     ostream&       os,
                     VisMF::How     how,
                     bool           dump_old = true);
    //
    //@ManDoc: Initializes state data from a checkpoint file.
    //
    void restart (istream&               is,
                  const StateDescriptor& d,
                  const aString&         restart_file);
    //
    //@ManDoc: Returns the StateDescriptor.
    //
    const StateDescriptor* descriptor () const { return desc; }
    //
    //@ManDoc: Returns the valid domain.
    //
    const Box& getDomain () const { return domain; }
    //
    //@ManDoc: Returns the BoxArray.
    //
    const BoxArray& boxArray () const { return grids; }
    //
    //@ManDoc: Returns the current time.
    //
    Real curTime () const { return 0.5*(new_time.start + new_time.stop); }
    //
    //@ManDoc: Returns the previous time.
    //
    Real prevTime () const { return 0.5*(old_time.start + old_time.stop); }
    //
    //@ManDoc: Returns the new data.
    //
    MultiFab& newData () { return *new_data; }
    //
    //@ManDoc: Returns the old data.
    //
    MultiFab& oldData () { return *old_data; }
    //
    //@ManDoc: Returns the FAB of new data at grid index `i'.
    //
    FArrayBox& newGrid (int i) { return (*new_data)[i]; }
    //
    //@ManDoc: Returns the FAB of old data at grid index `i'.
    //
    FArrayBox& oldGrid (int i) { return (*old_data)[i]; }
    //
    //@ManDoc: Returns the boundary conditions of specified component.
    //
    Array<BCRec>& getBCs (int comp) { return bc[comp]; }
    //
    //@ManDoc: Returns boundary conditions of specified component and direction.
    //
    const BCRec& getBC (int comp, int i) const { return bc[comp][i]; }
    //
    //@ManDoc: Prints out the time interval.
    //
    void printTimeInterval (ostream& os) const;
    //
    //@ManDoc: True if there is any old data available.
    //
    bool hasOldData () const { return old_data != 0; }

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
