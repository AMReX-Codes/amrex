//BL_COPYRIGHT_NOTICE

#ifndef _StateData_H_
#define _StateData_H_ 

//
// $Id: StateData.H,v 1.9 1998-03-30 21:57:34 lijewski Exp $
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
// Current and previous level-time data.
//
// StateData holds state data on a level for the current and previous time step
//
class StateData
{

public:

#ifdef __GNUG__
     bool operator== (const StateData& rhs) const { return true; }
    StateData(const StateData& rhs);
    StateData& operator= (const StateData& rhs);

#endif
    //
    // Bogus constructor.
    //
    StateData ();
    //
    // Constructor.
    //
    StateData (const Box&             p_domain,
               const BoxArray&        grds,
               const StateDescriptor* d,
               Real                   cur_time,
               Real                   dt);
    //
    // Destructor.
    //
    ~StateData ();
    //
    // Initialize data members.
    //
    void define (const Box&             p_domain,
                 const BoxArray&        grds,
                 const StateDescriptor& d, 
                 Real                   cur_time,
                 Real                   dt);

    void buildBC ();

    void allocOldData ();

    void removeOldData ();

    void reset ();

    void swapTimeLevels (Real dt);

    void setTimeLevel (Real t_new,
                       Real dt_old,
                       Real dt_new);

    void RegisterData (MultiFabCopyDescriptor& multiFabCopyDesc,
                       Array<MultiFabId>&      mfid);

    void linInterpAddBox (MultiFabCopyDescriptor& multiFabCopyDesc,
                          Array<MultiFabId>&      mfid,
                          BoxList&                returnedUnfillableBoxes,
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
    // Set physical bndry values for all state variables in range.
    //
    void FillBoundary (FArrayBox&     dest,
                       Real           time,
                       const Real*    dx,
                       const RealBox& prob_domain,
                       int            dest_comp,
                       int            src_comp,
                       int            num_comp = 1);
    //
    // Set physical bndry values for all state variables in range.
    //
    void FillBoundary (const Real*    dx,
                       const RealBox& prob_domain,
                       int            src_comp,
                       int            num_comp = 1,
                       int            do_new = 1);

    void checkPoint (const aString& name,
                     const aString& fullpathname,
                     ostream&       os,
                     VisMF::How     how,
                     bool           dump_old = true);

    void restart (istream&               is,
                  const StateDescriptor& d,
                  const aString&         restart_file);

    const StateDescriptor* descriptor () const { return desc; }

    const Box& getDomain () const { return domain; }

    const BoxArray& boxArray () const { return grids; }

    Real curTime () const { return 0.5*(new_time.start + new_time.stop); }

    Real prevTime () const { return 0.5*(old_time.start + old_time.stop); }

    MultiFab& newData () { return *new_data; }

    MultiFab& oldData () { return *old_data; }

    FArrayBox& newGrid (int i) { return (*new_data)[i]; }

    FArrayBox& oldGrid (int i) { return (*old_data)[i]; }

    Array<BCRec>& getBCs (int comp) { return bc[comp]; }

    const BCRec& getBC (int comp, int i) const { return bc[comp][i]; }

    void printTimeInterval (ostream& os) const;

    int hasOldData () const { return old_data != 0; }

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
