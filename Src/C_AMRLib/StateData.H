//BL_COPYRIGHT_NOTICE

#ifndef _StateData_H_
#define _StateData_H_ 

//
// $Id: StateData.H,v 1.2 1997-11-24 18:52:32 lijewski Exp $
//

#include <Box.H>
#include <RealBox.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <BCRec.H>
#include <Array.H>
#include <StateDescriptor.H>
#include <VisMF.H>

/// Current and previous level-time data
class StateData
{

/** StateData holds state data on a level
for the current and previous time step
*/

    struct TimeInterval {
	REAL start, stop;
    };

protected:

//  pointer to data descriptor
    const StateDescriptor  *desc;

//  problem domain
    BOX        domain;

//  grids defined at this level
    BoxArray           grids;

 // time variable assoc with new data
    TimeInterval     new_time;

//  time variable assoc with old data
    TimeInterval     old_time;

//  pointer to new-time data
    MultiFab           *new_data;

//  pointer to previous time data
    MultiFab           *old_data;

//  boundary conditions for the grids at this level
    Array< Array<BCRec> > bc;

public:
#ifdef __GNUG__
     bool operator== (const StateData& rhs) const { return true; }
    StateData(const StateData& rhs);
    StateData& operator= (const StateData& rhs);

#endif 
//  bogus constructor 
    StateData();

//  constructor
    StateData(const BOX& p_domain, const BoxArray& grds,
              const StateDescriptor* d,
              REAL cur_time, REAL dt);
//  destructor
    ~StateData();

//  initialize data members
    void define(const BOX& p_domain, const BoxArray& grds,
                const StateDescriptor &d, 
                REAL cur_time, REAL dt);
//
    void restart(istream& is, const StateDescriptor& d);

//
    void buildBC();
//
    void allocOldData();
//
    void removeOldData();
//
    void reset();
//
    void swapTimeLevels(REAL dt);
//
    void setTimeLevel(REAL t_new, REAL dt);

    void RegisterData(MultiFabCopyDescriptor &multiFabCopyDesc,
                      Array<MultiFabId> &mfid);

//  interpolation routine
    void linInterp(FARRAYBOX& dest, const BOX& subbox, REAL time,
                   int src_comp, int dest_comp, int num_comp,
                   bool extrap = false);

    void linInterpAddBox(MultiFabCopyDescriptor &multiFabCopyDesc,
                   Array<MultiFabId> &mfid,
                   BoxList &returnedUnfillableBoxes,
                   Array<FillBoxId> &returnedFillBoxIds,
                   const BOX &subbox,
                   REAL time, int src_comp, int dest_comp, int num_comp,
                   bool extrap = false);

    void linInterpFillFab(MultiFabCopyDescriptor &fabCopyDesc,
                   const Array<MultiFabId> &mfid,
                   const Array<FillBoxId> &fillBoxIds,
                   FARRAYBOX &dest,
                   REAL time,
                   int src_comp, int dest_comp, int num_comp,
                   bool extrap = false);

//  set physical bndry values for all state variables in range
    void FillBoundary(FARRAYBOX& dest, REAL time, const REAL* dx,
                   const RealBox& prob_domain, int dest_comp,
                   int src_comp, int num_comp = 1);

//  set physical bndry values for all state variables in range
    void FillBoundary(const REAL* dx, const RealBox& prob_domain,
                   int src_comp, int num_comp = 1, int do_new=1);

#ifdef BL_PARALLEL_IO
    void checkPoint (const aString& basename,
                     ostream&       os,
                     VisMF::How     how,
                     bool           dump_old = true);
#else
    void checkPoint (ostream& os,
                     bool     dump_old = true);
#endif

/// access functions
///
    const StateDescriptor* descriptor() {return desc; }
///
    const BOX& getDomain() { return domain; }
///
    const BoxArray& boxArray() { return grids; }
///
    REAL curTime() { return 0.5*(new_time.start + new_time.stop); }
///
    REAL prevTime() { return 0.5*(old_time.start + old_time.stop); }
///
    MultiFab& newData() { return *new_data; };
///
    MultiFab& oldData() { return *old_data; };
///
    FARRAYBOX& newGrid(int i) { return (*new_data)[i]; }
///
    FARRAYBOX& oldGrid(int i) { return (*old_data)[i]; }
///
    Array<BCRec>& getBCs(int comp) { return bc[comp]; }
///
    const BCRec& getBC( int comp, int i) { return bc[comp][i]; }
///
    void printTimeInterval(ostream &os);
///
    int hasOldData() { 	return (old_data != 0); }
};

#endif
