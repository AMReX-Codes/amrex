//BL_COPYRIGHT_NOTICE

//
// $Id: StationData.H,v 1.4 1998-11-29 00:21:41 lijewski Exp $
//

#ifndef _StationData_H_
#define _StationData_H_

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <fstream>
using std::ios;
#else
#include <iostream.h>
#include <fstream.h>
#endif

#include <Amr.H>
#include <Array.H>
#include <REAL.H>
#include <StateData.H>
#include <Utility.H>
#include <aString.H>

//
// Forward declaration.
//
class AmrLevel;

struct StationRec
{
    //
    // A constructor just to keep 3rd quiet.
    //
    StationRec()
    {
        D_TERM(pos[0],=pos[1],=pos[2]) = -1;
        id = level = grd = -1;
        own = false;
    }

    Real    pos[BL_SPACEDIM]; // Physical coordinates
    IntVect ix;               // Index position
    int     id;               // Identifier
    int     level;            // Finest level containing station
    int     grd;              // Index of BoxArray containing FAB
    bool    own;              // Does this CPU own the data at grd?
};

class StationData
{
public:

    StationData ();

    ~StationData ();
    //
    // Init from ParmParse.
    //
    // ParmParse variables:
    //
    //   StationData.vars     -- Array of names of StateData components
    //   StationData.coord    -- BL_SPACEDIM array of Reals
    //   StationData.coord    -- the next one
    //   StationData.coord    -- ditto ...
    //   StationData.rootname -- root name of output files.
    //
    // Default rootname is: "stn"
    //
    void init ();
    //
    // Collect data for all station points at level and write to station file.
    //
    void report (Real            time,
                 int             level,
                 const AmrLevel& amrlevel);
    //
    // Locate finest level grid for each station point.
    //
    void findGrid (const PArray<AmrLevel>& levels,
                   const Array<Geometry>&  geoms);
private:
    //
    // Open appropriately-named station file.
    //
    //  "rootname_CPU_NNNN"
    //
    // Default rootname is: "stn"
    //
    void openFile ();
    //
    // Writes out map of StationIDs to Real coordinates to "Station.List".
    //
    void listStations () const;

    Array<StationRec> m_stn;   // Array of stations.
    Array<aString>    m_vars;  // Names of StateData components to output.
    Array<int>        m_typ;   // The state_index corresponding to m_vars.
    Array<int>        m_ncomp; // The component of the state_index for m_typ.
    ofstream          m_ofile; // Output stream.
    aString           m_name;  // Root name of station file.
};

#endif /*_StationData_H_*/

