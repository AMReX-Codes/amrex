//BL_COPYRIGHT_NOTICE

//
// $Id: StationData.H,v 1.2 1998-11-25 23:58:22 lijewski Exp $
//

#ifndef _StationData_H_
#define _StationData_H_

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <fstream>
using std::ios;
else
#include <iostream.h>
#include <fstream.h>
#endif

#include "Amr.H"
#include "Array.H"
#include "REAL.H"
#include "aString.H"

struct StationRec
{
    StationRec ()
    {
        id    = -1;
        grd   = 0;
        level = -1;
    }

    Real    pos[BL_SPACEDIM]; // Physical coordinates
    int     id;               // Identifier
    int     level;            // Finest level containing station
    IntVect ix;               // Index position
    Grid*   grd;              // Pointer to grid containing point
};

class StationData
{
public:

    StationData ();

    ~StationData ();
    //
    // Open appropriately-named station file.
    //
    //  "rootname_CPU_NNNN_timestep"
    //
    void openFile (int timestep);
    //
    // Init from ParmParse.
    //
    // ParmParse variables:
    //
    //   StationData.vars     -- Array of name of StateData components
    //   StationData.coord    -- BL_SPACEDIM array of Reals
    //   StationData.coord    -- the next one
    //   StationData.coord    -- ditto ...
    //   StationData.rootname -- root name of output files.
    //
    void init ();
    //
    // Collect data for all station points located
    // at the given level and write to station file.
    //
    void report (Real time) const;
    //
    // Locate finest level grid for each station point.
    //
    void findGrid (const GridList* gl,
                   int             finest_level,
                   Real            dx_lev[MAX_LEV][BL_SPACEDIM],
                   Real            prob_lo[BL_SPACEDIM]);
private:
    //
    // Writes out map of StationIDs to Real coordinates to "Station.List".
    //
    void listStations () const;

    Array<StationRec> m_stn;   // Array of stations.
    Array<aString>    m_vars;  // Names of StateData components to output.
    ofstream          m_ofile; // Output stream.
    aString           m_name;  // Root name of station file.
};

#endif /*_StationData_H_*/

