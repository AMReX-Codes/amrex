//BL_COPYRIGHT_NOTICE

//
// $Id: StationData.H,v 1.3 1998-11-28 21:58:43 lijewski Exp $
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
    //   StationData.vars     -- Array of names of StateData components
    //   StationData.coord    -- BL_SPACEDIM array of Reals
    //   StationData.coord    -- the next one
    //   StationData.coord    -- ditto ...
    //   StationData.rootname -- root name of output files.
    //
    void init ();
    //
    // Collect data for all station points at level and write to station file.
    //
    void report (Real time,
                 int  level) const;
    //
    // Locate finest level grid for each station point.
    //
    void findGrid (const PArray<AmrLevel>& levels,
                   const Array<Geometry>&  geoms);

private:
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

