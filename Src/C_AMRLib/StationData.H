
//
// $Id: StationData.H,v 1.9 2001-07-23 20:06:06 vince Exp $
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

    Real pos[BL_SPACEDIM]; // Physical coordinates
    int  id;               // Identifier
    int  level;            // Finest level containing station
    int  grd;              // Index of BoxArray containing FAB
    bool own;              // Does this CPU own the data at grd?
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
    //   StationData.vars     -- Names of StateData components to output
    //   StationData.coord    -- BL_SPACEDIM array of Reals
    //   StationData.coord    -- the next one
    //   StationData.coord    -- ditto ...
    //
    // Data files have the form: "Station/stn_CPU_NNNN"
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

    Array<StationRec> m_stn;   // Array of stations.
    Array<aString>    m_vars;  // Names of StateData components to output.
    Array<bool>       m_IsDerived;  // true if the variable is a derived quantity
    Array<int>        m_typ;   // The state_index corresponding to m_vars.
    Array<int>        m_ncomp; // The component of the state_index for m_typ.
    ofstream          m_ofile; // Output stream.
};

#endif /*_StationData_H_*/

