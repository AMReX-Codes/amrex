//BL_COPYRIGHT_NOTICE

//
// $Id: StationData.H,v 1.1 1998-11-25 19:14:16 lijewski Exp $
//

#ifndef _StationData_H_
#define _StationData_H_

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
using std::vector
else
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <vector.h>
#endif

#include "Amr.H"
#include "REAL.H"
#include "aString.H"

struct StationRec
{
    StationRec ()
    {
        level = -1;
        grd   = 0;
        id    = -1;
    }

    int     id;               // Identifier
    Real    pos[BL_SPACEDIM]; // Physical coordinates
    int     level;            // Finest level containing station
    Grid*   grd;              // Pointer to grid containing point
    IntVect ix;               // Index position
};

class StationData
{
public:

    StationData ();

    ~StationData ();
    //
    // Open station file.
    //
    void openFile (int suffix);
    //
    // Alloc list of stations and open output stream.
    //
    void init (int nvar, const char* root_name);
    //
    // Add station at this physical location.
    //
    void addStation (Real* coord);
    //
    // Return current number of stations.
    //
    int numStation () const { return m_stn.size(); }
    //
    // List all current stations by id and position.
    //
    void listStations (ofstream& os);
    //
    // List all current stations by id and position.  Write to default file.
    //
    void listStations ();
    //
    // Write station info to a checkpoint file.
    //
    void dumpStation (ofstream& os, int step) ;
    //
    // Read station info from a restart file.
    //
    void readStation (ifstream& is, int step);
    //
    // Collect data for all station points located
    // at the given level and write to station file.
    //
    void report (int level, Real time, int nvar);
    //
    // Locate finest level grid for each station point.
    //
    void findGrid (const GridList* gl,
                   int             finest_level,
                   Real            dx_lev[MAX_LEV][BL_SPACEDIM],
                   Real            prob_lo[BL_SPACEDIM]);
private:

    vector<StationRec> m_stn;   // Array of stations.
    ofstream           m_ofile; // Output stream.
    aString            m_name;  // Root name of station file.
    int                m_nvar;  // Number of variables at data point.
};

#endif /*_StationData_H_*/

