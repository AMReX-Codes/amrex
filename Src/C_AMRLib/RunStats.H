#ifndef _RUNSTATS_H_
#define _RUNSTATS_H_ 

//
// $Id: RunStats.H,v 1.3 1997-11-21 03:50:11 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>
using std::ifstream;
using std::ofstream;
#else
#include <stdlib.h>
#include <fstream.h>
#include <string.h>
#include <math.h>
#endif

#include <aString.H>
#include <List.H>
#include <Array.H>

//
// Forward declaration.
//

class RunStatsData
{
    friend class RunStats;
    friend ostream& operator<< (ostream &, const RunStatsData&);
    friend istream& operator>> (istream &, RunStatsData&);
private:
    aString name;
    double  run_time;
    double  run_wtime;
    int     level;
    int     is_on;

    RunStatsData ()
	:
        run_time(0.0),
        run_wtime(0.0),
        level(0),
        is_on(0)
    {}
    int operator== (const RunStatsData &s1) const
    {
	return (name == s1.name) && (level == s1.level);
    }
    int operator!= (const RunStatsData &s1) const
    {
	return (name != s1.name) || (level != s1.level);
    }
    RunStatsData (const char* _name, int _level) 
	:
        name(_name),
        run_time(0.0),
        run_wtime(0.0),
        level(_level),
        is_on(0)
    {}
};

class RunStats
{
public:
    RunStats (const char *name, int level = 0);

    ~RunStats ();
    
    friend ostream& operator<< (ostream &os, const RunStats &r);
    //
    // Start timeing given variable.
    //
    void start ();
    //
    // Stop timing given variable and record.
    //
    void end ();
    //
    // Is the variable active?
    //
    int isOn () const;
    //
    // Record cells advanced.
    //
    static void addCells (int level, long count);
    //
    // Return number of cells advanced at given level.
    //
    static long getCells (int level);
    //
    // Turn on/off recording of this stat.
    //
    static void turnOn (const char * s, int level = -1);
    static void turnOff (const char * s, int level = -1);
    //
    // Init from ParmParse with "RunStats" prefix.
    //
    static void init ();
    //
    // Report stats in formatted form to output stream.
    //
    static void report (ostream &os);
    //
    // Unformatted write to output stream.
    //
    static void dumpStats (ofstream &os);
    //
    // Unformatted read from input stream.
    //
    static void readStats (ifstream &is);

private:
    static RunStatsData* find (const char* _name, int _level);

    void resize (int level);
    //
    // Local data.
    //
    RunStatsData* entry;
    RunStatsData* gentry;
    aString       name;
    int           level;
    double        time;
    double        wtime;
    //
    // Static data.
    //
    static double             TotalCPU; // Total CPU time
    static double             TotalWCT; // Total Wall Clock time
    static Array<long>        TheCells; // Cells advanced per level
    static List<RunStatsData> TheStats; // List of RunStatsData
};

//
// Inlines
//

inline
int
RunStats::isOn () const
{
    return gentry->is_on && entry->is_on;
}

inline
long
RunStats::getCells (int lev)
{
    return TheCells[lev];
}

inline
void
RunStats::turnOn (const char* s,
                  int         _level)
{
    find(s, _level)->is_on = true;
}

inline
void
RunStats::turnOff (const char* s,
                   int         _level)
{
    if (!(_level == -1))
    {
	find(s, _level)->is_on = false;
    }
}

#endif

