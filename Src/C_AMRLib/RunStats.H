#ifndef _RUNSTATS_H_
#define _RUNSTATS_H_ 

//
// $Id: RunStats.H,v 1.2 1997-11-19 16:25:55 lijewski Exp $
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

class RunStatsData
{
    friend class RunStats;
    friend ostream & operator << (ostream &, const RunStatsData&);
    friend istream & operator >> (istream &, RunStatsData&);
private:
    aString name;
    int level;
    int is_on;
    double run_time;
    double run_wtime;
    double max_time;
    double max_wtime;
public:
    RunStatsData()
	: run_time(0.0), run_wtime(0.0),
	  max_time(0.0), max_wtime(0.0)
    {}
    int operator == (const RunStatsData &s1) const {
	return (name == s1.name) && (level == s1.level);
    }
    int operator != (const RunStatsData &s1) const {
	return (name != s1.name) || (level != s1.level);
    }
    RunStatsData(aString _name, int _level) 
	: run_time(0.0), run_wtime(0.0),
	  max_time(0.0), max_wtime(0.0),
	  name(_name), level(_level)
    {}
};

class RunStats
{
public:
    RunStats (const char *name, int level = 0);
    ~RunStats () {}
    
    friend ostream& operator << (ostream &os, const RunStats &r);

    void start ();		// start timeing given variable
    void pause ();		// pause timing
    void resume ();		// resume after a pause
    void end ();		// stop timeing given variable and record

    int isOn () const;

      // record cells advanced
    static void addCells (int level, long count);

      // return number of cells advanced at given level 
    static long getCells (int level);

      // turn on recording of this stat
    static void turnOn (const char * s, int level = -1);
    static void turnOff (const char * s, int level = -1);

    static void init ();
      // report stats in formatted form to output stream
    static void report (ostream &os);
      // unformatted write to output stream
    static void dumpStats (ofstream &os);
      // unformatted read from input stream
    static void readStats (ifstream &is);

private:
    static RunStatsData* find (const char* _name, int _level);
    void resize (int level);

    RunStatsData *entry;
    RunStatsData *gentry;
    aString name;
    int level;
    double time;
    double wtime;

    static double total_run_time;
    static double total_run_wtime;
    static Array<long> cells;
    static List<RunStatsData> ld;
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
    return cells[lev];
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

