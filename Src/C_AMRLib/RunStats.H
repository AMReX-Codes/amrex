#ifndef _RUNSTATS_H_
#define _RUNSTATS_H_ 

//
// $Id: RunStats.H,v 1.1 1997-11-18 19:30:28 lijewski Exp $
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
    RunStats(const char *name, int level = 0);
    ~RunStats() {};
    
    friend ostream& operator << (ostream &os, const RunStats &r);

    void start();		// start timeing given variable
    void pause();		// pause timing
    void resume();		// resume after a pause
    void end();			// stop timeing given variable and record

    int isOn() const {
	return gentry->is_on && entry->is_on;
    }

      // record cells advanced
    static void addCells(int level, long count);

      // return number of cells advanced at given level 
    static long getCells(int level);

      // turn on recording of this stat
    static void turnOn(const char * s, int level = -1);
    static void turnOff(const char * s, int level = -1);

    static void init();
      // report stats in formatted form to output stream
    static void report(ostream &os);
      // unformatted write to output stream
    static void dumpStats(ofstream &os);
      // unformatted read from input stream
    static void readStats(ifstream &is);

private:
    static RunStatsData* find(aString _name, int _level);
    RunStatsData *entry;
    RunStatsData *gentry;
    aString name;
    int level;
    double time;
    double wtime;
    void resize(int level);
    static double total_run_time;
    static double total_run_wtime;
    static Array<long> cells;
    static List<RunStatsData> ld;
};

#endif

