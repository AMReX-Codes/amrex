//BL_COPYRIGHT_NOTICE

#ifndef _RUNSTATS_H_
#define _RUNSTATS_H_ 

//
// $Id: RunStats.H,v 1.8 1999-07-01 19:44:49 lijewski Exp $
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
#include <ParallelDescriptor.H>

class RunStatsData
{
    friend class RunStats;
    friend ostream& operator<< (ostream&, const RunStatsData&);
    friend istream& operator>> (istream&, RunStatsData&);

private:

    RunStatsData ()
        :
        run_time(0.0),
        run_wtime(0.0),
        level(0),
        is_on(0)
    {}

    RunStatsData (const aString& _name, int _level) 
        :
        name(_name),
        run_time(0.0),
        run_wtime(0.0),
        level(_level),
        is_on(0)
    {}

    aString name;
    Real    run_time;
    Real    run_wtime;
    int     level;
    int     is_on;

};

//
//@Man:
//@Memo: Collects run-time statistics.
/*@Doc:

  This class facilitates the generation, collection, and display of timing
  statistics for various parts of a program.
*/

class RunStats
{
public:
    //
    //@ManDoc: The constructor.
    //
    RunStats (const aString& name, int level = 0);
    //
    //@ManDoc: The destructor.
    //
    ~RunStats ();
    //
    //@ManDoc: A friendly output operator.
    //
    friend ostream& operator<< (ostream& os, const RunStats& r);
    //
    //@ManDoc: Start timing given variable.
    //
    void start ();
    //
    //@ManDoc: Stop timing given variable and record.
    //
    void end ();
    //
    //@ManDoc: Is the variable active?
    //
    int isOn () const;
    //
    //@ManDoc: Record number of cells advanced at given level.
    //
    static void addCells (int level, long count);
    //
    //@ManDoc: Return number of cells advanced at given level.
    //
    static long getCells (int level);
    //
    //@ManDoc: Record I/O bytes written to disk.
    //
    static void addBytes (long count);
    //
    //@ManDoc: Record numPts() of each FAB in each MultiFab on this CPU.
    //
    static void addNumPts (long count);
    //
    //@ManDoc: Report stats in formatted form to output stream.
    //
    static void report (ostream &os);
    //
    //@ManDoc: Unformatted write to output stream.
    //
    static void dumpStats (ofstream &os);
    //
    //@ManDoc: Unformatted read from input stream.
    //
    static void readStats (ifstream &is, bool restart = false);
    //
    //@ManDoc: Report the names of the statistics as an array of aString.
    //
    static void report_names (Array<aString>& stat_names);

    /*@ManDoc: Report the values of the statistics as arrays of REALS.
               This assumes `stat\_names' is result of previous call to
               report\_names().
    */
    static void report_values (const Array<aString>& stat_names,
                               Array<Real>&          stat_time,
                               Array<Real>&          stat_wtime,
                               Real&                 tot_run_time,
                               Real&                 tot_run_wtime,
                               long&                 tot_cells);
private:
    //
    // Init from ParmParse with "RunStats" prefix.
    //
    static void init ();
    //
    // Do we have an entry by the specified name?
    //
    static RunStatsData* find (const aString& _name, int _level);

    void resize (int level);

    static void ReduceIt (List<RunStatsData>& stats);

    static void Print (ostream&, const RunStatsData&, Real, Real);

    static void CollectNumPts ();
    //
    // Local data.
    //
    RunStatsData* entry;
    RunStatsData* gentry;
    aString       name;
    Real          time;
    Real          wtime;
    int           level;
    //
    // Static data.
    //
    static Real               TotalCPU;    // Total CPU time.
    static Real               TotalWCT;    // Total Wall Clock time.
    static Real               DiskBytes;   // Total bytes written to disk.
    static Array<long>        TheCells;    // Total cells advanced per level.
    static Array<Real>        TheNumPts;   // Total numPts() per CPU.
    static List<RunStatsData> TheStats;    // List of RunStatsData.
    static bool               Initialized; // init()'s been called?
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

#endif /*_RUNSTATS_H_*/
