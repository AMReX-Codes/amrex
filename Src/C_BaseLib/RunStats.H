#ifndef _RUNSTATS_H_
#define _RUNSTATS_H_ 
//
// $Id: RunStats.H,v 1.16 2001-07-26 20:08:46 lijewski Exp $
//
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>
#include <list>

#include <Array.H>
#include <ParallelDescriptor.H>

class RunStatsData
{
    friend class RunStats;
    friend std::ostream& operator<< (std::ostream&, const RunStatsData&);
    friend std::istream& operator>> (std::istream&, RunStatsData&);

private:

    RunStatsData ()
        :
        run_time(0.0),
        run_wtime(0.0),
        level(0),
        is_on(0)
    {}

    RunStatsData (const std::string& _name, int _level) 
        :
        name(_name),
        run_time(0.0),
        run_wtime(0.0),
        level(_level),
        is_on(0)
    {}

    std::string name;
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
    RunStats (const std::string& name, int level = 0);
    //
    //@ManDoc: The destructor.
    //
    ~RunStats ();
    //
    //@ManDoc: A friendly output operator.
    //
    friend std::ostream& operator<< (std::ostream& os, const RunStats& r);
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
    static void report (std::ostream &os);
    //
    //@ManDoc: Unformatted write to output stream.
    //
    static void dumpStats (std::ofstream &os);
    //
    //@ManDoc: Unformatted read from input stream.
    //
    static void readStats (std::ifstream &is, bool restart = false);
    //
    //@ManDoc: Report the names of the statistics as an array of std::string.
    //
    static void report_names (Array<std::string>& stat_names);

    /*@ManDoc: Report the values of the statistics as arrays of REALS.
               This assumes `stat\_names' is result of previous call to
               report\_names().
    */
    static void report_values (const Array<std::string>& stat_names,
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
    static RunStatsData* find (const std::string& _name, int _level);

    void resize (int level);

    static void ReduceIt (std::list<RunStatsData>& stats);

    static void Print (std::ostream&, const RunStatsData&, Real, Real);

    static void CollectNumPts ();
    //
    // Local data.
    //
    RunStatsData* entry;
    RunStatsData* gentry;
    std::string       name;
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
    static std::list<RunStatsData> TheStats;    // List of RunStatsData.
    static bool               Initialized; // init()'s been called?
};

#endif /*_RUNSTATS_H_*/
