//BL_COPYRIGHT_NOTICE

//
// $Id: RunStats.cpp,v 1.18 1999-07-01 16:36:00 lijewski Exp $
//

#include <Utility.H>
#include <Misc.H>
#include <ParmParse.H>
#include <RunStats.H>

#ifdef BL_USE_NEW_HFILES
using std::setw;
using std::ios;
using std::setprecision;
#endif

Real               RunStats::TotalCPU;
Real               RunStats::TotalWCT;
Real               RunStats::DiskBytes;
Array<long>        RunStats::TheCells;
Array<Real>        RunStats::TheNumPts;
List<RunStatsData> RunStats::TheStats;
bool               RunStats::Initialized = false;

void
RunStats::init ()
{
    RunStats::Initialized = true;

    ParmParse pp("RunStats");

    if (pp.contains("statvar"))
    {
        aString nm;

        for (int i = 0, n = pp.countval("statvar"); i < n; i++)
        {
            pp.get("statvar", nm, i);

            RunStats::find(nm,-1)->is_on = true;
        }
    }
}

RunStatsData*
RunStats::find (const aString& _name,
                int            _level)
{
    for (ListIterator<RunStatsData> it(RunStats::TheStats); it; ++it)
    {
        if (it().level == _level && it().name == _name)
        {
            return &RunStats::TheStats[it];
        }
    }
    RunStats::TheStats.append(RunStatsData(_name, _level));

    return &RunStats::TheStats.lastElement();
}

RunStats::RunStats (const aString& _name,
                    int            _level)
    :
    name(_name),
    level(_level)
{
    if (!RunStats::Initialized)
        RunStats::init();

    gentry = RunStats::find(_name,-1);
    entry  = RunStats::find(_name,_level);

    entry->is_on = true;

    time = wtime = 0;
}

RunStats::~RunStats () {}

void
RunStats::start ()
{
    if (isOn())
    {
        time  = -Utility::second();
        wtime = -Utility::wsecond();
    }
}

void
RunStats::end ()
{
    if (isOn())
    {
        time              += Utility::second();
        wtime             += Utility::wsecond();
        entry->run_time   += time;
        gentry->run_time  += time;
        entry->run_wtime  += wtime;
        gentry->run_wtime += wtime;
    }
}

void
RunStats::addCells (int  lev,
                    long count)
{
    if (lev >= RunStats::TheCells.length())
        RunStats::TheCells.resize(lev+1, 0);
    RunStats::TheCells[lev] += count;
}

//
// Holds incremental increase to RunStats::DiskBytes on this CPU.
//
static Real Incremental_Byte_Count;

void
RunStats::addBytes (long count)
{
    Incremental_Byte_Count += count;
}

//
// Holds incremental increase in RunStats::TheNumPts on this CPU.
//
static Real Incremental_Num_Pts;

void
RunStats::addNumPts (long count)
{
    Incremental_Num_Pts += count;
}

void
RunStats::Print (ostream&            os,
                 const RunStatsData& data,
                 Real                tot_run_time,
                 Real                tot_run_wtime)
{
    int old_prec = os.precision(4);

    os << "State "
       << data.name
       << " time:  "
       << data.run_time
       << "  ("
       << (tot_run_time?100*(data.run_time/tot_run_time):-1)
       << "%)  "
       << data.run_wtime
       << "  ("
       << (tot_run_wtime?(100*(data.run_wtime/tot_run_wtime)):-1) << ")%\n";

    os.precision(old_prec);
}

void
RunStats::report (ostream& os)
{
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real rtime  = Utility::second();
    Real rwtime = Utility::wsecond();

    ParallelDescriptor::ReduceRealSum(rtime,IOProc);
    ParallelDescriptor::ReduceRealMax(rwtime,IOProc);
    ParallelDescriptor::ReduceRealSum(Incremental_Byte_Count,IOProc);

    RunStats::DiskBytes += Incremental_Byte_Count;

    Incremental_Byte_Count = 0;
    //
    // Make a copy of the local RunStats::TheStats and reduce it to IOProc.
    //
    List<RunStatsData> TheTotals = RunStats::TheStats;

    for (ListIterator<RunStatsData> it(TheTotals); it; ++it)
    {
        ParallelDescriptor::ReduceRealSum(TheTotals[it].run_time, IOProc);
        ParallelDescriptor::ReduceRealMax(TheTotals[it].run_wtime,IOProc);
    }

    RunStats::CollectNumPts();

    if (ParallelDescriptor::IOProcessor())
    {
        const Real tot_run_time  = RunStats::TotalCPU + rtime;
        const Real tot_run_wtime = RunStats::TotalWCT + rwtime;

        os.setf(ios::showpoint);

        int old_prec = os.precision(15);

        long tot_cells = 0;
        for (int i = 0; i < RunStats::TheCells.length(); i++)
        {
            os << "Number of cells advanced at level "
               << i
               << ": "
               << RunStats::TheCells[i]
               << '\n';
            tot_cells += RunStats::TheCells[i];
        }
        os << "\nTotal cells advanced: " << tot_cells << '\n';
        os << "\nTotal bytes written to disk: " << RunStats::DiskBytes << '\n';

        Real tot_num_pts = 0;

        for (int i = 0; i < RunStats::TheNumPts.length(); i++)
            tot_num_pts += RunStats::TheNumPts[i];

        if (tot_num_pts > 0)
        {
            os << '\n';
            for (int i = 0; i < RunStats::TheNumPts.length(); i++)
            {
                if (RunStats::TheNumPts[i])
                {
                    os << "Percentage of FABs allocated on CPU #"
                       << i
                       << ":\t"
                       << 100*RunStats::TheNumPts[i]/tot_num_pts
                       << '\n';
                }
            }
        }

        int maxlev = 0;
        ListIterator<RunStatsData> it(TheTotals);
        for ( ; it; ++it)
            maxlev = Max(maxlev, it().level);

        for (int lev = 0; lev <= maxlev; ++lev)
        {
            os << "\nTimings for level " << lev << " ...\n\n";
            it.rewind();
            for ( ; it; ++it)
                if (it().level == lev && it().is_on)
                    Print(os,it(),tot_run_time,tot_run_wtime);
        }
        os << "\nTotals for all levels ...\n\n";

        it.rewind();

        for ( ; it; ++it)
            if (it().level == -1 && it().is_on )
                Print(os,it(),tot_run_time,tot_run_wtime);

        os << '\n'
           << "Total CPU time        : " << tot_run_time  << '\n'
           << "Total Wall Clock time : " << tot_run_wtime << '\n';

        if (ParallelDescriptor::NProcs() > 1 && tot_run_wtime)
        {
            os << "\nThe Parallel speedup  : "
               << tot_run_time/tot_run_wtime
               << '\n';
        }

        os.precision(old_prec);
    }
}

void
RunStats::report_names (Array<aString>& stat_names)
{
    //
    // Compute the number of distinct statistics.
    //
    int nstat = 0;
    
    for (ListIterator<RunStatsData> ldi(RunStats::TheStats); ldi; ++ldi)
        if (ldi().level == -1)
            nstat++;

    stat_names.resize(nstat);
    //
    // Get their names.
    //
    int ind = 0;

    for (ListIterator<RunStatsData> ldi(RunStats::TheStats); ldi; ++ldi)
    {
        if (ldi().level == -1)
        {
            stat_names[ind] = ldi().name;
            cout << " found statistic = " << stat_names[ind] << '\n';
            ind++;
        }
    }

    BL_ASSERT(ind == nstat);
}

void
RunStats::report_values (const Array<aString>& stat_names,
                         Array<Real>&          stat_time,
                         Array<Real>&          stat_wtime,
                         Real&                 tot_run_time,
                         Real&                 tot_run_wtime,
                         long&                 tot_cells)
{
    const int NStats = stat_names.length();

    stat_time.resize(NStats);
    stat_wtime.resize(NStats);

    Real rtime  = Utility::second();
    Real rwtime = Utility::wsecond();

    ParallelDescriptor::ReduceRealSum(rtime);
    ParallelDescriptor::ReduceRealMax(rwtime);

    tot_run_time  = RunStats::TotalCPU + rtime;
    tot_run_wtime = RunStats::TotalWCT + rwtime;
    //
    // Make a copy of the local RunStats::TheStats and reduce it.
    //
    List<RunStatsData> TheTotals = RunStats::TheStats;

    for (ListIterator<RunStatsData> it(TheTotals); it; ++it)
    {
        ParallelDescriptor::ReduceRealSum(TheTotals[it].run_time);
        ParallelDescriptor::ReduceRealMax(TheTotals[it].run_wtime);
    }

    tot_cells = 0;
    for (int i = 0; i < RunStats::TheCells.length(); i++)
    {
        tot_cells += RunStats::TheCells[i];
    }

    for (int i = 0; i < NStats ; i++)
    {
        RunStatsData* e = RunStats::find(stat_names[i], -1);
        BL_ASSERT(e != 0);
        stat_time[i]  = e->run_time;
        stat_wtime[i] = e->run_wtime;
    }
}

void
RunStats::dumpStats (ofstream& os)
{
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real rtime  = Utility::second();
    Real rwtime = Utility::wsecond();

    ParallelDescriptor::ReduceRealSum(rtime,IOProc);
    ParallelDescriptor::ReduceRealMax(rwtime,IOProc);
    ParallelDescriptor::ReduceRealSum(Incremental_Byte_Count,IOProc);

    RunStats::DiskBytes += Incremental_Byte_Count;

    Incremental_Byte_Count = 0;
    //
    // Make a copy of the local RunStats::TheStats and reduce to IOProc.
    //
    List<RunStatsData> TheTotals = RunStats::TheStats;

    for (ListIterator<RunStatsData> it(TheTotals); it; ++it)
    {
        ParallelDescriptor::ReduceRealSum(TheTotals[it].run_time, IOProc);
        ParallelDescriptor::ReduceRealMax(TheTotals[it].run_wtime,IOProc);
    }

    RunStats::CollectNumPts();

    if (ParallelDescriptor::IOProcessor())
    {
        os << "(ListRunStats "
           << TheTotals.length()            << '\n'
           << (RunStats::TotalCPU + rtime)  << '\n'
           << (RunStats::TotalWCT + rwtime) << '\n'
           << RunStats::DiskBytes           << '\n';

        for (ListIterator<RunStatsData> it(TheTotals); it; ++it)
            os << it();

        long nlev = RunStats::TheNumPts.length();

        os << nlev;

        for (int i = 0; i < nlev; i++)
            os << ' ' << RunStats::TheNumPts[i];

        os << '\n';

        nlev = RunStats::TheCells.length();

        os << nlev;

        for (int i = 0; i < nlev; i++)
            os << ' ' << RunStats::TheCells[i];

        os << ")\n";
    }
}

void
RunStats::readStats (ifstream& is,
                     bool      restart)
{
    is.ignore(BL_IGNORE_MAX,'(');
    aString s;
    is >> s;
    if (s != "ListRunStats")
    {
        cerr << "unexpected token " << s << '\n';
        BoxLib::Error();
    }
    int n;
    is >> n;
    is >> RunStats::TotalCPU;
    is >> RunStats::TotalWCT;
    is >> RunStats::DiskBytes;
    while (n--)
    {
        RunStatsData rd;
        is >> rd;

        if (restart && ParallelDescriptor::NProcs() > 1)
        {
            //
            // We divide the run_time by NProcs() so that the sum
            // over all the processors comes out correct.
            //
            rd.run_time /= ParallelDescriptor::NProcs();
        }

        RunStatsData* d = find(rd.name,rd.level);

        d->run_time  += rd.run_time;
        d->run_wtime += rd.run_wtime;
    }
    long nlen;
    is >> nlen;
    RunStats::TheNumPts.resize(nlen);
    for (int i = 0; i < nlen; i++)
        is >> RunStats::TheNumPts[i];
    is >> nlen;
    RunStats::TheCells.resize(nlen);
    for (int i = 0; i < nlen; i++)
        is >> RunStats::TheCells[i];
    is.ignore(BL_IGNORE_MAX,')');
}

void
RunStats::CollectNumPts ()
{
    //
    // Gather the processor specific Incremental_Num_Pts to IOProcessor().
    //
    Array<Real> numpts(ParallelDescriptor::NProcs());

    ParallelDescriptor::Gather(&Incremental_Num_Pts,
                               1,
                               numpts.dataPtr(),
                               ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        if (ParallelDescriptor::NProcs() > RunStats::TheNumPts.length())
            //
            // Never decrease the size of RunStats::TheNumPts.
            //
            RunStats::TheNumPts.resize(ParallelDescriptor::NProcs(), 0);

        for (int i = 0; i < numpts.length(); i++)
            RunStats::TheNumPts[i] += numpts[i];
    }

    Incremental_Num_Pts = 0;
}

ostream&
operator<< (ostream&            os,
            const RunStatsData& rd)
{
    os << "(RunStatsData "
       << rd.name      << ' '
       << rd.level     << ' '
       << rd.is_on     << ' '
       << rd.run_time  << ' '
       << rd.run_wtime << ")\n";
    return os;
}

istream &
operator>> (istream&      is,
            RunStatsData& rd)
{
    is.ignore(BL_IGNORE_MAX, '(');
    aString s;
    is >> s;
    if (s != "RunStatsData")
    {
        cerr << "unexpected token " << s << '\n';
        BoxLib::Abort();
    }
    is >> rd.name;
    is >> rd.level;
    is >> rd.is_on;
    is >> rd.run_time;
    is >> rd.run_wtime;
    is.ignore(BL_IGNORE_MAX,')');
    return is;
}

ostream&
operator<< (ostream&        os,
            const RunStats& r)
{
    os << "(RunStats "
       << r.name
       << " level "
       << r.level
       << (r.isOn() ? "" : "in")
       << "active)"
       << '\n';
    return os;
}
