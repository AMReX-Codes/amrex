
//
// $Id: RunStats.cpp,v 1.7 1997-11-23 00:10:56 lijewski Exp $
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

double RunStats::TotalCPU;

double RunStats::TotalWCT;

Array<long> RunStats::TheCells;

List<RunStatsData> RunStats::TheStats;

RunStats::RunStats (const char* _name,
                    int         _level)
    :
    name(_name),
    level(_level)
{
    gentry = find(_name, -1);
    entry  = find(_name, _level);
    entry->is_on  = true;
}

RunStats::~RunStats () {}

void
RunStats::addCells (int  lev,
                    long count)
{
    if (lev >= TheCells.length())
    {
	TheCells.resize(lev+1);
	TheCells[lev] = 0;
    }
    TheCells[lev] += count;
}

RunStatsData *
RunStats::find (const char* _name,
                int         _level)
{
    
    for (ListIterator<RunStatsData> it(RunStats::TheStats); it; ++it)
	if (it().level == _level && it().name == _name)
	    return &RunStats::TheStats[it];

    RunStats::TheStats.append(RunStatsData(_name, _level));

    return &RunStats::TheStats.lastElement();
}

void
RunStats::init ()
{
    ParmParse pp("RunStats");

    if (pp.contains("statvar"))
    {
        aString nm;

	for (int i = 0, n = pp.countval("statvar"); i < n; i++)
        {
	    pp.get("statvar", nm, i);
	    turnOn(nm.c_str());
	}
    }
}

void
RunStats::Print (ostream&            os,
                 const RunStatsData& data,
                 double              tot_run_time,
                 double              tot_run_wtime)
{
    os << "State " << data.name << " time:  " << data.run_time << "  (";

    if (tot_run_time)
    {
       os << 100*(data.run_time/tot_run_time);
    }
    else
    {
        os << "XXXXXX";
    }
 
    os << "%)  " << data.run_wtime << "  (";

    if (tot_run_wtime)
    {
        os << 100*(data.run_wtime/tot_run_wtime);
    }
    else
    {
        os << "XXXXXX";
    }

    os << ")%\n";
}

void
RunStats::report (ostream& os)
{
    double rtime  = Utility::second();
    double rwtime = ParallelDescriptor::second();

    ParallelDescriptor::ReduceRealSum(rtime);
    ParallelDescriptor::ReduceRealMax(rwtime);

    double tot_run_time  = RunStats::TotalCPU + rtime;
    double tot_run_wtime = RunStats::TotalWCT + rwtime;
    //
    // Make a copy of the local RunStats::TheStats and reduce it.
    //
    List<RunStatsData> TheTotals = RunStats::TheStats;

    for (ListIterator<RunStatsData> it(TheTotals); it; ++it)
    {
        ParallelDescriptor::ReduceRealSum(TheTotals[it].run_time);
        ParallelDescriptor::ReduceRealMax(TheTotals[it].run_wtime);
    }

    if (ParallelDescriptor::IOProcessor())
    {
	long tot_cells = 0;
	for (int i = 0; i < TheCells.length(); i++)
        {
	    os << "Number of cells advanced at level "
               << i
               << " = "
	       << TheCells[i]
               << '\n';
	    tot_cells += TheCells[i];
	}
	os << "Total cells advanced = " << tot_cells << '\n';

	int maxlev = 0;
        ListIterator<RunStatsData> it(TheTotals);
	for ( ; it; ++it)
        {
	    maxlev = Max(maxlev, it().level);
        }

	os.setf(ios::showpoint);
        os << setprecision(4);

	for (int lev = 0; lev <= maxlev; ++lev)
        {
	    os << "\nTimings for level " << lev << ":\n\n";
	    it.rewind();
	    for ( ; it; ++it)
            {
		if (it().level == lev)
                {
                    ListIterator<RunStatsData> iti(TheTotals);
		    for ( ; iti; ++iti)
			if (iti().name == it().name && iti().level == -1)
			    break;
		    if (it().is_on)
                    {
                        Print(os, it(), tot_run_time, tot_run_wtime);
		    }
		}
	    }
	}
        os << "\nTotals for all levels:\n\n";

	it.rewind();

	for ( ; it; ++it)
        {
	    if (it().level == -1 && it().is_on == true)
            {
                Print(os, it(), tot_run_time, tot_run_wtime);
	    }
	}
        os << setprecision(8);
        os << '\n'
           << "Total CPU time        = " << tot_run_time  << '\n'
           << "Total Wall Clock time = " << tot_run_wtime << '\n';

        if (ParallelDescriptor::NProcs() > 1 && tot_run_wtime)
        {
            os << "The Parallel speedup  = "
               << tot_run_time/tot_run_wtime
               << '\n';
        }
    }
}

void
RunStats::dumpStats (ofstream& os)
{
    double rtime  = Utility::second();
    double rwtime = ParallelDescriptor::second();

    ParallelDescriptor::ReduceRealSum(rtime);
    ParallelDescriptor::ReduceRealMax(rwtime);
    //
    // Make a copy of the local RunStats::TheStats and sum the run_time's.
    //
    List<RunStatsData> TheTotals = RunStats::TheStats;

    for (ListIterator<RunStatsData> it(TheTotals); it; ++it)
    {
        ParallelDescriptor::ReduceRealSum(TheTotals[it].run_time);
    }

    if (ParallelDescriptor::IOProcessor())
    {
	os << "(ListRunStats "
           << TheTotals.length()            << '\n'
           << (RunStats::TotalCPU + rtime)  << '\n'
           << (RunStats::TotalWCT + rwtime) << '\n';

	for (ListIterator<RunStatsData> it(TheTotals); it; ++it)
        {
	    os << it();
        }
	int nlev = TheCells.length();
        os << nlev;
	for (int i = 0; i < nlev; i++)
	    os << ' ' << TheCells[i];
	os << ")\n";
    }
}

void
RunStats::readStats (ifstream& is,
                     bool      restart)
{
    RunStats::TheStats.clear();
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
	RunStats::TheStats.append(rd);
    }
    int nlev;
    is >> nlev;
    TheCells.resize(nlev);
    int i; 
    for (i = 0; i < nlev; i++)
    {
	is >> TheCells[i];
    }
    is.ignore(BL_IGNORE_MAX,')');
}

ostream &
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
