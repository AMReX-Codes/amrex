//BL_COPYRIGHT_NOTICE

//
// $Id: RunStats.cpp,v 1.14 1997-11-26 04:21:47 lijewski Exp $
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

double RunStats::DiskBytes;

Array<long> RunStats::TheCells;

Array<double> RunStats::TheNumPts;

List<RunStatsData> RunStats::TheStats;

bool RunStats::Initialized = false;

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
	    turnOn(nm.c_str());
	}
    }
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

RunStats::RunStats (const char* _name,
                    int         _level)
    :
    name(_name),
    level(_level)
{
    if (!RunStats::Initialized)
        RunStats::init();

    gentry = find(_name, -1);
    entry  = find(_name, _level);

    entry->is_on  = true;
}

RunStats::~RunStats () {}

void
RunStats::start ()
{
    if (gentry->is_on && entry->is_on)
    {
        time  = -Utility::second();
        wtime = -Utility::wsecond();
    }
}

void
RunStats::end ()
{
    if (gentry->is_on && entry->is_on)
    {
        time  += Utility::second();
        wtime += Utility::wsecond();

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
    {
	RunStats::TheCells.resize(lev+1, 0);
    }
    RunStats::TheCells[lev] += count;
}

//
// Holds incremental increase to RunStats::DiskBytes on this CPU.
//
static double Incremental_Byte_Count;

void
RunStats::addBytes (long count)
{
    Incremental_Byte_Count += count;
}

//
// Holds incremental increase in RunStats::TheNumPts on this CPU.
//
static double Incremental_Num_Pts;

void
RunStats::addNumPts (long count)
{
    Incremental_Num_Pts += count;
}

void
RunStats::Print (ostream&            os,
                 const RunStatsData& data,
                 double              tot_run_time,
                 double              tot_run_wtime)
{
    int old_prec = os.precision(4);

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

    os.precision(old_prec);
}

void
RunStats::report (ostream& os)
{
    double rtime  = Utility::second();
    double rwtime = Utility::wsecond();

    ParallelDescriptor::ReduceRealSum(rtime);
    ParallelDescriptor::ReduceRealMax(rwtime);

    ParallelDescriptor::ReduceRealSum(Incremental_Byte_Count);
    RunStats::DiskBytes += Incremental_Byte_Count;
    Incremental_Byte_Count = 0;

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

    RunStats::CollectNumPts();

    if (ParallelDescriptor::IOProcessor())
    {
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

        double tot_num_pts = 0;

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
        {
	    maxlev = Max(maxlev, it().level);
        }

	for (int lev = 0; lev <= maxlev; ++lev)
        {
	    os << "\nTimings for level " << lev << " ...\n\n";
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
        os << "\nTotals for all levels ...\n\n";

	it.rewind();

	for ( ; it; ++it)
        {
	    if (it().level == -1 && it().is_on == true)
            {
                Print(os, it(), tot_run_time, tot_run_wtime);
	    }
	}

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
RunStats::dumpStats (ofstream& os)
{
    double rtime  = Utility::second();
    double rwtime = Utility::wsecond();

    ParallelDescriptor::ReduceRealSum(rtime);
    ParallelDescriptor::ReduceRealMax(rwtime);

    ParallelDescriptor::ReduceRealSum(Incremental_Byte_Count);
    RunStats::DiskBytes += Incremental_Byte_Count;
    Incremental_Byte_Count = 0;
    //
    // Make a copy of the local RunStats::TheStats and sum the run_time's.
    //
    List<RunStatsData> TheTotals = RunStats::TheStats;

    for (ListIterator<RunStatsData> it(TheTotals); it; ++it)
    {
        ParallelDescriptor::ReduceRealSum(TheTotals[it].run_time);
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
        {
	    os << it();
        }
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
	RunStats::TheStats.append(rd);
    }
    long nlen;
    is >> nlen;
    RunStats::TheNumPts.resize(nlen);
    int i;
    for (i = 0; i < nlen; i++)
    {
	is >> RunStats::TheNumPts[i];
    }
    is >> nlen;
    RunStats::TheCells.resize(nlen);
    for (i = 0; i < nlen; i++)
    {
	is >> RunStats::TheCells[i];
    }
    is.ignore(BL_IGNORE_MAX,')');
}

void
RunStats::CollectNumPts ()
{
    struct
    {
        int    m_cpu; // CPU #
        double m_pts; // Total numPts() of FABs on that CPU.
    }
    msg_hdr;

    ParallelDescriptor::SetMessageHeaderSize(sizeof(msg_hdr));

    int MyProc = ParallelDescriptor::MyProc();

    if (!ParallelDescriptor::IOProcessor())
    {
        msg_hdr.m_cpu = MyProc;
        msg_hdr.m_pts = Incremental_Num_Pts;
        ParallelDescriptor::SendData(ParallelDescriptor::IOProcessor(),
                                     &msg_hdr,
                                     0,
                                     0);
    }

    ParallelDescriptor::Synchronize();

    if (ParallelDescriptor::IOProcessor())
    {
        if (ParallelDescriptor::NProcs() > RunStats::TheNumPts.length())
        {
            //
            // Never decrease the size of RunStats::TheNumPts.
            //
            RunStats::TheNumPts.resize(ParallelDescriptor::NProcs(),0);
        }
        RunStats::TheNumPts[MyProc] += Incremental_Num_Pts;

        for (int len; ParallelDescriptor::GetMessageHeader(len, &msg_hdr); )
        {
            assert(len == 0);
            RunStats::TheNumPts[msg_hdr.m_cpu] += msg_hdr.m_pts;
            ParallelDescriptor::ReceiveData(0, 0);
        }
    }

    Incremental_Num_Pts = 0;
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
