
//
// $Id: RunStats.cpp,v 1.6 1997-11-21 03:50:11 lijewski Exp $
//

#include <Utility.H>
#include <Misc.H>
#include <ParmParse.H>
#include <RunStats.H>
#include <ParallelDescriptor.H>

#ifdef BL_USE_NEW_HFILES
using std::setw;
using std::ios;
using std::setprecision;
#endif

List<RunStatsData> RunStats::TheStats;

double RunStats::TotalCPU;

double RunStats::TotalWCT;

Array<long> RunStats::TheCells;

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
    
    for (ListIterator<RunStatsData> it(TheStats); it; ++it)
	if (it().level == _level && it().name == _name)
	    return &TheStats[it];

    TheStats.append(RunStatsData(_name, _level));

    return &TheStats.lastElement();
}

void
RunStats::init ()
{
    ParmParse pp("RunStats");

    if (pp.contains("statvar"))
    {
	int n = pp.countval("statvar");

        aString nm;

	for (int i = 0; i < n; i++)
        {
	    pp.get("statvar", nm, i);
	    turnOn(nm.c_str());
	}
    }
}

ostream &
operator<< (ostream&            os,
            const RunStatsData& rd)
{
	os << "(RunStatsData "
	   << rd.name      << ' '
	   << rd.level     << ' '
	   << rd.is_on     <<  ' '
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

void
RunStats::report (ostream& os)
{
    double rtime  = Utility::second();
    double rwtime = ParallelDescriptor::second();

    ParallelDescriptor::ReduceRealSum(rtime);
    ParallelDescriptor::ReduceRealMax(rwtime);

    double tot_run_time  = TotalCPU  + rtime;
    double tot_run_wtime = TotalWCT + rwtime;

    double inv_t_r_time  = tot_run_time == 0.0  ? 0.0 : 1.0/tot_run_time;
    double inv_t_r_wtime = tot_run_wtime == 0.0 ? 0.0 : 1.0/tot_run_wtime;

    if (ParallelDescriptor::IOProcessor())
    {
	long tot_cells = 0;
        int i; 
	for (i = 0; i < TheCells.length(); i++)
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
        ListIterator<RunStatsData> it(TheStats);
	for ( ; it; ++it)
	    maxlev = Max(maxlev, it().level);

        int lev; 

	os.setf(ios::showpoint);
	os << setprecision(4);
	for (lev = 0; lev <= maxlev; ++lev)
        {
	    os << "\ntimings for level " << lev << ":\n\n";
	    it.rewind();
	    for ( ; it; ++it)
            {
		if (it().level == lev)
                {
		    for (ListIterator<RunStatsData> iti(TheStats); iti; ++iti)
			if (iti().name == it().name && iti().level == -1)
			    break;
		    if (it().is_on)
                    {
			os << "State " << it().name;
			os << " time  = "
			   << setw(6) << it().run_time << ' '
			   << setw(6) << 100.0*(it().run_time*inv_t_r_time) << "% "
			   << setw(6) << it().run_wtime << ' '
			   << setw(6) << 100.0*(it().run_wtime*inv_t_r_wtime) << "%\n";
		    }
		}
	    }
	}
        os << '\n';

	it.rewind();

	for ( ; it; ++it)
        {
	    if (it().level == -1 && it().is_on == true)
            {
                os << "total " << it().name << " time";
                os << "  = "
                   << setw(6) << it().run_time << "  "
                   << setw(6) << 100.0*(it().run_time*inv_t_r_time) << "% "
                   << setw(6) << it().run_wtime << "  "
                   << setw(6) << 100.0*(it().run_wtime*inv_t_r_wtime) << "%\n";
                os << '\n';
	    }
	}
	os << "total CPU time          = " << tot_run_time  << '\n';
	os << "total Wall Clock time   = " << tot_run_wtime << '\n';
    }
}

void
RunStats::dumpStats (ofstream& os)
{
    double rtime  = Utility::second();
    double rwtime = ParallelDescriptor::second();

    ParallelDescriptor::ReduceRealSum(rtime);
    ParallelDescriptor::ReduceRealMax(rwtime);

    if (ParallelDescriptor::IOProcessor())
    {
	os << "(ListRunStats " << TheStats.length() << '\n'
           << rtime  + TotalCPU         << '\n'
           << rwtime + TotalWCT        << '\n';

	for (ListIterator<RunStatsData> it(TheStats); it; ++it)
	    os << it();
	int nlev = TheCells.length();
        os << nlev;
        int i; 
	for (i = 0; i < nlev; i++)
	    os << ' ' << TheCells[i];
	os << ")\n";
    }
}

void
RunStats::readStats (ifstream& is)
{
    TheStats.clear();
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
    is >> TotalCPU;
    is >> TotalWCT;
    while (n--)
    {
	RunStatsData rd;
	is >> rd;
	TheStats.append(rd);
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

void
RunStats::start ()
{
    if (gentry->is_on && entry->is_on)
    {
	time  = -Utility::second();
	wtime = -ParallelDescriptor::second();
    }
}

void
RunStats::end ()
{
    if (gentry->is_on && entry->is_on)
    {
	time  += Utility::second();
	wtime += ParallelDescriptor::second();

	if (ParallelDescriptor::NProcs() == 1)
        {
	    entry->run_time   += time;
	    entry->run_wtime  += wtime;
	    gentry->run_time  += time;
	    gentry->run_wtime += wtime;
	}
        else
        {
	    double tmp[2] = { time, wtime };
	    ParallelDescriptor::ReduceRealSum(tmp[0]);
	    ParallelDescriptor::ReduceRealSum(tmp[1]);
	    entry->run_time   += tmp[0];
	    entry->run_wtime  += tmp[1];
	    gentry->run_time  += tmp[0];
	    gentry->run_wtime += tmp[1];
	}
    }
}
