
//
// $Id: RunStats.cpp,v 1.5 1997-11-20 00:19:30 lijewski Exp $
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

List<RunStatsData> RunStats::ld;

double RunStats::total_run_time;

double RunStats::total_run_wtime;

Array<long> RunStats::cells;

RunStats::RunStats (const char* _name,
                    int         _level)
    :
    name(_name),
    level(_level)
{
    gentry = find(_name, -1);
    entry  = find(_name, _level);
    entry->is_on = true;
}

void
RunStats::addCells (int  lev,
                    long count)
{
    if (lev >= cells.length())
    {
	cells.resize(lev+1);
	cells[lev] = 0;
    }
    cells[lev] += count;
}

RunStatsData *
RunStats::find (const char* _name,
                int         _level)
{
    
    for (ListIterator<RunStatsData> ldi(ld); ldi; ++ldi)
	if (ldi().level == _level && ldi().name == _name)
	    return &ld[ldi];

    ld.append(RunStatsData(_name, _level));

    return &ld.lastElement();
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
	   << rd.run_wtime << ' '
	   << rd.max_time  << ' '
	   << rd.max_wtime << ")\n";
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
	is >> rd.max_time;
	is >> rd.max_wtime;
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

    double tot_run_time  = total_run_time  + rtime;
    double tot_run_wtime = total_run_wtime + rwtime;

    double inv_t_r_time  = tot_run_time == 0.0  ? 0.0 : 1.0/tot_run_time;
    double inv_t_r_wtime = tot_run_wtime == 0.0 ? 0.0 : 1.0/tot_run_wtime;

    if (ParallelDescriptor::IOProcessor())
    {
	long tot_cells = 0;
        int i; 
	for (i = 0; i < cells.length(); i++)
        {
	    os << "Number of cells advanced at level "
               << i
               << " = "
	       << cells[i]
               << '\n';
	    tot_cells += cells[i];
	}
	os << "Total cells advanced = " << tot_cells << '\n';

	int maxlev = 0;
        ListIterator<RunStatsData> ldi(ld);
	for ( ; ldi; ++ldi)
	    maxlev = Max(maxlev, ldi().level);

        int lev; 

	os.setf(ios::showpoint);
	os << setprecision(4);
	for (lev = 0; lev <= maxlev; ++lev)
        {
	    os << "\ntimings for level " << lev << ":\n\n";
	    ldi.rewind();
	    for ( ; ldi; ++ldi)
            {
		if (ldi().level == lev)
                {
		    for (ListIterator<RunStatsData> ldii(ld); ldii; ++ldii)
			if (ldii().name == ldi().name && ldii().level == -1)
			    break;
		    if (ldi().is_on)
                    {
			os << "State " << ldi().name;
			os << " time  = "
			   << setw(6) << ldi().run_time << ' '
			   << setw(6) << 100.0*(ldi().run_time*inv_t_r_time) << "% "
			   << setw(6) << ldi().run_wtime << ' '
			   << setw(6) << 100.0*(ldi().run_wtime*inv_t_r_wtime) << "% "
			   << setw(6) << ldi().max_time << ' '
			   << setw(6) << ldi().max_wtime;
			os << '\n';
		    }
		}
	    }
	}
        os << '\n';

	ldi.rewind();

	for ( ; ldi; ++ldi)
        {
	    if (ldi().level == -1 && ldi().is_on == true)
            {
                os << "total " << ldi().name << " time";
                os << "  = "
                   << setw(6) << ldi().run_time << "  "
                   << setw(6) << 100.0*(ldi().run_time*inv_t_r_time) << "% "
                   << setw(6) << ldi().run_wtime << "  "
                   << setw(6) << 100.0*(ldi().run_wtime*inv_t_r_wtime) << "% ";
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
	os << "(ListRunStats " << ld.length() << '\n'
           << rtime  + total_run_time         << '\n'
           << rwtime + total_run_wtime        << '\n';

	for (ListIterator<RunStatsData> ldi(ld); ldi; ++ldi)
	    os << ldi();
	int nlev = cells.length();
        os << nlev;
        int i; 
	for (i = 0; i < nlev; i++)
	    os << ' ' << cells[i];
	os << ")\n";
    }
}

void
RunStats::readStats (ifstream& is)
{
    ld.clear();
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
    is >> total_run_time;
    is >> total_run_wtime;
    while (n--)
    {
	RunStatsData rd;
	is >> rd;
	ld.append(rd);
    }
    int nlev;
    is >> nlev;
    cells.resize(nlev);
    int i; 
    for (i = 0; i < nlev; i++)
    {
	is >> cells[i];
    }
    is.ignore(BL_IGNORE_MAX,')');
}

void
RunStats::start ()
{
    if (gentry->is_on && entry->is_on)
    {
        ParallelDescriptor::Synchronize();
	time  = -Utility::second();
	wtime = -ParallelDescriptor::second();
    }
}

void
RunStats::pause ()
{
    if (gentry->is_on && entry->is_on)
    {
	time  += Utility::second();
	wtime += ParallelDescriptor::second();
        ParallelDescriptor::Synchronize();
    }
}

void
RunStats::resume ()
{
    if (gentry->is_on && entry->is_on)
    {
        ParallelDescriptor::Synchronize();
	time  += Utility::second();
	wtime += ParallelDescriptor::second();
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
	    entry->max_time   += time;
	    entry->max_wtime  += wtime;
	    gentry->run_time  += time;
	    gentry->run_wtime += wtime;
	    gentry->max_time  += time;
	    gentry->max_wtime += wtime;
	}
        else
        {
	    double tmp[2];
	    tmp[0] = time;
	    tmp[1] = wtime;
	    ParallelDescriptor::ReduceRealSum(tmp[0]);
	    ParallelDescriptor::ReduceRealSum(tmp[1]);
	    entry->run_time   += tmp[0];
	    entry->run_wtime  += tmp[1];
	    gentry->run_time  += tmp[0];
	    gentry->run_wtime += tmp[1];
	    tmp[0] = time;
	    tmp[1] = wtime;
	    ParallelDescriptor::ReduceRealMax(tmp[0]);
	    ParallelDescriptor::ReduceRealMax(tmp[1]);
	    entry->max_time   += tmp[0];
	    entry->max_wtime  += tmp[1];
	    gentry->max_time  += tmp[0];
	    gentry->max_wtime += tmp[1];
	}
    }
}
