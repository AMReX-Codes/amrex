//
// $Id: BLProfiler.cpp,v 1.8 2001-07-21 17:41:43 car Exp $
//

#include <winstd.H>

#ifdef __GNUC__
#include <cstdio>
#endif
#include <cctype>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <vector>
#include <stack>
#include <limits>

#include <Profiler.H>
#include <Thread.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>

struct timer_packet
{
    enum sort_criterion
    {
	sort_time,
	sort_self,
	sort_count,
	sort_name,
	sort_avg_time,
	sort_var_time,
	sort_max_time,
	sort_min_time
    };
    timer_packet()
	: count(0), time(0), self(0), avg_time(0), var_time(0),
	  max_time(std::numeric_limits<double>::min()),
	  min_time(std::numeric_limits<double>::max())
    {}
    timer_packet(int count_, double time_, double self_, double avg_time_, double var_time_, double max_time_, double min_time_)
	: count(count_), time(time_), self(self_), avg_time(avg_time_), var_time(var_time_), max_time(max_time_), min_time(min_time_)
    {}
    bool operator<(const timer_packet& r) const
    {
	if ( increasing )
	{
	    switch ( sort_by )
	    {
	    case sort_time     : return time     > r.time;
	    case sort_self     : return self     > r.self;
	    case sort_count    : return count    > r.count;
	    case sort_avg_time : return avg_time > r.avg_time;
	    case sort_var_time : return var_time > r.var_time;
	    case sort_max_time : return max_time > r.max_time;
	    case sort_min_time : return min_time > r.min_time;
	    case sort_name     : return name     < r.name;
	    default            : return self     > r.self;
	    }
	}
	else
	{
	    switch ( sort_by )
	    {
	    case sort_time     : return time     < r.time;
	    case sort_self     : return self     < r.self;
	    case sort_count    : return count    < r.count;
	    case sort_avg_time : return avg_time < r.avg_time;
	    case sort_var_time : return var_time < r.var_time;
	    case sort_max_time : return max_time < r.max_time;
	    case sort_min_time : return min_time < r.min_time;
	    case sort_name     : return name     > r.name;
	    default            : return self     < r.self;
	    }
	}
    }
    timer_packet& operator+=(const timer_packet& tp)
    {
	if ( name == "" )
	{
	    name = tp.name;
	}
	assert( name == tp.name );
	count   += tp.count;
	time    += tp.time;
	self    += tp.self;
	max_time = std::max(max_time, tp.max_time);
	min_time = std::min(min_time, tp.min_time);
	avg_time = time/count;
	var_time = 0.0;
	return *this;
    }
    struct by_name
    {
	bool operator()(const timer_packet& l, const timer_packet& r)
	{
	    return l.name < r.name;
	}
    };
private:
    static sort_criterion sort_by;
    static bool increasing;
public:
    int count;
    double time;
    double self;
    double avg_time;
    double var_time;
    double max_time;
    double min_time;
    std::string name;		// wont be messaged
};

bool timer_packet::increasing = true;
timer_packet::sort_criterion timer_packet::sort_by = timer_packet::sort_self;

#ifdef BL_USE_MPI
template <> MPI_Datatype ParallelDescriptor::Mpi_typemap<timer_packet>::type()
{
    static MPI_Datatype mine(MPI_DATATYPE_NULL);
    if ( mine == MPI_DATATYPE_NULL )
    {
	timer_packet tp[2];	// Used to construct the data types
	MPI_Datatype types[] = {
	    MPI_LB,
	    MPI_INT,
	    MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
	    MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
	    MPI_UB};
	int blocklens[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1};
	MPI_Aint disp[9];
	int n = 0;
	BL_MPI_REQUIRE( MPI_Address(&tp[0],          &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&tp[0].count,    &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&tp[0].time,     &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&tp[0].self,     &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&tp[0].avg_time, &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&tp[0].var_time, &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&tp[0].max_time, &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&tp[0].min_time, &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&tp[1],          &disp[n++]) );
	for ( int i = n-1; i >= 0; i-- )
	{
	    disp[i] -= disp[0];
	}
	assert( n == sizeof(types)/sizeof(MPI_Datatype) );
	assert( n == sizeof(disp)/sizeof(MPI_Aint) );
	assert( n == sizeof(blocklens)/sizeof(int) );
	BL_MPI_REQUIRE( MPI_Type_struct(n, blocklens, disp, types, &mine) );
	BL_MPI_REQUIRE( MPI_Type_commit( &mine ) );
    }
    return mine;
}
#endif
namespace
{
std::map<std::string, int> tagr;
std::map<int, std::string> inv_tagr;
Mutex profiling_mutex;
}

class ThreadTimerNode
{
    friend std::ostream& operator<<(std::ostream& os, const ThreadTimerNode& tree);
public:
    typedef std::map<std::string, ThreadTimerNode*>::iterator iterator;
    typedef std::map<std::string, ThreadTimerNode*>::const_iterator const_iterator;
    ThreadTimerNode()
	: wtimer(new BoxLib::WallTimer),
	  avg_time(0), var_time(0), max_time(0), min_time(0)
    {}
    ThreadTimerNode* operator[](const std::string& str);
    ThreadTimerNode* find(const std::string& str);
    const_iterator begin() const
    {
	return cntr.begin();
    }
    const_iterator end() const
    {
	return cntr.end();
    }
    void stop();
    void start()
    {
	wtimer->start();
    }
    double accum_time() const;
    int accum_count() const;
    double self_time() const;
    double max_timer() const
    {
	return max_time;
    }
    double min_timer() const
    {
	return min_time;
    }
    double child_time() const;
    void print(std::ostream& os, const std::string& str, int level);
    void mma_print(std::ostream& os, const std::string& str, int level);
    const std::map<std::string, ThreadTimerNode*>& nodes() const;
private:
    ThreadTimerNode(const ThreadTimerNode&);
    ThreadTimerNode& operator=(const ThreadTimerNode&);
    BoxLib::WallTimer* wtimer;
    double avg_time;
    double var_time;
    double max_time;
    double min_time;
    std::map<std::string, ThreadTimerNode*> cntr;
};

inline
void
ThreadTimerNode::stop()
{
    wtimer->stop();
    double X = wtimer->time();
    int N = wtimer->count();
    if ( N == 1 )
    {
	avg_time = X;
	max_time = X;
	min_time = X;
	var_time = 0.0;
    }
    else
    {
	max_time = std::max(X, max_time);
	min_time = std::min(X, min_time);
	double oavg = avg_time;
	avg_time = ( X + avg_time*(N-1))/N;
	var_time = ( var_time*(N-2) + X*X - N*avg_time*avg_time + (N-1)*oavg*oavg )/ (N - 1);
	if ( var_time < 0.0 ) var_time = 0.0;
    }
}

const std::map<std::string, ThreadTimerNode*>&
ThreadTimerNode::nodes() const
{
    return cntr;
}

double
ThreadTimerNode::accum_time() const
{
    return wtimer->accum_time();
}

int
ThreadTimerNode::accum_count() const
{
    return wtimer->count();
}

double
ThreadTimerNode::self_time() const
{
    return accum_time() - child_time();
}

double
ThreadTimerNode::child_time() const
{
    double t = 0;
    for ( const_iterator it = begin(); it != end(); ++it )
    {
	t += it->second->wtimer->accum_time();
    }
    return t;
}

inline
ThreadTimerNode*
ThreadTimerNode::operator[](const std::string& str)
{
    iterator it = cntr.find(str);
    if ( it == cntr.end() )
    {
	return cntr[str] = new ThreadTimerNode;
    }
    else
    {
	return it->second;
    }
}

inline
ThreadTimerNode*
ThreadTimerNode::find(const std::string& str)
{
    iterator it = cntr.find(str);
    if ( it == cntr.end() )
    {
	return cntr[str] = new ThreadTimerNode;
    }
    else
    {
	return it->second;
    }
}

class ThreadTimerTree
{
public:
    ThreadTimerTree()
    {
	current = new ThreadTimerNode;
    }
    ThreadTimerNode* push(const std::string& str)
    {
	ThreadTimerNode* deeper = current->find(str);
	ttd_stack.push(current);
	current = deeper;
	BL_ASSERT(current);
	return current;
    }
    ThreadTimerNode* pop()
    {
	ThreadTimerNode* t = current;
	BL_ASSERT(ttd_stack.size() > 0);
	current = ttd_stack.top();
	ttd_stack.pop();
	BL_ASSERT(t);
	return t;
    }
    ThreadTimerNode* head() const
    {
	return current;
    }
    void print(std::ostream& os) const;
    void mma_print(std::ostream& os) const;
private:
    ThreadTimerTree(const ThreadTimerTree&);
    ThreadTimerTree& operator=(const ThreadTimerTree&);
    ThreadTimerNode* current;
    std::stack<ThreadTimerNode*> ttd_stack;
};

namespace
{
const int spcmult = 2;
const int name_field_width = 40;
const int name_field_bumper = 2;
const int name_field_ellipse = 3;

const int timer_field_width = 24;
const int time_field_width  = 14;
const int count_field_width =  8;

void
spacer(std::ostream& os, int n, char c = ' ')
{
    for ( int i = 0; i < n; ++i )
    {
	os << c;
    }
}
void show_name_field(std::ostream&, int level, const std::string& str);
void show_time_count(std::ostream&, int wd, double acc, int cnt);
void show_time(std::ostream& os, double time, int scale = 1000);
void show_count(std::ostream& os, int count);
void aggregate_field_title(std::ostream& os);
}

struct ttn_packet
{
    std::string name;
    int count;
    double time;
    double self;
    double max_time;
    double min_time;
    double avg_time;
    double var_time;
    ThreadTimerNode* ttn;
    struct by_name
    {
	bool operator()(const ttn_packet& l, const ttn_packet& r)
	{
	    return l.name < r.name;
	}
    };
    struct by_self
    {
	bool operator()(const ttn_packet& l, const ttn_packet& r)
	{
	    return l.self > r.self;
	}
    };
};

void
ThreadTimerNode::print(std::ostream& os, const std::string& str, int level)
{
    show_time_count(os, timer_field_width, wtimer->accum_time(), wtimer->count());
    show_time(os, self_time()); os << ' ';
    show_time(os, max_time); os << ' ';
    show_time(os, min_time); os << ' ';
    show_time(os, avg_time); os << ' ';
    show_time(os, std::sqrt(var_time), 1000000); os << "  ";
    show_name_field(os, level, str);
    os << '\n';
    std::vector<ttn_packet> ttns;
    for ( const_iterator it = begin(); it != end(); ++it )
    {
	ttn_packet ttn;
	ttn.name = it->first;
	ttn.ttn = it->second;
	ttn.count = it->second->wtimer->count();
	ttn.time = it->second->wtimer->accum_time();
	ttn.self = it->second->self_time();
	ttn.max_time = it->second->max_time;
	ttn.min_time = it->second->min_time;
	ttn.avg_time = it->second->avg_time;
	ttn.var_time = it->second->var_time;
	ttns.push_back(ttn);
    }
#if 0
    for ( const_iterator it = begin(); it != end(); ++it )
    {
	it->second->print(os, it->first, level+1);
    }
#else
    std::sort(ttns.begin(), ttns.end(), ttn_packet::by_self());
    for ( std::vector<ttn_packet>::const_iterator it = ttns.begin(); it != ttns.end(); ++it )
    {
	it->ttn->print(os, it->name, level+1);
    }
#endif
}

void
ThreadTimerTree::print(std::ostream& os) const
{
    for ( std::map<std::string, ThreadTimerNode*>::const_iterator it = head()->begin(); it != head()->end(); ++it )
    {
	it->second->print(os, it->first, 0);
    }
}

void
ThreadTimerNode::mma_print(std::ostream& os, const std::string& str, int level)
{
    os << "\n";
    spacer(os, 2*(1+level), ' ');
    os << "bl3TimerNode[" << "\"" << str << "\"" << ", "
       << wtimer->accum_time() << ", "
       << wtimer->count() << ", "
       << self_time() << ", "
       << max_time << ", "
       << min_time << ", "
       << avg_time << ", "
       << std::sqrt(var_time) << ", ";
    os << "{";
    for ( const_iterator it = begin(); it != end(); )
    {
	it->second->mma_print(os, it->first, level+1);
	if ( ++it != end() ) os << ", ";
    }
    os << "}";
    os << "]";
}

void
ThreadTimerTree::mma_print(std::ostream& os) const
{
    os << "{";
    for ( std::map<std::string, ThreadTimerNode*>::const_iterator it = head()->begin(); it != head()->end(); )
    {
	it->second->mma_print(os, it->first, 0);
	if ( ++it != head()->end() ) os << ", ";
    }
    os << "}";
}

void
grovel(const ThreadTimerNode* nodes, const std::string& str, timer_packet& t)
{
    for ( ThreadTimerNode::const_iterator it = nodes->begin(); it != nodes->end(); ++it )
    {
	if ( it->first == str )
	{
	    t.time    += it->second->accum_time();
	    t.count   += it->second->accum_count();
	    t.self    += it->second->self_time();
	    t.max_time = std::max(it->second->max_timer(), t.max_time);
	    t.min_time = std::min(it->second->min_timer(), t.min_time);
	    t.avg_time = t.time/t.count;
	}
	grovel(it->second, str, t);
    }
}

ThreadSpecificData<int> tt_i;
Mutex tt_mutex;
std::vector<ThreadTimerTree*> tt_data;

bool Profiler::profiling = true;
int Profiler::Tag::next_itag = 0;

Profiler::Tag::Tag(const std::string& tag_)
    : tag(tag_)
{
    if ( is_profiling() )
    {
	Lock<Mutex> lock(profiling_mutex);
	if ( tagr.find(tag) != tagr.end() )
	{
	    BoxLib::Error("name already registred: ");
	}
	itag = next_itag++;
	tagr[tag] = itag;
	inv_tagr[itag] = tag;
    }
}


const std::string&
Profiler::Tag::name() const
{
    return tag;
}

bool Profiler::initialized = false;

namespace
{
std::string filename("bl3_prof");
}

void
Profiler::Initialize(int& argc, char**& argv)
{
    if ( initialized )
    {
	return;
    }
    initialized = true;
    ParmParse pp("profiler");
}

std::ostream&
operator<<(std::ostream& os, const timer_packet& tp)
{
    show_count(os, tp.count); os << ' ';
    show_time(os, tp.time); os << ' ';
    show_time(os, tp.self); os << ' ';
    show_time(os, tp.max_time); os << ' ';
    show_time(os, tp.min_time); os << ' ';
    show_time(os, tp.avg_time); os << ' ';
    show_name_field(os, 0, tp.name);
    return os;
}

std::string
Profiler::clean_name(const std::string& str)
{
#ifdef __GNUC__
    std::string result;
    unsigned int i = 0;
    unsigned cnt = 0;
    while ( i < str.length() )
    {
	if ( !isdigit(str[i]) ) break;
	cnt += (str[i] - '0') + cnt*10;
	i++;
    }
    for (; i < str.length(); ++i)
    {
	result += str[i];
    }
    return result;
#else
    return str;
#endif
}

void
Profiler::Finalize()
{
    // Try to measure overhead:
    for ( int i = 0; i < 100; ++i )
    {
	BL_PROFILE("Profiler::Finalize():load");
    }
    if ( profiling ) off();

    glean();
}

Profiler::Profiler(const Tag& tag_, bool hold)
    : tag(tag_), started(false)
{
    int* a = tt_i.get();
    if ( a == 0 )
    {
	tt_i.set(a = new int);
	Lock<Mutex> lock(tt_mutex);
	*a = tt_data.size();
	tt_data.push_back(new ThreadTimerTree);
    }
    if ( !hold )
    {
	start();
    }
}

Profiler::~Profiler()
{
    if ( started ) stop();
}

void
Profiler::start()
{
    assert( !started );
    if ( profiling ) tt_data[*tt_i.get()]->push(tag.name())->start();
    started = true;
}

void
Profiler::stop()
{
    assert( started );
    if ( profiling ) tt_data[*tt_i.get()]->pop()->stop();
    started = false;
}

#ifdef WIN32
void
Profiler::glean()
{}
#else
void
Profiler::glean()
{
    // for all threads on this processor build a cumulative timer structure using grovel.
    std::vector<timer_packet> tps;
    for ( std::map<std::string, int>::const_iterator jt = tagr.begin(); jt != tagr.end(); ++jt )
    {
	timer_packet t;
	t.name = jt->first;
	for ( std::vector<ThreadTimerTree*>::const_iterator it = tt_data.begin(); it != tt_data.end(); ++it)
	{
	    grovel((*it)->head(), jt->first, t);
	}
	tps.push_back(t);
    }

    std::vector< std::vector< timer_packet> > t_packets;
    if ( ParallelDescriptor::IOProcessor() )
    {
	t_packets.resize(ParallelDescriptor::NProcs());
	t_packets[0] = tps;
    }

    std::vector<size_t> ntags = ParallelDescriptor::Gather(tagr.size(), ParallelDescriptor::IOProcessor());
    for ( int i = 1; i < ParallelDescriptor::NProcs(); ++i )
    {
	if ( ParallelDescriptor::IOProcessor() )
	{
	    std::vector<size_t> lngths(ntags[i]);
	    ParallelDescriptor::Recv(lngths, i, 101);
	    std::vector< std::string > strngs;
	    for ( unsigned int j = 0; j < lngths.size(); ++j )
	    {
		std::vector<char> a(lngths[j]);
		ParallelDescriptor::Recv(a, i, 102);
		strngs.push_back(std::string(a.begin(), a.end()));
	    }
	    t_packets[i].resize(ntags[i]);
	    ParallelDescriptor::Recv(t_packets[i], i, 103);
	    for ( unsigned int j = 0; j < ntags[i]; ++j )
	    {
		t_packets[i][j].name = strngs[j];
	    }
	}
	else if ( ParallelDescriptor::MyProc() == i )
	{
	    std::vector<size_t> lngths;
	    for ( std::vector<timer_packet>::const_iterator it = tps.begin(); it != tps.end(); ++it )
	    {
		lngths.push_back(it->name.size());
	    }
	    ParallelDescriptor::Send(lngths, 0, 101);
	    for ( std::vector<timer_packet >::const_iterator it = tps.begin(); it != tps.end(); ++it )
	    {
		const char* name = it->name.c_str();
		ParallelDescriptor::Send(name, it->name.size(), 0, 102);
	    }
	    ParallelDescriptor::Send(tps, 0, 103);
	}
    }

    if ( ParallelDescriptor::IOProcessor() )
    {
	std::vector< timer_packet > tp_total;
	for ( int i = 0; i < ParallelDescriptor::NProcs(); ++i )
	{
	    std::copy(t_packets[i].begin(), t_packets[i].end(), std::back_inserter(tp_total));
	}
	std::sort(tp_total.begin(), tp_total.end(), timer_packet::by_name());
	std::vector<timer_packet> tp_summary;
	timer_packet t;
	for ( std::vector<timer_packet>::iterator it = tp_total.begin(); it != tp_total.end(); ++it )
	{
	    t += *it;
	    std::vector<timer_packet>::const_iterator it1(it); ++it1;
	    if ( it1 == tp_total.end() || it->name != it1->name )
	    {
		// finish the current one
		tp_summary.push_back(t);
		t = timer_packet();
		continue;
	    }
	}
	std::sort(tp_summary.begin(), tp_summary.end());

	// SUmmary
	std::ofstream os(filename.c_str());
	if ( !os )
	{
	    std::cerr << "filename = " << filename << std::endl;
	    BoxLib::Error("failed to open prof file");
	}

	os << "------------------------------------------------------------------------\n\n";
	os << "Profiling report\n\n";
	os << "------------------------------------------------------------------------\n\n";
	os << "Timer resolution is "; show_time(os, BoxLib::WallTimer::tick(), 1000000); os << " (us)\n";
	os << "Number of Processors: " << ParallelDescriptor::NProcs() << std::endl;

	spacer(os,  2, '\n');
	spacer(os, 72, '-'); os << '\n';
	spacer(os, 72, '-'); os << '\n';
	os << "Aggregate report\n\n";
	spacer(os, 16, '-'); os << '\n';
	aggregate_field_title(os);
	std::copy(tp_summary.begin(), tp_summary.end(), std::ostream_iterator<timer_packet>(os, "\n"));

	spacer(os,  2, '\n');
	spacer(os, 72, '-'); os << '\n';
	spacer(os, 72, '-'); os << '\n';
	os << "Per-Processor Report" << '\n';
	spacer(os, 20, '-'); os << '\n';
	os << "Number of processes Reporting " << t_packets.size() << std::endl;
	for ( int i = 0; i < ParallelDescriptor::NProcs(); ++i )
	{
	    std::vector< timer_packet > tps( t_packets[i] );
	    std::sort(tps.begin(), tps.end());
	    os << "\nProcessor :" << i << std::endl;
	    spacer(os, 20, '-'); os << '\n';
	    aggregate_field_title(os);
	    std::copy(tps.begin(), tps.end(), std::ostream_iterator<timer_packet>(os, "\n"));
	}

	spacer(os,  2, '\n');
	spacer(os, 72, '-'); os << '\n';
	os << "Details Profiling report\n";
	spacer(os, 72, '-'); os << '\n';
    }

    for ( int i = 0; i < ParallelDescriptor::NProcs(); ++i )
    {
	if ( i == ParallelDescriptor::MyProc() )
	{
	    std::ofstream os(filename.c_str(), std::ios::app);
	    os << "\nProcessor Number " << std::setw(4) << i << std::endl;
	    os <<   "-----------------" "----\n";
	    os << "\nNumber of threads = " << std::setw(4) << tt_data.size() << std::endl;
	    os <<   "--------------------" "----\n";
	    int cnt = 0;
	    for ( std::vector<ThreadTimerTree*>::const_iterator it = tt_data.begin(); it != tt_data.end(); ++it)
	    {
		os << "\n\n";
		os << "Thread " << std::setw(4) << cnt++ << std::endl;
		os << "-------" "----\n";
		std::ios::fmtflags fldadjust = os.setf(std::ios::right, std::ios::adjustfield);
		os << std::setw(timer_field_width) << "Time (ms)/  Count " << ' ';
		os << std::setw(time_field_width) << "Self (ms)" << ' ';
		os << std::setw(time_field_width) << "Max (ms)" << ' ';
		os << std::setw(time_field_width) << "Min (ms)" << ' ';
		os << std::setw(time_field_width) << "Avg (ms)" << ' ';
		os << std::setw(time_field_width) << "STD (us)" << ' ';
		os.setf(std::ios::left, std::ios::adjustfield);
		os << std::setw(name_field_width) << "Timer: ";
		os << std::endl << std::flush;
		os.setf(fldadjust, std::ios::adjustfield);
		spacer(os, name_field_width + timer_field_width + 5*(time_field_width+1), '-'); os << '\n';
		(*it)->print(os);
		os << '\n';
		aggregate_field_title(os);
		std::vector<timer_packet> tps;
		for ( std::map<std::string, int>::const_iterator jt = tagr.begin(); jt != tagr.end(); ++jt )
		{
		    timer_packet t;
		    t.name = jt->first;
		    grovel((*it)->head(), jt->first, t);
		    tps.push_back(t);
		}
		std::sort(tps.begin(), tps.end());
		for ( std::vector<timer_packet>::const_iterator jt = tps.begin(); jt != tps.end(); ++jt )
		{
		    if ( jt->time != 0 )
		    {
			os << *jt << std::endl;
		    }
		}
	    }
	    os << std::endl;
	}
	ParallelDescriptor::Barrier();
    }
    // MMA dump
    std::string mma_fname = filename + ".m";
    if ( ParallelDescriptor::IOProcessor() )
    {
	std::ofstream os(mma_fname.c_str());
	os << "{";
    }
    for ( int i = 0; i < ParallelDescriptor::NProcs(); ++i )
    {
	if ( i == ParallelDescriptor::MyProc() )
	{
	    std::ofstream os(mma_fname.c_str(), std::ios::app);
	    std::ios::fmtflags oldFlags = os.flags();
	    os.setf(std::ios::fixed, std::ios::floatfield);
	    int oprec = os.precision(8);
	    os << "{";
	    for ( std::vector<ThreadTimerTree*>::const_iterator it = tt_data.begin(); it != tt_data.end(); )
	    {
		(*it)->mma_print(os);
		if ( ++it != tt_data.end() ) os << ", ";
	    }
	    os.precision(oprec);
	    os.flags(oldFlags);
	    os << "}\n";
	    if ( i != ParallelDescriptor::NProcs()-1 ) os << ", ";
	}
	ParallelDescriptor::Barrier();
    }
    if ( ParallelDescriptor::IOProcessor() )
    {
	std::ofstream os(mma_fname.c_str(), std::ios::app);
	os << "}";
    }
}
#endif

bool
Profiler::is_profiling()
{
    Lock<Mutex> lock(profiling_mutex);
    return profiling;
}

void
Profiler::on()
{
    Lock<Mutex> lock(profiling_mutex);
    profiling = true;
}

void
Profiler::off()
{
    Lock<Mutex> lock(profiling_mutex);
    profiling = false;
}


namespace
{

void
show_time_count(std::ostream& os, int wd, double time, int count)
{
    std::ios::fmtflags oldFlags = os.flags();
    if ( wd >= 3 )
    {
	wd -= 3;
    }
    int fwd = 2*wd/3;
    int oprec = os.precision();
    int prec = fwd-1-5-3;
    if ( prec < 0 ) prec = oprec;
#ifdef __GNUC__
    char buff[128];
    sprintf(buff, "[%*.*f/%*d]", 2*wd/3, prec, time*1000, 1*wd/3, count);
    os << buff;
#else
    os << std::fixed;
    os << std::showpoint;
    os << std::setprecision(prec);
    os << std::setw(0)	// Reset the width to zero..."
       << "[" << std::setw(2*wd/3) << time*1000
       << "/" << std::setw(1*wd/3) << count << "]";

    //
    os << std::setprecision(oprec);
#endif
    os.flags(oldFlags);
}

void
show_time(std::ostream& os, double time, int scale)
{
    int precred  = 0;
    int tscale = scale;
    while ( tscale /= 10 )
    {
	precred++;
    }
    std::ios::fmtflags oldFlags = os.flags();
    int fwd = time_field_width;
    int oprec = os.precision();
    int prec = fwd-1-5-precred;
    if ( prec < 0 ) prec = oprec;
#ifdef __GNUC__
    char buff[1024];
    sprintf(buff, "%*.*f", fwd, prec, time*scale);
    os << buff;
#else
    os << std::fixed;
    os << std::showpoint;
    os << std::setprecision(prec) << std::setw(fwd) << time*scale;
    os << std::setprecision(oprec);
#endif
    os.flags(oldFlags);
}

void
show_count(std::ostream& os, int count)
{
    os << std::setw(count_field_width) << count;
}

void
show_name_field(std::ostream& os, int level, const std::string& str)
{
    spacer(os, level*spcmult);
    os << str;
    if ( false )
    {
	int len = str.size();
	int wdth = name_field_width - level*spcmult - name_field_bumper;
	if ( len > wdth )
	{
	    int plen = len - wdth + name_field_ellipse;
	    spacer(os, name_field_ellipse, '.');
	    os << std::setw(wdth-name_field_ellipse) << str.substr(plen);
	}
	else
	{
	    os << std::setw(wdth) << str.substr(0, wdth);
	}
	spacer(os, name_field_bumper);
    }
}

void
aggregate_field_title(std::ostream& os)
{
    std::ios::fmtflags fldadjust = os.setf(std::ios::right, std::ios::adjustfield);
    os << std::setw(count_field_width) << "Count" << ' ';
    os << std::setw(time_field_width) << "Total (ms)" << ' ';
    os << std::setw(time_field_width) << "Self (ms)" << ' ';
    os << std::setw(time_field_width) << "Max (ms)" << ' ';
    os << std::setw(time_field_width) << "Min (ms)" << ' ';
    os << std::setw(time_field_width) << "Avg (ms)" << ' ';
    os.setf(std::ios::left, std::ios::adjustfield);
    os << std::setw(name_field_width) << "Registered profilers:";
    os.setf(fldadjust, std::ios::adjustfield);
    os << '\n';
    spacer(os, name_field_width + count_field_width + 5*(time_field_width+1), '-'); os << '\n';
}

}
