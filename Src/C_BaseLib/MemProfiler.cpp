
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <fstream>

#ifdef __linux
#include <sys/sysinfo.h>
#endif

#include <MemProfiler.H>
#include <ParallelDescriptor.H>
#include <BoxLib.H>
#include <ParmParse.H>

void 
MemProfiler::add(const std::string& name, std::function<MemInfo()>&& f)
{
    MemProfiler& mprofiler = getInstance();
    auto it = std::find(mprofiler.the_names.begin(), mprofiler.the_names.end(), name);
    if (it != mprofiler.the_names.end()) {
        std::string s = "MemProfiler::add failed because " + name + " already existed";
        BoxLib::Abort(s.c_str());
    }
    mprofiler.the_names.push_back(name);
    mprofiler.the_funcs.push_back(std::forward<std::function<MemInfo()> >(f));
}

MemProfiler& 
MemProfiler::getInstance ()
{
    static MemProfiler the_instance;
    return the_instance;
}

void
MemProfiler::report (const std::string& prefix)
{
    static std::string memory_log_name;
    if (memory_log_name.empty()) {
	ParmParse pp("boxlib");
	pp.query("memory_log", memory_log_name);
	if (memory_log_name.empty())
	    memory_log_name = "memlog";
    }

    getInstance().report_(prefix, memory_log_name);
}

void
MemProfiler::report_ (const std::string& prefix, const std::string& memory_log_name) const
{
    std::vector<long> cur_min;
    std::vector<long> hwm_min;
    for (auto&& f: the_funcs) {
	const MemInfo& minfo = f();
	cur_min.push_back(minfo.current_bytes);
	hwm_min.push_back(minfo.hwm_bytes);
    }

    std::vector<long> cur_max = cur_min;
    std::vector<long> hwm_max = hwm_min;

#ifdef __linux
    const int N = 5;
#else
    const int N = 1;
#endif

    std::vector<long> mymin(N, 0L);
    std::vector<long> mymax(N, 0L);

    mymin[0] = mymax[0] = std::accumulate(cur_min.begin(), cur_min.end(), 0L);

#ifdef __linux
    int ierr_sysinfo;
    int ierr_cached;
    {
	struct sysinfo info;
	ierr_sysinfo = sysinfo(&info);
	if (ierr_sysinfo == 0) {
	    mymin[1] = mymax[1] = info.mem_unit * info.totalram;
	    mymin[2] = mymax[2] = info.mem_unit * info.freeram;
	    mymin[3] = mymax[3] = info.mem_unit * (info.bufferram + info.freeram);
	    mymin[4] = mymax[4] = info.mem_unit * info.sharedram;
	}

	// Have to read /proc/meminfo to get Cached memory.
	{
	    std::ifstream ifs("/proc/meminfo");
	    std::string token;
	    long cached;
	    std::string unit;
	    while (ifs >> token) {
		if (token == "Cached:") {
		    ifs >> cached >> unit;
		    break;
		}
		ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	    }
	    if (unit == "kB") {
		ierr_cached = 0;
		mymin[3] += cached*1024L;
		mymax[3] = mymin[3];
	    } else {
		ierr_cached = -1;
	    }
	}
    }
#endif

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceLongMin(&cur_min[0], cur_min.size(), IOProc);
    ParallelDescriptor::ReduceLongMax(&cur_max[0], cur_max.size(), IOProc);
    ParallelDescriptor::ReduceLongMin(&hwm_min[0], hwm_min.size(), IOProc);
    ParallelDescriptor::ReduceLongMax(&hwm_max[0], hwm_max.size(), IOProc);
    ParallelDescriptor::ReduceLongMin(&mymin[0], N, IOProc);
    ParallelDescriptor::ReduceLongMax(&mymax[0], N, IOProc);

    if (ParallelDescriptor::IOProcessor()) {

	std::ofstream memlog(memory_log_name.c_str(), 
			     std::ofstream::out|std::ofstream::app);
	if (!memlog.good()) return;

	static int width_name = 0;
	if (width_name == 0) {
	    for (auto& x: the_names)
		width_name = std::max(width_name, int(x.size()));
	}
	const int width_bytes = 18;

	const std::string dash_name(width_name,'-');
	const std::string dash_bytes(width_bytes,'-');
	const std::string ident1(3,' ');
	const std::string ident2(6,' ');

	if (!prefix.empty())
	    memlog << prefix << " ";

	memlog << "Memory Profile Report Across Processes:\n";

	memlog << std::setfill(' ');
	memlog << ident2;
	memlog << "| " << std::setw(width_name) << std::left << "Name" << " | "
	       << std::setw(width_bytes) << std::right << "Current     " << " | "
	       << std::setw(width_bytes) << "High Water Mark " << " |\n";
	std::setw(0);

	memlog << ident2;
	memlog << "|-" << dash_name << "-+-" << dash_bytes << "-+-" << dash_bytes << "-|\n";

	for (int i = 0; i < the_names.size(); ++i) {
	    memlog << ident2;
	    memlog << "| " << std::setw(width_name) << std::left << the_names[i] << " | ";
	    memlog << Bytes{cur_min[i],cur_max[i]} << " | ";
	    memlog << Bytes{hwm_min[i],hwm_max[i]} << " |\n";
	}

	memlog << ident2;
	memlog << "|-" << dash_name << "-+-" << dash_bytes << "-+-" << dash_bytes << "-|\n";

	memlog << ident2;
	memlog << "| " << std::setw(width_name) << std::left << "Total" << " | ";
	memlog << Bytes{mymin[0],mymax[0]} << " | " << std::setw(width_bytes) << " " << " |\n";
	memlog << std::setw(0);

#ifdef __linux
	if (ierr_sysinfo == 0) {
	    memlog << ident1 << std::setw(width_bytes) << std::left << "Node total"
		   << "   " << std::setw(width_bytes) << "free";
	    if (ierr_cached == 0)
		memlog << "   " << std::setw(width_bytes) << "free+buffers+cached" << "  ";
	    else
		memlog << "   " << std::setw(width_bytes) << "free+buffers" << "   ";
	    memlog << std::setw(width_bytes) << "shared" << "\n";
	    memlog << ident1 << "[" << Bytes{mymin[1], mymax[1]} << "]"
		   << " " << "[" << Bytes{mymin[2], mymax[2]} << "]"
		   << " " << "[" << Bytes{mymin[3], mymax[3]} << "]"
		   << " " << "[" << Bytes{mymin[4], mymax[4]} << "]"
		   << "\n";
	}
#endif

	memlog << std::endl;

	memlog.close();
    }
}

std::ostream& 
operator<< (std::ostream& os, const MemProfiler::Bytes& bytes)
{
    static const long G = 1024L*1024L*1024L;
    static const long M = 1024L*1024L;
    static const long K = 1024L;

    long fac;
    std::string unit;
    if (bytes.mn >= 10L*G) {
	fac  =  G; 
	unit = "GB";
    } else if (bytes.mn >= 10L*M) {
	fac  =  M; 
	unit = "MB";
    } else if (bytes.mn >= 10L*K) {
	fac  =  K; 
	unit = "KB";
    } else {
	fac  = 1L; 
	unit = "B";
    }

    os << std::setw(5) << std::right << bytes.mn/fac << " ... "
       << std::setw(5) << std::left  << bytes.mx/fac << " " 
       << std::setw(2) << std::right << unit;
    os << std::setw(0);
    return os;
}
