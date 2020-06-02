
#include <limits>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <fstream>

#ifdef __linux__
#include <unistd.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#endif

#include <AMReX_MemProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>

namespace amrex {

std::unique_ptr<MemProfiler> MemProfiler::the_instance = nullptr;

void 
MemProfiler::add (const std::string& name, std::function<MemInfo()>&& f)
{
    MemProfiler& mprofiler = getInstance();
    auto it = std::find(mprofiler.the_names.begin(), mprofiler.the_names.end(), name);
    if (it != mprofiler.the_names.end()) {
        std::string s = "MemProfiler::add (MemInfo) failed because " + name + " already existed";
        amrex::Abort(s.c_str());
    }
    mprofiler.the_names.push_back(name);
    mprofiler.the_funcs.push_back(std::move(f));
}

void 
MemProfiler::add (const std::string& name, std::function<NBuildsInfo()>&& f)
{
    MemProfiler& mprofiler = getInstance();
    auto it = std::find(mprofiler.the_names_builds.begin(), mprofiler.the_names_builds.end(), name);
    if (it != mprofiler.the_names_builds.end()) {
        std::string s = "MemProfiler::add (NBuildsInfo) failed because " + name + " already existed";
        amrex::Abort(s.c_str());
    }
    mprofiler.the_names_builds.push_back(name);
    mprofiler.the_funcs_builds.push_back(std::move(f));
}

MemProfiler& 
MemProfiler::getInstance ()
{
    if (the_instance == nullptr) {
        the_instance.reset(new MemProfiler());
    }
    return *the_instance;
}

void
MemProfiler::Finalize ()
{
    the_instance.reset();
}

void
MemProfiler::report (const std::string& prefix)
{
    static std::string memory_log_name;
    if (memory_log_name.empty()) {
	ParmParse pp("amrex");
	pp.query("memory_log", memory_log_name);
	if (memory_log_name.empty())
	    memory_log_name = "memlog";
    }

    getInstance().report_(prefix, memory_log_name);
}

void
MemProfiler::report_ (const std::string& prefix, const std::string& memory_log_name) const
{
    std::vector<Long> cur_min;
    std::vector<Long> hwm_min;
    for (auto&& f: the_funcs) {
	const MemInfo& minfo = f();
	cur_min.push_back(minfo.current_bytes);
	hwm_min.push_back(minfo.hwm_bytes);
    }
    std::vector<Long> cur_max = cur_min;
    std::vector<Long> hwm_max = hwm_min;

    std::vector<int>  num_builds_min;
    std::vector<int>  hwm_builds_min;
    for (auto&& f: the_funcs_builds) {
	const NBuildsInfo& binfo = f();
	num_builds_min.push_back(binfo.current_builds);
	hwm_builds_min.push_back(binfo.hwm_builds);
    }
    std::vector<int>  num_builds_max = num_builds_min;
    std::vector<int>  hwm_builds_max = hwm_builds_min;

#ifdef __linux__
    const int N = 9;
#else
    const int N = 1;
#endif

    std::vector<Long> mymin(N, 0L);
    std::vector<Long> mymax(N, 0L);

    mymin[0] = mymax[0] = std::accumulate(cur_min.begin(), cur_min.end(), 0L);

#ifdef __linux__
    int ierr_proc_status = 0;
    const int ipstat = 1;
    const int npstat = 4;
    {
	static pid_t pid = getpid();
	std::string fname = "/proc/"+std::to_string(pid) + "/status";
	std::ifstream ifs(fname.c_str());
	std::string token, unit;
	Long n;
	int nfound = 0;
	while (ifs >> token) {
	    if (token == "VmPeak:") {
		nfound++;
		ifs >> n >> unit;
		if (unit == "kB") {
		    mymin[ipstat+0] = mymax[ipstat+0] = n * 1024L;
		} else {
		    ierr_proc_status = -1;
		    break;
		}
	    } else if (token == "VmSize:") {
		nfound++;
		ifs >> n >> unit;
		if (unit == "kB") {
		    mymin[ipstat+1] = mymax[ipstat+1] = n * 1024L;
		} else {
		    ierr_proc_status = -1;
		    break;
		}
	    } else if (token == "VmHWM:") {
		nfound++;
		ifs >> n >> unit;
		if (unit == "kB") {
		    mymin[ipstat+2] = mymax[ipstat+2] = n * 1024L;
		} else {
		    ierr_proc_status = -1;
		    break;
		}
	    } else if (token == "VmRSS:") {
		nfound++;
		ifs >> n >> unit;
		if (unit == "kB") {
		    mymin[ipstat+3] = mymax[ipstat+3] = n * 1024L;
		} else {
		    ierr_proc_status = -1;
		    break;
		}
	    }
	    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	    if (nfound == npstat) break;
	}
    }

    int ierr_sysinfo;
    int ierr_cached;
    const int isinfo = ipstat + npstat;
    const int nsinfo = 4;
    {
	struct sysinfo info;
	ierr_sysinfo = sysinfo(&info);
	if (ierr_sysinfo == 0) {
	    mymin[isinfo+0] = mymax[isinfo+0] = info.mem_unit * info.totalram;
	    mymin[isinfo+1] = mymax[isinfo+1] = info.mem_unit * info.freeram;
	    mymin[isinfo+2] = mymax[isinfo+2] = info.mem_unit * (info.bufferram + info.freeram);
	    mymin[isinfo+3] = mymax[isinfo+3] = info.mem_unit * info.sharedram;
	}

	// Have to read /proc/meminfo to get Cached memory.
	{
	    std::ifstream ifs("/proc/meminfo");
	    std::string token, unit;
	    Long cached;
	    while (ifs >> token) {
		if (token == "Cached:") {
		    ifs >> cached >> unit;
		    break;
		}
		ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	    }
	    if (unit == "kB") {
		ierr_cached = 0;
		mymin[isinfo+2] += cached*1024L;
		mymax[isinfo+2] = mymin[isinfo+2];
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

    ParallelDescriptor::ReduceIntMin (&num_builds_min[0], num_builds_min.size(), IOProc);
    ParallelDescriptor::ReduceIntMax (&num_builds_max[0], num_builds_max.size(), IOProc);
    ParallelDescriptor::ReduceIntMin (&hwm_builds_min[0], hwm_builds_min.size(), IOProc);
    ParallelDescriptor::ReduceIntMax (&hwm_builds_max[0], hwm_builds_max.size(), IOProc);

    if (ParallelDescriptor::IOProcessor()) {

	std::ofstream memlog(memory_log_name.c_str(), 
			     std::ofstream::out|std::ofstream::app);
	if (!memlog.good()) return;

	static int width_name = 0;
	if (width_name == 0) {
	    for (auto& x: the_names)
		width_name = std::max(width_name, int(x.size()));
	    for (auto& x: the_names_builds)
		width_name = std::max(width_name, int(x.size()));
	}
	const int width_bytes = 18;

	const std::string dash_name(width_name,'-');
	const std::string dash_bytes(width_bytes,'-');
	const std::string ident(6,' ');

	if (!prefix.empty())
	    memlog << prefix << " ";

	memlog << "Memory Profile Report Across Processes:\n";

	memlog << std::setfill(' ');
	memlog << ident;
	memlog << "| " << std::setw(width_name) << std::left << "Name" << " | "
	       << std::setw(width_bytes) << std::right << "Current     " << " | "
	       << std::setw(width_bytes) << "High Water Mark " << " |\n";
	std::setw(0);

	memlog << ident;
	memlog << "|-" << dash_name << "-+-" << dash_bytes << "-+-" << dash_bytes << "-|\n";

	std::vector<int> idxs(the_names.size());
	std::iota(idxs.begin(), idxs.end(), 0);
	std::sort(idxs.begin(), idxs.end(), [&](int i, int j)
		  { return hwm_max[i] > hwm_max[j]; });

	for (int ii = 0; ii < idxs.size(); ++ii) {
	    int i = idxs[ii];
	    if (hwm_max[i] > 0) {
		memlog << ident;
		memlog << "| " << std::setw(width_name) << std::left << the_names[i] << " | ";
		memlog << Bytes{cur_min[i],cur_max[i]} << " | ";
		memlog << Bytes{hwm_min[i],hwm_max[i]} << " |\n";
	    }
	}

	memlog << ident;
	memlog << "|-" << dash_name << "-+-" << dash_bytes << "-+-" << dash_bytes << "-|\n";

	memlog << ident;
	memlog << "| " << std::setw(width_name) << std::left << "Total" << " | ";
	memlog << Bytes{mymin[0],mymax[0]} << " | " << std::setw(width_bytes) << " " << " |\n";
	memlog << std::setw(0);

	// Number of builds
	{
	    memlog << "\n";
	    memlog << ident;
	    memlog << "| " << std::setw(width_name) << std::left << "Name" << " | "
		   << std::setw(width_bytes) << std::right << "Current #    " << " | "
		   << std::setw(width_bytes) << "High Water Mark #" << " |\n";
	    std::setw(0);
	    
	    memlog << ident;
	    memlog << "|-" << dash_name << "-+-" << dash_bytes << "-+-" << dash_bytes << "-|\n";
	    
	    for (int i = 0; i < the_names_builds.size(); ++i) {
		if (hwm_builds_max[i] > 0) {
		    memlog << ident;
		    memlog << "| " << std::setw(width_name) << std::left << the_names_builds[i] << " | ";
		    memlog << Builds{num_builds_min[i],num_builds_max[i]} << " | ";
		    memlog << Builds{hwm_builds_min[i],hwm_builds_max[i]} << " |\n";
		}
	    }
	}

#ifdef __linux__
	if (ierr_proc_status == 0) {
	    memlog << "\n";
	    memlog << " * " << std::setw(width_bytes) << std::left << "Proc VmPeak"
		   << "   " << std::setw(width_bytes) << "VmSize"	    
		   << "   " << std::setw(width_bytes) << "VmHWM"	    
		   << "   " << std::setw(width_bytes) << "VmRSS" << "\n";
	    memlog << "  ";
	    for (int i = 0; i < npstat; ++i)
		memlog << " [" << Bytes{mymin[ipstat+i], mymax[ipstat+i]} << "]";
	    memlog << "\n";
	}

	if (ierr_sysinfo == 0) {
	    memlog << "\n";
	    memlog << " * " << std::setw(width_bytes) << std::left << "Node total"
		   << "   " << std::setw(width_bytes) << "free";
	    if (ierr_cached == 0)
		memlog << "   " << std::setw(width_bytes) << "free+buffers+cached" << "  ";
	    else
		memlog << "   " << std::setw(width_bytes) << "free+buffers" << "   ";
	    memlog << std::setw(width_bytes) << "shared" << "\n";
	    memlog << "  ";
	    for (int i = 0; i < nsinfo; ++i)
		memlog << " [" << Bytes{mymin[isinfo+i], mymax[isinfo+i]} << "]";
	    memlog << "\n";
	}
#endif

	memlog << std::endl;

	memlog.close();
    }
}

std::ostream& 
operator<< (std::ostream& os, const MemProfiler::Bytes& bytes)
{
    constexpr Long G = 1024L*1024L*1024L;
    constexpr Long M = 1024L*1024L;
    constexpr Long K = 1024L;

    Long fac;
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

std::ostream& 
operator<< (std::ostream& os, const MemProfiler::Builds& builds)
{
    os << std::setw(6) << std::right << builds.mn << " ... "
       << std::setw(7) << std::left  << builds.mx; 
    os << std::setw(0);
    return os;
}

}
