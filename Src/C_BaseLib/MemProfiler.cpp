
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <fstream>

#include <unistd.h>

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

    long mymin[3], mymax[3];
    mymin[0] = mymax[0] = std::accumulate(cur_min.begin(), cur_min.end(), 0L);

    static const long page_size = sysconf(_SC_PAGESIZE);
    long sys_free  = sysconf(_SC_AVPHYS_PAGES) * page_size;
    long sys_total = sysconf(_SC_PHYS_PAGES) * page_size;
    mymin[1] = mymax[1] = sys_free;
    mymin[2] = mymax[2] = sys_total;

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceLongMin(&cur_min[0], cur_min.size(), IOProc);
    ParallelDescriptor::ReduceLongMax(&cur_max[0], cur_max.size(), IOProc);
    ParallelDescriptor::ReduceLongMin(&hwm_min[0], hwm_min.size(), IOProc);
    ParallelDescriptor::ReduceLongMax(&hwm_max[0], hwm_max.size(), IOProc);
    ParallelDescriptor::ReduceLongMin(mymin, 3, IOProc);
    ParallelDescriptor::ReduceLongMax(mymax, 3, IOProc);

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

	memlog << ident1;
	memlog << "Node Free  : " << "[" << Bytes{mymin[1],mymax[1]} << "]\n";
	memlog << ident1;
	memlog << "Node Total : " << "[" << Bytes{mymin[2],mymax[2]} << "]\n";

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
