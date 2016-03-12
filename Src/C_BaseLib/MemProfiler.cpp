
#include <numeric>
#include <algorithm>

#include <unistd.h>

#include <MemProfiler.H>
#include <ParallelDescriptor.H>
#include <BoxLib.H>

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
    getInstance().report_(prefix);
}

void
MemProfiler::report_ (const std::string& prefix) const
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
	if (!prefix.empty())
	    std::cout << prefix << " ";
	std::cout << "Memory Profile Report Across Processes:\n";
	for (int i = 0; i < the_names.size(); ++i) {
	    std::cout << "      " << the_names[i] << " ::"
		      << "  current : " << Bytes{cur_min[i],cur_max[i]}
	              << "  high water mark : " << Bytes{hwm_min[i],hwm_max[i]}
	              << "\n";
	}
	std::cout << "   Process Uses : " << Bytes{mymin[0],mymax[0]} << "\n"
	          << "   Node Free    : " << Bytes{mymin[1],mymax[1]} << "\n"
		  << "   Node Total   : " << Bytes{mymin[2],mymax[2]};
	std::cout << std::endl;
    }
}

std::ostream& 
operator<< (std::ostream& os, const MemProfiler::Bytes& bytes)
{
    static const long GB = 10L*1024L*1024L*1024L;
    static const long MB = 10L*1024L*1024L;
    static const long KB = 10L*1024L;

    long fac;
    std::string unit;
    if (bytes.mn >= 10L*GB) {
	fac  =  GB; 
	unit = "GB";
    } else if (bytes.mn >= 10L*MB) {
	fac  =  MB; 
	unit = "MB";
    } else if (bytes.mn >= 10L*KB) {
	fac  =  KB; 
	unit = "KB";
    } else {
	fac  = 1L; 
	unit = "B";
    }

    os << "[" << bytes.mn/fac << " ... " << bytes.mx/fac << " " << unit << "]";
    if (os.fail())
        BoxLib::Error("operator<<(ostream&,const MemProfiler::Bytes&) failed");
    return os;
}
