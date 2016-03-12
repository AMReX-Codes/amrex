
#include <iostream>
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
	    std::cout << "      " << the_names[i] << " ::";
	    const std::string& uc = convertUnit(cur_min[i], cur_max[i]);
	    const std::string& uh = convertUnit(hwm_min[i], hwm_max[i]);
	    std::cout << "  current : [" << cur_min[i] << " ... " << cur_max[i] << " "<<uc<<"]";
	    std::cout << "  high water mark : [" << hwm_min[i] << " ... " << hwm_max[i] <<" "<<uh<< "]";
	    std::cout << "\n";
	}
	const std::string& u0 = convertUnit(mymin[0], mymax[0]);
	const std::string& u1 = convertUnit(mymin[1], mymax[1]);
	const std::string& u2 = convertUnit(mymin[2], mymax[2]);
	std::cout << "   Using : [" << mymin[0] << " ... " << mymax[0] << " "<<u0<<"]" << std::endl;
	std::cout << "   Avail : [" << mymin[1] << " ... " << mymax[1] << " "<<u1<<"]" << std::endl;
	std::cout << "   Total : [" << mymin[2] << " ... " << mymax[2] << " "<<u2<<"]" << std::endl;
    }
}

std::string
MemProfiler::convertUnit(long& mn, long& mx)
{
    static const long GB = 10L*1024L*1024L*1024L;
    static const long MB = 10L*1024L*1024L;
    static const long KB = 10L*1024L;

    std::string unit;

    if (mn >= 10L*GB) {
	mn /= GB;
	mx /= GB;
	unit = "GB";
    } else if (mn >= 10L*MB) {
	mn /= MB;
	mx /= MB;
	unit = "MB";
    } else if (mn >= 10L*KB) {
	mn /= KB;
	mx /= KB;
	unit = "KB";
    } else {
	unit = "B";
    }

    return unit;
}

