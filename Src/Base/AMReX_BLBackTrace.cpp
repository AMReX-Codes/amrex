#include <iostream>
#include <sstream>
#include <fstream>

#include <unistd.h>

#include <AMReX_BLBackTrace.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX.H>

namespace amrex {

#ifdef BL_BACKTRACING
std::stack<std::pair<std::string, std::string> >  BLBackTrace::bt_stack;
#endif

void
BLBackTrace::handler(int s)
{
    signal(s, SIG_DFL);

    switch (s) {
    case SIGSEGV:
	amrex::ErrorStream() << "Segfault\n";
	break;
    case SIGFPE:
	amrex::ErrorStream() << "Erroneous arithmetic operation\n";
	break;
    case SIGTERM:
	amrex::ErrorStream() << "SIGTERM\n";
	break;
    case SIGINT:
	amrex::ErrorStream() << "SIGINT\n";
	break;
    case SIGABRT:
	amrex::ErrorStream() << "SIGABRT\n";
	break;
    }

#if defined(__linux__) && !defined(__NEC__)

    std::string errfilename;
    {
	std::ostringstream ss;
	ss << "Backtrace." << ParallelDescriptor::MyProc();
#ifdef _OPENMP
 	ss << "." << omp_get_thread_num();
#endif
	errfilename = ss.str();
    }

    if (FILE* p = fopen(errfilename.c_str(), "w")) {
	BLBackTrace::print_backtrace_info(p);
	fclose(p);
    }
    
    amrex::ErrorStream() << "See " << errfilename << " file for details" << std::endl;

#ifdef BL_BACKTRACING
    if (!bt_stack.empty()) {
	std::ofstream errfile;
	errfile.open(errfilename.c_str(), std::ofstream::out | std::ofstream::app);
	if (errfile.is_open()) {
	    errfile << std::endl;
	    while (!bt_stack.empty()) {
		errfile << "== BACKTRACE == " << bt_stack.top().first
			<<", " << bt_stack.top().second << "\n";
		bt_stack.pop();
	    }
	    errfile << std::endl;
	}
    }
#endif

    if (ParallelDescriptor::NProcs() > 1)
	sleep(3);

#endif // __linux__

    ParallelDescriptor::Abort(s, false);
}

#if defined(__linux__) && !defined(__NEC__)
void
BLBackTrace::print_backtrace_info (const std::string& filename)
{
    if (FILE* p = fopen(filename.c_str(), "w"))
    {
      BLBackTrace::print_backtrace_info(p);
      fclose(p);
    }
    else
    {
        amrex::Print() << "Warning @ BLBackTrace::print_backtrace_info: " 
                       << filename << " is not a valid output file."
                       << std::endl;
    }
}

void
BLBackTrace::print_backtrace_info (FILE* f)
{
    const int nbuf = 32;
    char **strings = NULL;
    void *buffer[nbuf];
    int nptrs = backtrace(buffer, nbuf);
    strings = backtrace_symbols(buffer, nptrs);
    if (strings != NULL) {
	int have_addr2line = 0;
	std::string cmd = "/usr/bin/addr2line";
	if (FILE *fp = fopen(cmd.c_str(), "r")) {
	    fclose(fp);
	    have_addr2line = 1;
	}
	cmd += " -Cpfie " + amrex::system::exename;

	fprintf(f, "=== If no file names and line numbers are shown below, one can run\n");
	fprintf(f, "            addr2line -Cpfie my_exefile my_line_address\n");
	fprintf(f, "    to convert `my_line_address` (e.g., 0x4a6b) into file name and line number.\n\n");
	fprintf(f, "=== Please note that the line number reported by addr2line may not be accurate.\n");
	fprintf(f, "    One can use\n");
	fprintf(f, "            readelf -wl my_exefile | grep my_line_address'\n");
	fprintf(f, "    to find out the offset for that line.\n\n");

	for (int i = 0; i < nptrs; ++i) {
	    std::string line = strings[i];
	    line += "\n";
#if !defined(_OPENMP) || !defined(__INTEL_COMPILER)
	    if (amrex::system::call_addr2line && have_addr2line && !amrex::system::exename.empty()) {
		std::size_t found1 = line.rfind('[');
		std::size_t found2 = line.rfind(']');
		if (found1 != std::string::npos && found2 != std::string::npos) {
		    std::string addr = line.substr(found1+1, found2-found1-1);
		    std::string full_cmd = cmd + " " + addr;
		    if (FILE * ps = popen(full_cmd.c_str(), "r")) {
			char buff[512];
			while (fgets(buff, sizeof(buff), ps)) {
			    line += "    ";
			    line += buff;
			}
			pclose(ps);
		    }
		}
	    }
#endif
	    fprintf(f, "%2d: %s\n", i, line.c_str());
	}
    }
}
#endif

#ifdef BL_BACKTRACING

BLBTer::BLBTer(const std::string& s, const char* file, int line)
{
    std::ostringstream ss;
    ss << "Line " << line << ", File " << file;
    line_file = ss.str();
    
#ifdef _OPENMP
    if (omp_in_parallel()) {
	std::ostringstream ss0;
	ss0 << "Proc. " << ParallelDescriptor::MyProc()
	    << ", Thread " << omp_get_thread_num()
	    << ": \"" << s << "\"";
	BLBackTrace::bt_stack.push(std::make_pair(ss0.str(), line_file));
    }
    else {
        #pragma omp parallel
	{
	    std::ostringstream ss0;
	    ss0 << "Proc. " << ParallelDescriptor::MyProc()
		<< ", Master Thread"
		<< ": \"" << s << "\"";
	    BLBackTrace::bt_stack.push(std::make_pair(ss0.str(), line_file));
	}
    }
#else
    std::ostringstream ss0;
    ss0 << "Proc. " << ParallelDescriptor::MyProc()
	<< ": \"" << s << "\"";
    BLBackTrace::bt_stack.push(std::make_pair(ss0.str(), line_file));
#endif    
}

BLBTer::~BLBTer()
{
#ifdef _OPENMP
    if (omp_in_parallel()) {
	pop_bt_stack();
    }
    else {
        #pragma omp parallel
	{
	    pop_bt_stack();
	}	
    }
#else
    pop_bt_stack();
#endif
}

void
BLBTer::pop_bt_stack()
{
    if (!BLBackTrace::bt_stack.empty()) {
	if (BLBackTrace::bt_stack.top().second.compare(line_file) == 0) {
	    BLBackTrace::bt_stack.pop();
	}
    }
}

#endif

}
