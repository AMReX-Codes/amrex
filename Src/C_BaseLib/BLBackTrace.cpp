#ifdef BL_BACKTRACE

#include <iostream>
#include <sstream>

#ifdef __linux__
#include <unistd.h>
#include <execinfo.h>
#endif

#include <BoxLib.H>
#include <BLBackTrace.H>
#include <ParallelDescriptor.H>

std::stack<std::pair<std::string, std::string> >  BLBackTrace::bt_stack;

void
BLBackTrace::handler(int s)
{
#ifdef __linux__
    void *buffer[10];
    int nptrs = backtrace(buffer, 10);
    backtrace_symbols_fd(buffer, nptrs, STDERR_FILENO);
#endif

#pragma omp parallel 
    {
#pragma omp critical(print_bt_stack)
	{
	    std::cout << std::endl;
	    while (!bt_stack.empty()) {
		std::cout << "== BACKTRACE == proc. " << ParallelDescriptor::MyProc()
#ifdef _OPENMP
			  << " thread " << omp_get_thread_num() 
#endif
			  << " : " << bt_stack.top().first
			  << ", " << bt_stack.top().second << "\n";
		bt_stack.pop();
	    }
	    std::cout << std::endl;
	}
    }

    BoxLib::Abort("exiting due to segfault");
}

BLBTer::BLBTer(const std::string& s, const char* file, int line)
{
    std::ostringstream ss;
    ss << "File " << file << ", Line " << line;
    file_line = ss.str();
    BLBackTrace::bt_stack.push(std::make_pair(s,file_line));
}

BLBTer::~BLBTer()
{
    if (!BLBackTrace::bt_stack.empty()) {
	if (BLBackTrace::bt_stack.top().second.compare(file_line) == 0) {
	    BLBackTrace::bt_stack.pop();
	}
    }
}

#endif
