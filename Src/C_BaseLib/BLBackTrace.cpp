#ifdef BL_BACKTRACING

#include <iostream>
#include <sstream>

#include <unistd.h>
#include <execinfo.h>
#include <signal.h>

#include <BoxLib.H>
#include <BLBackTrace.H>
#include <ParallelDescriptor.H>

std::stack<std::pair<std::string, std::string> >  BLBackTrace::bt_stack;

void
BLBackTrace::handler(int s)
{
    const int nbuf = 16;
    void *buffer[nbuf];
    int nptrs = backtrace(buffer, nbuf);
    backtrace_symbols_fd(buffer, nptrs, STDERR_FILENO);

#ifdef _OPENMP
#pragma omp critical(print_bt_stack)
#endif
    {
	std::cout << std::endl;
	while (!bt_stack.empty()) {
	    std::cout << "== BACKTRACE == " << bt_stack.top().first
		      <<", " << bt_stack.top().second << "\n";
	    bt_stack.pop();
	}
	std::cout << std::endl;
    }

    switch (s) {
    case SIGSEGV:
	BoxLib::Abort("Segfault");
	break;
    case SIGFPE:
	BoxLib::Abort("Erroneous arithmetic operation");
	break;
    default:
	ParallelDescriptor::Abort();
    }
}

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
