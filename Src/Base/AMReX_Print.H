#ifndef AMREX_PRINT_H_
#define AMREX_PRINT_H_

#include <sstream>
#include <AMReX_ParallelDescriptor.H>

namespace amrex
{

    // Examples:
    //   (1) print on I/O Processor of default communicator
    //       Print() << " x = " << x << std::endl;  
    //
    //   (2) print on all processors of default communicator
    //       Print(Print::AllProcs) << " x = " << x << std::endl;
    //
    //   (3) print on process rank_ of communicator comm_
    //       Print(rank_, comm_) << " x = " << x << std::endl;
    //
    //   (4) print on all processes of communicator comm_
    //       Print(Print::AllProcs, comm_) << " x = " << x << std::endl;
    
    class Print
    {
    public:

	static constexpr int AllProcs = -1;

	Print (std::ostream& os_ = std::cout) 
	    : rank(ParallelDescriptor::IOProcessorNumber())
	    , comm(ParallelDescriptor::Communicator())
	    , os(os_)
	    {}

	Print (int rank_, std::ostream& os_ = std::cout) 
	    : rank(rank_)
	    , comm(ParallelDescriptor::Communicator())
	    , os(os_)
	    {}

	Print (int rank_, MPI_Comm comm_, std::ostream& os_ = std::cout)
	    : rank(rank_)
	    , comm(comm_)
	    , os(os_)
	    {}

	~Print () {
	    if (rank == AllProcs || rank == ParallelDescriptor::MyProc(comm)) {
		os << ss.str();
	    }
	}

	template <typename T>
	Print& operator<< (const T& x) {
	    ss << x;
	    return *this;
	}

    private:
	int rank;
	MPI_Comm comm;
	std::ostream& os;
	std::ostringstream ss;
    };
}

#endif
