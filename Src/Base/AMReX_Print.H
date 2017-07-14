#ifndef AMREX_PRINT_H_
#define AMREX_PRINT_H_

#include <sstream>

#include <AMReX_ParallelDescriptor.H>

namespace amrex
{
    //! This class provides the user with a few print options
    class Print
    {
    public:

	static constexpr int AllProcs = -1;

    /**
    * \brief Print on I/O Processor of the default communicator
    * Example: Print() << " x = " << x << std::endl;  
    */
	Print (std::ostream& os_ = std::cout) 
	    : rank(ParallelDescriptor::IOProcessorNumber())
	    , comm(ParallelDescriptor::Communicator())
	    , os(os_)
	    { ss.precision(os.precision()); }

    /**
    * \brief Print on all processors of the default communicator
    * Example: Print(Print::AllProcs) << " x = " << x << std::endl;
    */
	Print (int rank_, std::ostream& os_ = std::cout) 
	    : rank(rank_)
	    , comm(ParallelDescriptor::Communicator())
	    , os(os_)
	    { ss.precision(os.precision()); }

    /**   
    * \brief Print on process rank_ of communicator comm_
    * Example: Print(rank_, comm_) << " x = " << x << std::endl;
    */
	Print (int rank_, MPI_Comm comm_, std::ostream& os_ = std::cout)
	    : rank(rank_)
	    , comm(comm_)
	    , os(os_)
	    { ss.precision(os.precision()); }

	~Print () {
	    if (rank == AllProcs || rank == ParallelDescriptor::MyProc(comm)) {
		os.flush();
		os << ss.str();
		os.flush();
	    }
	}

        Print& SetPrecision(int p) {
            ss.precision(p);
            return *this;
	}

	template <typename T>
	Print& operator<< (const T& x) {
            ss << x;
	    return *this;
	}

        Print& operator<< (        std::basic_ostream<char, std::char_traits<char> >& 
                           (*func)(std::basic_ostream<char, std::char_traits<char> >&))
        {
            ss << func;
            return *this;
        }

    private:
	int rank;
	MPI_Comm comm;
	std::ostream& os;
	std::ostringstream ss;
    };

    //! Print on all processors of the default communicator
    class AllPrint
	: public Print
    {
    public:
        //! Example: AllPrint() << " x = " << x << std::endl;
	AllPrint (std::ostream& os_ = std::cout)
	    : Print(Print::AllProcs, os_) 
	    {}
    };
}

#endif
