#ifndef AMREX_PRINT_H_
#define AMREX_PRINT_H_

#include <sstream>
#include <fstream>

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
	    , os(&os_, [](std::ostream*){})
	    { ss.precision(os->precision()); }

    /**
    * \brief Print on all processors of the default communicator
    * Example: Print(Print::AllProcs) << " x = " << x << std::endl;
    */
	Print (int rank_, std::ostream& os_ = std::cout) 
	    : rank(rank_)
	    , comm(ParallelDescriptor::Communicator())
	    , os(&os_, [](std::ostream*){})
	    { ss.precision(os->precision()); }

    /**   
    * \brief Print on process rank_ of communicator comm_
    * Example: Print(rank_, comm_) << " x = " << x << std::endl;
    */
	Print (int rank_, MPI_Comm comm_, std::ostream& os_ = std::cout)
	    : rank(rank_)
	    , comm(comm_)
	    , os(&os_, [](std::ostream*){})
	    { ss.precision(os->precision()); }

        /**
         * \brief Print to a file on process rank_ of communicator comm_.
         * If called from an OMP parallel region, each thread will output
         * to a different file.
         *
         */
        Print (const std::string& file_name_,
               int rank_ = ParallelDescriptor::IOProcessorNumber(),
               MPI_Comm comm_ = ParallelDescriptor::Communicator())
            : rank(rank_)
            , comm(comm_)
            , os(NULL, std::default_delete<std::ostream>())
        {
            int my_proc = ParallelDescriptor::MyProc(comm);
	    if (rank == AllProcs || rank == my_proc) {
                std::string proc_file_name = file_name_ + "." + std::to_string(my_proc);
#ifdef _OPENMP
                proc_file_name += "." + std::to_string(omp_get_thread_num());
#endif
                os.reset(new std::ofstream(proc_file_name, std::ios_base::app));
                if (!dynamic_cast<std::ofstream&>(*os).is_open()) {
                    amrex::Error("Could not open file for appending in amrex::Print()");
                }
                ss.precision(os->precision());
            }
        }

	~Print () {
	    if (rank == AllProcs || rank == ParallelDescriptor::MyProc(comm)) {
		os->flush();
		*os << ss.str();
		os->flush();
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
        std::unique_ptr<std::ostream, std::function<void(std::ostream*)>> os;
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

        //! All procs output to a unique file (or files if threading is on)
        AllPrint (const std::string& file_name_)
            : Print(file_name_, Print::AllProcs)
        {}
    };


}

#endif
