#ifndef AMREX_PRINT_H_
#define AMREX_PRINT_H_

#include <sstream>
#include <fstream>
#include <iomanip>

#include <AMReX.H>
#include <AMReX_ParallelContext.H>
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
	explicit Print (std::ostream& os_ = amrex::OutStream()) 
	    : rank(ParallelContext::IOProcessorNumberSub())
	    , comm(ParallelContext::CommunicatorSub())
	    , os(os_)
	    { ss.precision(os.precision()); }

    /**
    * \brief Print on all processors of the default communicator
    * Example: Print(Print::AllProcs) << " x = " << x << std::endl;
    */
	Print (int rank_, std::ostream& os_ = amrex::OutStream())
	    : rank(rank_)
	    , comm(ParallelContext::CommunicatorSub())
	    , os(os_)
	    { ss.precision(os.precision()); }

    /**   
    * \brief Print on process rank_ of communicator comm_
    * Example: Print(rank_, comm_) << " x = " << x << std::endl;
    */
	Print (int rank_, MPI_Comm comm_, std::ostream& os_ = amrex::OutStream())
	    : rank(rank_)
	    , comm(comm_)
	    , os(os_)
	    { ss.precision(os.precision()); }

	~Print () {
	    if (rank == AllProcs || rank == ParallelContext::MyProcSub()) {
                std::ostream * my_os = ParallelContext::OFSPtrSub();
                if (my_os) {
                    my_os->flush();
                    (*my_os) << ss.str();
                    my_os->flush();
                }
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
        std::ostream &os;
	std::ostringstream ss;
    };

    //! Print on all processors of the default communicator
    class AllPrint
	: public Print
    {
    public:
        //! Example: AllPrint() << " x = " << x << std::endl;
	explicit AllPrint (std::ostream& os_ = amrex::OutStream())
	    : Print(Print::AllProcs, os_) 
	    {}

    };

    //! This class prints to a file with a given base name
    class PrintToFile
    {
    public:

	static constexpr int AllProcs = -1;

	explicit PrintToFile (const std::string& file_name_) 
	    : file_name(file_name_)
            , rank(ParallelDescriptor::IOProcessorNumber())
	    , comm(ParallelDescriptor::Communicator())
        { Initialize(); }
        
	PrintToFile (const std::string& file_name_, int rank_ ) 
	    : file_name(file_name_)
	    , rank(rank_)
	    , comm(ParallelDescriptor::Communicator())
        { Initialize(); }
        
	PrintToFile (const std::string& file_name_, int rank_, MPI_Comm comm_)
	    : file_name(file_name_)
	    , rank(rank_)
	    , comm(comm_)
        { Initialize(); }

	~PrintToFile () {
	    if (rank == AllProcs || rank == ParallelDescriptor::MyProc(comm)) {
		ofs.flush();
		ofs << ss.str();
		ofs.flush();
	    }
	}

        PrintToFile& SetPrecision(int p) {
            ss.precision(p);
            return *this;
	}

	template <typename T>
	PrintToFile& operator<< (const T& x) {
            ss << x;
	    return *this;
	}

        PrintToFile& operator<< (        std::basic_ostream<char, std::char_traits<char> >& 
                                 (*func)(std::basic_ostream<char, std::char_traits<char> >&))
        {
            ss << func;
            return *this;
        }

    private:

        void Initialize() {
            int my_proc = ParallelDescriptor::MyProc(comm);
            if (rank == AllProcs || rank == my_proc) {
                std::string proc_file_name = file_name + "." + std::to_string(my_proc);
#ifdef _OPENMP
                proc_file_name += "." + std::to_string(omp_get_thread_num());
#endif
                ofs.open(proc_file_name, std::ios_base::app);
                if (!ofs.is_open()) {
                    amrex::Error("Could not open file for appending in amrex::Print()");
                }
                ss.precision(ofs.precision());
            }
        }

        std::string file_name;
	int rank;
	MPI_Comm comm;
        std::ofstream ofs;
	std::ostringstream ss;
    };

    //! Print on all processors of the default communicator
    class AllPrintToFile
	: public PrintToFile
    {
    public:
        //! Example: AllPrint() << " x = " << x << std::endl;
	explicit AllPrintToFile (const std::string& file_name_)
	    : PrintToFile(file_name_, Print::AllProcs)
	    {}

    };
}

#endif
