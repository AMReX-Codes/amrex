#ifndef AMREX_PRINT_H_
#define AMREX_PRINT_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ANSIEscCode.H>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <utility>

namespace amrex
{
    template <typename T>
    std::ostream& operator<< (std::ostream& os, Array<T,AMREX_SPACEDIM> const& a)
    {
        os << AMREX_D_TERM( '(' << a[0] , <<
                            ',' << a[1] , <<
                            ',' << a[2])  << ')';
        return os;
    }

    template <typename T, typename S>
    std::ostream& operator<<(std::ostream& os, const std::pair<T, S>& v)
    {
        os << "(" << v.first << ", " << v.second << ")";
        return os;
    }

    //! This class provides the user with a few print options
    class Print
    {
    public:

        static constexpr int AllProcs = -1;

    /**
    * \brief Print on I/O Processor of the default communicator
    * Example: Print() << " x = " << x << '\n';
    */
        explicit Print (std::ostream& os_ = amrex::OutStream())
            : to_print(ParallelContext::IOProcessorSub())
            , os(os_)
            { ss.precision(os.precision()); }

    /**
    * \brief Print on given rank of the default communicator
    * Example: Print(2) << " x = " << x << '\n'; // Print on rank 2
    */
        Print (int rank_, std::ostream& os_ = amrex::OutStream())
            : to_print(rank_ == AllProcs || rank_ == ParallelContext::MyProcSub())
            , os(os_)
            { ss.precision(os.precision()); }

    /**
    * \brief Print on process rank_ of communicator comm_
    * Example: Print(rank_, comm_) << " x = " << x << '\n';
    */
        Print (int rank_, MPI_Comm comm_, std::ostream& os_ = amrex::OutStream())
            : to_print((rank_ == AllProcs) || (
#ifdef AMREX_USE_MPI
                        comm_ != MPI_COMM_NULL &&
#endif
                        rank_ == ParallelDescriptor::MyProc(comm_)))
            , os(os_)
        {
            ss.precision(os.precision());
        }

        ~Print () {
            if (to_print) {
                os.flush();
                os << ss.str();
                os.flush();
            }
        }

        Print (Print const&) = delete;
        Print (Print &&) = delete;
        Print& operator= (Print const&) = delete;
        Print& operator= (Print &&) = delete;

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
        bool to_print;
        std::ostream &os;
        std::ostringstream ss;
    };

    //! Print on all processors of the default communicator
    class AllPrint
        : public Print
    {
    public:
        //! Example: AllPrint() << " x = " << x << '\n';
        explicit AllPrint (std::ostream& os_ = amrex::OutStream())
            : Print(Print::AllProcs, os_)
            {}

    };

    //! This class prints to a file with a given base name
    class PrintToFile
    {
    public:

        static constexpr int AllProcs = -1;

        explicit PrintToFile (std::string file_name_)
            : file_name(std::move(file_name_))
            , to_print(ParallelContext::IOProcessorSub())
        { Initialize(); }

        PrintToFile (std::string file_name_, int rank_ )
            : file_name(std::move(file_name_))
            , to_print(rank_ == AllProcs || rank_ == ParallelContext::MyProcSub())
        { Initialize(); }

        PrintToFile (std::string file_name_, int rank_, MPI_Comm comm_)
            : file_name(std::move(file_name_))
            , to_print((rank_ == AllProcs) || (
#ifdef AMREX_USE_MPI
                        comm_ != MPI_COMM_NULL &&
#endif
                        rank_ == ParallelDescriptor::MyProc(comm_)))
        { Initialize(); }

        ~PrintToFile () {
            if (to_print) {
                ofs.flush();
                ofs << ss.str();
                ofs.flush();
            }
        }

        PrintToFile (PrintToFile const&) = delete;
        PrintToFile (PrintToFile &&) = delete;
        PrintToFile& operator= (PrintToFile const&) = delete;
        PrintToFile& operator= (PrintToFile &&) = delete;

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
            if (to_print) {
                int my_proc_global = ParallelDescriptor::MyProc();
                std::string proc_file_name = file_name + "." + std::to_string(my_proc_global);
#ifdef AMREX_USE_OMP
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
        bool to_print;
        std::ofstream ofs;
        std::ostringstream ss;
    };

    //! Print on all processors of the default communicator
    class AllPrintToFile
        : public PrintToFile
    {
    public:
        //! Example: AllPrint() << " x = " << x << '\n';
        explicit AllPrintToFile (std::string file_name_)
            : PrintToFile(std::move(file_name_), Print::AllProcs)
            {}

    };
}

#endif
