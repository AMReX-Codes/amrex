#ifndef AMREX_PRINT_H_
#define AMREX_PRINT_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Print.H
 * \brief MPI-aware print utility classes.
 */

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
    //! Write an AMREX_SPACEDIM-element Array to an output stream.
    template <typename T>
    std::ostream& operator<< (std::ostream& os, Array<T,AMREX_SPACEDIM> const& a)
    {
        os << AMREX_D_TERM( '(' << a[0] , <<
                            ',' << a[1] , <<
                            ',' << a[2])  << ')';
        return os;
    }

    //! Write a std::pair to an output stream as \c (first, second).
    template <typename T, typename S>
    std::ostream& operator<<(std::ostream& os, const std::pair<T, S>& v)
    {
        os << "(" << v.first << ", " << v.second << ")";
        return os;
    }

    /**
     * \class amrex::Print
     * \ingroup amrex_utilities
     * \brief Print to an output stream from a single MPI rank.
     *
     * Output is buffered and flushed to the stream only from the designated
     * rank when the object is destroyed (typically at the end of the
     * statement).
     *
     * Example:
     * \code{.cpp}
     * amrex::Print() << "x = " << x << '\n';          // I/O processor only
     * amrex::Print(2) << "rank 2 says " << x << '\n'; // rank 2 only
     * \endcode
     */
    class Print
    {
    public:

        static constexpr int AllProcs = -1;

        //! Print on the I/O processor of the default communicator.
        explicit Print (std::ostream& os_ = amrex::OutStream())
            : rank(ParallelContext::IOProcessorNumberSub())
            , comm(ParallelContext::CommunicatorSub())
            , os(os_)
            { ss.precision(os.precision()); }

        //! Print on rank \p rank_ of the default communicator.
        Print (int rank_, std::ostream& os_ = amrex::OutStream())
            : rank(rank_)
            , comm(ParallelContext::CommunicatorSub())
            , os(os_)
            { ss.precision(os.precision()); }

        //! Print on rank \p rank_ of communicator \p comm_.
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

        Print (Print const&) = delete;
        Print (Print &&) = delete;
        Print& operator= (Print const&) = delete;
        Print& operator= (Print &&) = delete;

        //! Set the output precision to \p p digits and return \c *this.
        Print& SetPrecision(int p) {
            ss.precision(p);
            return *this;
        }

        //! Append \p x to the output buffer and return \c *this.
        template <typename T>
        Print& operator<< (const T& x) {
            ss << x;
            return *this;
        }

        //! Support stream manipulators such as \c std::endl and \c std::flush.
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

    /**
     * \class amrex::AllPrint
     * \ingroup amrex_utilities
     * \brief Print to an output stream from all MPI ranks.
     *
     * Example:
     * \code{.cpp}
     * amrex::AllPrint() << "rank " << amrex::ParallelDescriptor::MyProc() << '\n';
     * \endcode
     */
    class AllPrint
        : public Print
    {
    public:
        //! Print on all processors of the default communicator.
        explicit AllPrint (std::ostream& os_ = amrex::OutStream())
            : Print(Print::AllProcs, os_)
            {}

    };

    /**
     * \class amrex::PrintToFile
     * \ingroup amrex_utilities
     * \brief Print to a per-process file from a single MPI rank.
     *
     * Output is appended to a file named \c file_name.rank (and
     * \c file_name.rank.thread when OpenMP is enabled).
     */
    class PrintToFile
    {
    public:

        static constexpr int AllProcs = -1;

        //! Print on the I/O processor of the default communicator.
        explicit PrintToFile (std::string file_name_)
            : file_name(std::move(file_name_))
            , rank(ParallelContext::IOProcessorNumberSub())
            , comm(ParallelContext::CommunicatorSub())
        { Initialize(); }

        //! Print on rank \p rank_ of the default communicator.
        PrintToFile (std::string file_name_, int rank_ )
            : file_name(std::move(file_name_))
            , rank(rank_)
            , comm(ParallelContext::CommunicatorSub())
        { Initialize(); }

        //! Print on rank \p rank_ of communicator \p comm_.
        PrintToFile (std::string file_name_, int rank_, MPI_Comm comm_)
            : file_name(std::move(file_name_))
            , rank(rank_)
            , comm(comm_)
        { Initialize(); }

        ~PrintToFile () {
            if (rank == AllProcs || rank == ParallelContext::MyProcSub()) {
                ofs.flush();
                ofs << ss.str();
                ofs.flush();
            }
        }

        PrintToFile (PrintToFile const&) = delete;
        PrintToFile (PrintToFile &&) = delete;
        PrintToFile& operator= (PrintToFile const&) = delete;
        PrintToFile& operator= (PrintToFile &&) = delete;

        //! Set the output precision to \p p digits and return \c *this.
        PrintToFile& SetPrecision(int p) {
            ss.precision(p);
            return *this;
        }

        //! Append \p x to the output buffer and return \c *this.
        template <typename T>
        PrintToFile& operator<< (const T& x) {
            ss << x;
            return *this;
        }

        //! Support stream manipulators such as \c std::endl and \c std::flush.
        PrintToFile& operator<< (        std::basic_ostream<char, std::char_traits<char> >&
                                 (*func)(std::basic_ostream<char, std::char_traits<char> >&))
        {
            ss << func;
            return *this;
        }

    private:

        void Initialize() {
            int my_proc = ParallelContext::MyProcSub();
            if (rank == AllProcs || rank == my_proc) {
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
        int rank;
        MPI_Comm comm;
        std::ofstream ofs;
        std::ostringstream ss;
    };

    /**
     * \class amrex::AllPrintToFile
     * \ingroup amrex_utilities
     * \brief Print to a per-process file from all MPI ranks.
     */
    class AllPrintToFile
        : public PrintToFile
    {
    public:
        //! Print on all processors of the default communicator.
        explicit AllPrintToFile (std::string file_name_)
            : PrintToFile(std::move(file_name_), Print::AllProcs)
            {}

    };
}

#endif
