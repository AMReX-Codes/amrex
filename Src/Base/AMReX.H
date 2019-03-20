
#ifndef BL_AMREX_H
#define BL_AMREX_H

#include <AMReX_GpuQualifiers.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_Exception.H>

#include <iostream>
#include <functional>
#include <string>
#include <cstdio>

//
// Initialize, Finalize, Error Reporting, and Version String Functions

/*
 This class consists of initialize, finalize, error-reporting, and version
 functions that are used throughout AMReX.  Note that all the output functions
 add a terminating exclamation mark, so there is no need to add any punctuation
 to sentences that get output explicitly.
*/

namespace amrex
{

    typedef void (*PTR_TO_VOID_FUNC)(void);
    typedef void (*ErrorHandler) (const char* msg);

    namespace system
    {
#ifndef AMREX_DEBUG
	constexpr bool NDebug = true;
        constexpr bool  Debug = false;
#else
	constexpr bool NDebug = false;
        constexpr bool  Debug = true;
#endif

	extern std::string exename;
	
	extern int verbose;

        extern int signal_handling;
        extern int call_addr2line;
        extern int throw_exception;

        extern int regtest_reduction;

        extern std::ostream* osout;
        extern std::ostream* oserr;
        
        extern ErrorHandler error_handler;
    }

    std::string Version ();

    void Initialize (MPI_Comm mpi_comm,
                     std::ostream& a_osout = std::cout,
                     std::ostream& a_oserr = std::cerr,
                     ErrorHandler a_errhandler = nullptr);

    void Initialize (int& argc, char**& argv, bool build_parm_parse=true,
                     MPI_Comm mpi_comm = MPI_COMM_WORLD,
                     const std::function<void()>& func_parm_parse = {},
                     std::ostream& a_osout = std::cout,
                     std::ostream& a_oserr = std::cerr,
                     ErrorHandler a_errorhandler = nullptr);

    void Finalize(bool finalize_parallel = true);
    /**
    * \brief We maintain a stack of functions that need to be called in Finalize().
    * The functions are called in LIFO order.  The idea here is to allow
    * classes to clean up any "global" state that they maintain when we're
    * exiting from AMReX.
    */
    void ExecOnFinalize (PTR_TO_VOID_FUNC);
    void ExecOnInitialize (PTR_TO_VOID_FUNC);

    //! Print out message to cerr and exit via amrex::Abort().
    AMREX_GPU_HOST_DEVICE
    void Error (const char * msg = 0);
    void Error (const std::string& msg);

    //! Print out warning message to cerr.
    AMREX_GPU_HOST_DEVICE
    void Warning (const char * msg);
    void Warning (const std::string& msg);

    //! Print out message to cerr and exit via abort().
    AMREX_GPU_HOST_DEVICE
    void Abort (const char * msg = 0);
    void Abort (const std::string& msg);

    /**
    * \brief Prints assertion failed messages to cerr and exits
    * via abort().  Intended for use by the BL_ASSERT() macro
    * in <AMReX_BLassert.H>.
    */
    AMREX_GPU_HOST_DEVICE
    void Assert (const char* EX,
                 const char* file,
                 int         line,
                 const char* msg = nullptr);
    /**
    * \brief Prints out an out-of-memory message and aborts.  It is
    *         called by various routines when a call to new fails.
    *
    *        Called as amrex::OutOfMemory(__FILE__, __LINE__);
    */
    void OutOfMemory (const char* file,
                      int         line);

    /**
    * \brief This is used by amrex::Error(), amrex::Abort(), and amrex::Assert()
    * to ensure that when writing the message to stderr, that no additional
    * heap-based memory is allocated.
    */
    void write_to_stderr_without_buffering (const char* str);

    void SetErrorHandler (ErrorHandler errhandler);

    std::ostream& OutStream ();
    std::ostream& ErrorStream ();

    int Verbose () noexcept;
    void SetVerbose (int v) noexcept;

    // ! Get the entire command line including the executable
    std::string get_command ();

    // ! Get number of command line arguments after the executable
    int command_argument_count ();

    /**
    * \brief Get command line arguments. The executable name is the
    * zero-th argument.  Return empty string if there are not that
    * many arguments.  std::string
    */
    std::string get_command_argument (int number);

    inline void GccPlacater ()
    {
        std::allocator<bool> a_b;
        std::allocator<char> a_c;
        std::allocator<int> a_i;
        std::allocator<long> a_l;
        std::allocator<unsigned char> a_uc;
        std::allocator<unsigned int> a_ui;
        std::allocator<unsigned long> a_ul;
        std::allocator<float> a_f;
        std::allocator<double> a_d;
        std::allocator<std::string> a_s;
    }
}

#endif /*BL_AMREX_H*/
