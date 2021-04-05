
#ifndef BL_AMREX_H
#define BL_AMREX_H
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_GpuAssert.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_Exception.H>
#include <AMReX_Extension.H>
#include <AMReX_INT.H>
#include <AMReX_REAL.H>
#include <AMReX_Math.H>

#include <cstdio>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

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
    class AMReX;

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

        extern AMREX_EXPORT std::string exename;

        extern AMREX_EXPORT int verbose;

        extern AMREX_EXPORT int signal_handling;
        extern AMREX_EXPORT int call_addr2line;
        extern AMREX_EXPORT int throw_exception;

        extern AMREX_EXPORT int regtest_reduction;

        extern AMREX_EXPORT std::ostream* osout;
        extern AMREX_EXPORT std::ostream* oserr;

        extern AMREX_EXPORT ErrorHandler error_handler;
        extern AMREX_EXPORT int abort_on_unused_inputs;
    }

    std::string Version ();

    // The returned AMReX* is non-owning! To delete it, call Finalize(AMReX*).
    AMReX* Initialize (MPI_Comm mpi_comm,
                       std::ostream& a_osout = std::cout,
                       std::ostream& a_oserr = std::cerr,
                       ErrorHandler a_errhandler = nullptr);

    // The returned AMReX* is non-owning! To delete it, call Finalize(AMReX*).
    AMReX* Initialize (int& argc, char**& argv, bool build_parm_parse=true,
                       MPI_Comm mpi_comm = MPI_COMM_WORLD,
                       const std::function<void()>& func_parm_parse = {},
                       std::ostream& a_osout = std::cout,
                       std::ostream& a_oserr = std::cerr,
                       ErrorHandler a_errorhandler = nullptr);

    void Finalize (AMReX* pamrex);
    void Finalize (); // Finalize the current top
    /**
    * \brief We maintain a stack of functions that need to be called in Finalize().
    * The functions are called in LIFO order.  The idea here is to allow
    * classes to clean up any "global" state that they maintain when we're
    * exiting from AMReX.
    */
    void ExecOnFinalize (PTR_TO_VOID_FUNC);
    void ExecOnInitialize (PTR_TO_VOID_FUNC);

    //! This shuts up the compiler about unused variables
    template <class... Ts>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void ignore_unused (const Ts&...) {}

    //! Print out message to cerr and exit via amrex::Abort().
    void Error (const std::string& msg);

    void Error_host (const char* msg);

#if defined(AMREX_USE_GPU)
    AMREX_GPU_EXTERNAL AMREX_GPU_DEVICE AMREX_NO_INLINE
    void Error_device (const char * msg);
#endif

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Error (const char* msg = 0) {
#if AMREX_DEVICE_COMPILE
#ifdef NDEBUG
        amrex::ignore_unused(msg);
#else
        Error_device(msg);
#endif
#else
        Error_host(msg);
#endif
    }

    //! Print out warning message to cerr.
    void Warning (const std::string& msg);

    void Warning_host (const char * msg);

#if defined(AMREX_USE_GPU)
    AMREX_GPU_EXTERNAL AMREX_GPU_DEVICE AMREX_NO_INLINE
    void Warning_device (const char * msg);
#endif

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Warning (const char * msg) {
#if AMREX_DEVICE_COMPILE
#ifdef NDEBUG
        amrex::ignore_unused(msg);
#else
        Warning_device(msg);
#endif
#else
        Warning_host(msg);
#endif
    }

    //! Print out message to cerr and exit via abort().
    void Abort (const std::string& msg);

    void Abort_host (const char * msg);

#if defined(AMREX_USE_GPU)
    AMREX_GPU_EXTERNAL AMREX_GPU_DEVICE AMREX_NO_INLINE
    void Abort_device (const char * msg);
#endif

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Abort (const char * msg = 0) {
#if AMREX_DEVICE_COMPILE
#ifdef NDEBUG
        amrex::ignore_unused(msg);
#else
        Abort_device(msg);
#endif
#else
        Abort_host(msg);
#endif
    }

    /**
    * \brief Prints assertion failed messages to cerr and exits
    * via abort().  Intended for use by the BL_ASSERT() macro
    * in <AMReX_BLassert.H>.
    */

    void Assert_host (const char* EX, const char* file, int line, const char* msg);

#if defined(AMREX_USE_GPU)
    AMREX_GPU_EXTERNAL AMREX_GPU_DEVICE AMREX_NO_INLINE
    void Assert_device (const char* EX, const char* file, int line, const char* msg);
#endif

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Assert (const char* EX, const char* file, int line, const char* msg = nullptr) {
#if AMREX_DEVICE_COMPILE
#ifdef NDEBUG
        amrex::ignore_unused(EX,file,line,msg);
#else
        Assert_device(EX,file,line,msg);
#endif
#else
        Assert_host(EX,file,line,msg);
#endif
    }

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

#ifndef _MSC_VER
    inline void GccPlacater ()
    {
        std::allocator<bool> a_b;
        std::allocator<char> a_c;
        std::allocator<int> a_i;
        std::allocator<long> a_l;
        std::allocator<long long> a_ll;
        std::allocator<unsigned char> a_uc;
        std::allocator<unsigned int> a_ui;
        std::allocator<unsigned long> a_ul;
        std::allocator<unsigned long long> a_ull;
        std::allocator<float> a_f;
        std::allocator<double> a_d;
        std::allocator<std::string> a_s;
    }
#endif

    class Geometry;

    class AMReX
    {
    public:
        AMReX ();
        ~AMReX ();

        static bool empty () noexcept { return m_instance.empty(); }

        static int size () noexcept { return m_instance.size(); }

        static AMReX* top () noexcept { return m_instance.back().get(); }

        // Thisfunction will take the ownership of the AMReX pointer,
        // and put it on the top of the stack (i.e., back of the
        // vector).  If the pointer is already in the stack, it will
        // be moved to the top.
        static void push (AMReX* pamrex);

        // This erases `pamrex` from the stack.
        static void erase (AMReX* pamrex);

        Geometry* getDefaultGeometry () noexcept { return m_geom; }

    private:

        static std::vector<std::unique_ptr<AMReX> > m_instance;

        Geometry* m_geom = nullptr;
    };
}

#endif /*BL_AMREX_H*/
