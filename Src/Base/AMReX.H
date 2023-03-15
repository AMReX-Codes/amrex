
#ifndef BL_AMREX_H
#define BL_AMREX_H
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_GpuAssert.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_Exception.H>
#include <AMReX_Extension.H>

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

    typedef void (*PTR_TO_VOID_FUNC)();
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
        extern AMREX_EXPORT int handle_sigsegv;
        extern AMREX_EXPORT int handle_sigterm;
        extern AMREX_EXPORT int handle_sigint;
        extern AMREX_EXPORT int handle_sigabrt;
        extern AMREX_EXPORT int handle_sigfpe;

        extern AMREX_EXPORT int call_addr2line;
        extern AMREX_EXPORT int throw_exception;

        extern AMREX_EXPORT int regtest_reduction;

        extern AMREX_EXPORT std::ostream* osout;
        extern AMREX_EXPORT std::ostream* oserr;

        extern AMREX_EXPORT ErrorHandler error_handler;
        extern AMREX_EXPORT int abort_on_unused_inputs;
    }

    /** the AMReX "git describe" version */
    [[nodiscard]] std::string Version ();

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
                       ErrorHandler a_errhandler = nullptr);

    /**
       \brief Returns true if there are any currently-active and initialized
       AMReX instances (i.e. one for which amrex::Initialize has been called,
       and amrex::Finalize has not). Otherwise false.
     */
    [[nodiscard]] bool Initialized ();

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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Error (const char* msg = nullptr) {
#if AMREX_DEVICE_COMPILE
#if defined(NDEBUG)
        amrex::ignore_unused(msg);
#else
        if (msg) AMREX_DEVICE_PRINTF("Error %s\n", msg);
        AMREX_DEVICE_ASSERT(0);
#endif
#else
        Error_host(msg);
#endif
    }

    //! Print out warning message to cerr.
    void Warning (const std::string& msg);

    void Warning_host (const char * msg);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Warning (const char * msg) {
#if AMREX_DEVICE_COMPILE
#if defined(NDEBUG)
        amrex::ignore_unused(msg);
#else
        if (msg) AMREX_DEVICE_PRINTF("Warning %s\n", msg);
#endif
#else
        Warning_host(msg);
#endif
    }

    //! Print out message to cerr and exit via abort().
    void Abort (const std::string& msg);

    void Abort_host (const char * msg);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Abort (const char * msg = nullptr) {
#if AMREX_DEVICE_COMPILE
#if defined(NDEBUG)
        amrex::ignore_unused(msg);
#else
        if (msg) AMREX_DEVICE_PRINTF("Abort %s\n", msg);
        AMREX_DEVICE_ASSERT(0);
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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Assert (const char* EX, const char* file, int line, const char* msg = nullptr) {
#if AMREX_DEVICE_COMPILE
#if defined(NDEBUG)
        amrex::ignore_unused(EX,file,line,msg);
#else
        if (msg) {
            AMREX_DEVICE_PRINTF("Assertion `%s' failed, file \"%s\", line %d, Msg: %s",
                                EX, file, line, msg);
        } else {
            AMREX_DEVICE_PRINTF("Assertion `%s' failed, file \"%s\", line %d",
                                EX, file, line);
        }
        AMREX_DEVICE_ASSERT(0);
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

    void SetErrorHandler (ErrorHandler f);

    std::ostream& OutStream ();
    std::ostream& ErrorStream ();

    [[nodiscard]] int Verbose () noexcept;
    void SetVerbose (int v) noexcept;

    // ! Get the entire command line including the executable
    [[nodiscard]] std::string get_command ();

    // ! Get number of command line arguments after the executable
    [[nodiscard]] int command_argument_count ();

    /**
    * \brief Get command line arguments. The executable name is the
    * zero-th argument.  Return empty string if there are not that
    * many arguments.  std::string
    */
    [[nodiscard]] std::string get_command_argument (int number);

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

        amrex::ignore_unused(a_b);
        amrex::ignore_unused(a_c);
        amrex::ignore_unused(a_i);
        amrex::ignore_unused(a_l);
        amrex::ignore_unused(a_ll);
        amrex::ignore_unused(a_uc);
        amrex::ignore_unused(a_ui);
        amrex::ignore_unused(a_ul);
        amrex::ignore_unused(a_ull);
        amrex::ignore_unused(a_ll);
        amrex::ignore_unused(a_f);
        amrex::ignore_unused(a_d);
        amrex::ignore_unused(a_s);
    }
#endif

    class Geometry;

    class AMReX
    {
    public:
        AMReX ();
        ~AMReX ();
        AMReX (AMReX const&) = delete;
        AMReX (AMReX &&) = delete;
        AMReX& operator= (AMReX const&) = delete;
        AMReX& operator= (AMReX &&) = delete;

        static bool empty () noexcept { return m_instance.empty(); }

        static int size () noexcept { return static_cast<int>(m_instance.size()); }

        static AMReX* top () noexcept { return m_instance.back().get(); }

        // Thisfunction will take the ownership of the AMReX pointer,
        // and put it on the top of the stack (i.e., back of the
        // vector).  If the pointer is already in the stack, it will
        // be moved to the top.
        static void push (AMReX* pamrex);

        // This erases `pamrex` from the stack.
        static void erase (AMReX* pamrex);

        [[nodiscard]] Geometry* getDefaultGeometry () noexcept { return m_geom; }

    private:

        static AMREX_EXPORT std::vector<std::unique_ptr<AMReX> > m_instance;

        Geometry* m_geom = nullptr;
    };
}

#endif /*BL_AMREX_H*/
