
#ifndef BL_AMREX_H
#define BL_AMREX_H
#include <AMReX_Config.H>

/**
 * \file AMReX.H
 * \brief Runtime initialization/finalization helpers and global diagnostics.
 */

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_GpuAssert.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_Exception.H>
#include <AMReX_Extension.H>

#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
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

    using PTR_TO_VOID_FUNC = void (*)();
    using ErrorHandler = void (*)(const char*);

    /** \brief Compile- and run-time configuration flags; defaults are set during \c Initialize() and may be updated later (e.g., via ParmParse or setters). */
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

        extern AMREX_EXPORT bool signal_handling;
        extern AMREX_EXPORT bool handle_sigsegv;
        extern AMREX_EXPORT bool handle_sigterm;
        extern AMREX_EXPORT bool handle_sigint;
        extern AMREX_EXPORT bool handle_sigabrt;
        extern AMREX_EXPORT bool handle_sigfpe;
        extern AMREX_EXPORT bool handle_sigill;

        extern AMREX_EXPORT bool call_addr2line;
        extern AMREX_EXPORT bool throw_exception;

        extern AMREX_EXPORT bool regtest_reduction;

        extern AMREX_EXPORT std::ostream* osout;
        extern AMREX_EXPORT std::ostream* oserr;

        extern AMREX_EXPORT ErrorHandler error_handler;
        extern AMREX_EXPORT bool abort_on_unused_inputs;

        extern AMREX_EXPORT bool init_snan;
    }

    /**
     * \brief Return the AMReX version string reported by `git describe`.
     *
     * Useful for embedding build provenance into logs or metadata.
     */
    [[nodiscard]] std::string Version ();

    /**
     * \brief Initialize AMReX using an existing MPI communicator.
     *
     * This wires up stream redirection, records an optional custom error handler,
     * initializes GPU resources, and pushes a new AMReX instance onto the global stack.
     *
     * If MPI has not yet been initialized, AMReX initializes it and uses
     * MPI_COMM_WORLD.  If MPI is already initialized, AMReX duplicates
     * \p mpi_comm and uses the copy.
     *
     * The returned pointer is non-owning; pass it to Finalize(AMReX*) when done.
     * If only one AMReX instance is active, the no-argument Finalize() overload
     * can also be used.
     *
     * \param mpi_comm    Communicator used for AMReX (ignored if AMReX initializes MPI).
     * \param a_osout     Stream receiving informational output (defaults to std::cout).
     * \param a_oserr     Stream receiving diagnostic output (defaults to std::cerr).
     * \param a_errhandler Optional custom error handler (nullptr = default handler).
     * \param a_device_id GPU device index to select (-1 chooses the default policy).
     * \return Non-owning pointer to the initialized AMReX instance.
     */
    AMReX* Initialize (MPI_Comm mpi_comm,
                       std::ostream& a_osout = std::cout,
                       std::ostream& a_oserr = std::cerr,
                       ErrorHandler a_errhandler = nullptr,
                       int a_device_id = -1);

    /**
     * \brief Initialize AMReX from argc/argv with custom ParmParse handling.
     *
     * Some AMReX runtime parameters are consumed during \c Initialize().
     * These can be overridden via an inputs file or command-line arguments,
     * but if you prefer to set them programmatically, you should use the
     * \p func_parm_parse callback.  The callback is invoked after ParmParse
     * is constructed but before those parameters are read, giving the caller
     * a chance to override defaults via \c ParmParse::add.  Calling
     * \c ParmParse::add after \c Initialize() returns is too late for
     * parameters that have already been consumed.  There is no harm in
     * setting parameters that are consumed later the same way.
     *
     * If MPI has not yet been initialized, AMReX initializes it and uses
     * MPI_COMM_WORLD.  If MPI is already initialized, AMReX duplicates
     * MPI_COMM_WORLD and uses the copy.
     *
     * The returned pointer is non-owning; pass it to Finalize(AMReX*) when done.
     * If only one AMReX instance is active, the no-argument Finalize() overload
     * can also be used.
     *
     * \param argc Reference to argc (may be modified as AMReX consumes arguments).
     * \param argv Reference to argv.
     * \param func_parm_parse Callback invoked after AMReX constructs ParmParse.
     * \param a_osout Stream receiving informational output.
     * \param a_oserr Stream receiving diagnostic output.
     * \param a_errhandler Optional custom error handler.
     * \param a_device_id GPU device index to select (-1 chooses the default policy).
     * \return Non-owning pointer to the initialized AMReX instance.
     */
    AMReX* Initialize (int& argc, char**& argv,
                       const std::function<void()>& func_parm_parse,
                       std::ostream& a_osout = std::cout,
                       std::ostream& a_oserr = std::cerr,
                       ErrorHandler a_errhandler = nullptr,
                       int a_device_id = -1);

    /**
     * \brief Initialize AMReX from argc/argv, optionally constructing ParmParse.
     *
     * Some AMReX runtime parameters are consumed during \c Initialize().
     * These can be overridden via an inputs file or command-line arguments,
     * but if you prefer to set them programmatically, you should use the
     * \p func_parm_parse callback.  The callback is invoked after ParmParse
     * is constructed but before those parameters are read.  Calling
     * \c ParmParse::add after \c Initialize() returns is too late for
     * parameters that have already been consumed.  There is no harm in
     * setting parameters that are consumed later the same way.
     *
     * If MPI has not yet been initialized, AMReX initializes it and uses
     * MPI_COMM_WORLD.  If MPI is already initialized, AMReX duplicates
     * \p mpi_comm and uses the copy.
     *
     * The returned pointer is non-owning; pass it to Finalize(AMReX*) when done.
     * If only one AMReX instance is active, the no-argument Finalize() overload
     * can also be used.
     *
     * \param argc Reference to argc (may be modified as arguments are consumed).
     * \param argv Reference to argv.
     * \param build_parm_parse Whether \c argc/\c argv are parsed into ParmParse
     *                       (an empty ParmParse is still created when false).
     *                       Note that \c argv may include an inputs file
     *                       specified as \c argv[1].
     * \param mpi_comm Communicator used for AMReX (ignored if AMReX initializes MPI).
     * \param func_parm_parse Optional callback invoked after ParmParse is constructed
     *                      but before early-consumed runtime parameters are read.
     * \param a_osout Stream receiving informational output.
     * \param a_oserr Stream receiving diagnostic output.
     * \param a_errhandler Optional custom error handler.
     * \param a_device_id GPU device index to select (-1 chooses the default policy).
     * \return Non-owning pointer to the initialized AMReX instance.
     */
    AMReX* Initialize (int& argc, char**& argv, bool build_parm_parse=true,
                       MPI_Comm mpi_comm = MPI_COMM_WORLD,
                       const std::function<void()>& func_parm_parse = {},
                       std::ostream& a_osout = std::cout,
                       std::ostream& a_oserr = std::cerr,
                       ErrorHandler a_errhandler = nullptr,
                       int a_device_id = -1);

    /**
     * \brief Minimal initialization for utilities that need only part of AMReX.
     *
     * Skips signal handlers, RNG initialization, async I/O, and the default
     * multi-GPU device-mapping policy. Arena setup still runs; minimal mode
     * disables the pinned-arena warm-up and avoids the large default GPU arena
     * warm-up. FFT and other components still initialize normally.
     * \c Init_minimal calls \c Initialize(); MPI is started via
     * \c StartParallel if it is not already up.  On multi-GPU systems the caller
     * may need to select devices before calling if the default minimal policy is
     * insufficient.
     *
     * \param mpi_comm Communicator.
     */
    void Init_minimal (MPI_Comm mpi_comm = MPI_COMM_WORLD);

    /**
       \brief Returns true if there are any currently-active and initialized
       AMReX instances (i.e. one for which amrex::Initialize has been called,
       and amrex::Finalize has not). Otherwise false.
     */
    [[nodiscard]] bool Initialized ();

    /**
     * \brief Finalize AMReX state associated with the supplied instance pointer.
     *
     * Synchronizes outstanding GPU work, removes \p pamrex from the global stack,
     * runs registered \c ExecOnFinalize callbacks, tears down arenas and other
     * subsystems, restores signal handlers, and shuts down MPI when appropriate.
     *
     * \param pamrex Non-owning pointer previously returned by Initialize().
     */
    void Finalize (AMReX* pamrex);
    //! Convenience overload that finalizes the instance at the top of the stack.
    void Finalize ();
    //! Finalize a run that was initialized via Init_minimal().
    void Finalize_minimal();

    /**
     * \brief Register a callback invoked during \c Finalize() (LIFO order).
     *
     * Finalize callbacks execute in LIFO order so subsystems can release
     * global resources cleanly.  They run after the \c AMReX instance is
     * removed from the stack but before arena teardown and MPI shutdown.
     *
     * \param func Callback registered via \c ExecOnFinalize().
     */
    void ExecOnFinalize (std::function<void()> func);
    /**
     * \brief Register a callback invoked early in \c Initialize().
     *
     * Callbacks run after \c ParallelDescriptor::StartParallel() but before
     * ParmParse setup, GPU initialization, and arena setup.
     *
     * \param func Callback registered via \c ExecOnInitialize().
     */
    void ExecOnInitialize (std::function<void()> func);

    /**
     * \brief No-op helper that marks variables as intentionally unused.
     * \tparam Ts Argument types to ignore.
     */
    template <class... Ts>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void ignore_unused (const Ts&...) {}

    /**
     * \brief Print a message to stderr and abort the program.
     * \param msg Description explaining why the run is terminating.
     */
    void Error (const std::string& msg);

    /**
     * \brief Host-side error sink invoked by GPU-friendly overloads.
     *
     * The helper either dispatches to the registered error handler, throws
     * RuntimeError (when `system::throw_exception` is true), or prints the
     * message and aborts. Device wrappers call this once they return to the host.
     *
     * \param type Short label describing the failure class (e.g., "Error").
     * \param msg  Null-terminated string explaining the failure (may be nullptr).
     */
    void Error_host (const char* type, const char* msg);

    //! \copydoc amrex::Error(const std::string&)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Error (const char* msg = nullptr) {
#if defined(NDEBUG)
        AMREX_IF_ON_DEVICE((amrex::ignore_unused(msg);))
#else
        AMREX_IF_ON_DEVICE((
                if (msg) { AMREX_DEVICE_PRINTF("Error %s\n", msg); }
                AMREX_DEVICE_ASSERT(0);
        ))
#endif
        AMREX_IF_ON_HOST((Error_host("Error", msg);))
    }

    /**
     * \brief Print a warning message to the diagnostic stream and keep running.
     * \param msg Human-readable explanation of the non-fatal condition.
     */
    void Warning (const std::string& msg);

    /**
     * \brief Print a warning message to the diagnostic stream (host only).
     * \param msg Message (may be nullptr).
     */
    void Warning_host (const char * msg);

    //! \copydoc amrex::Warning(const std::string&)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Warning (const char * msg) {
#if defined(NDEBUG)
        AMREX_IF_ON_DEVICE((amrex::ignore_unused(msg);))
#else
        AMREX_IF_ON_DEVICE((if (msg) { AMREX_DEVICE_PRINTF("Warning %s\n", msg); }))
#endif
        AMREX_IF_ON_HOST((Warning_host(msg);))
    }

    /**
     * \brief Print a fatal-error message to stderr and abort execution.
     * \param msg Explanation emitted before terminating the program.
     */
    void Abort (const std::string& msg);

    //! \copydoc amrex::Abort(const std::string&)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Abort (const char * msg = nullptr) {
#if defined(NDEBUG)
        AMREX_IF_ON_DEVICE((amrex::ignore_unused(msg);))
#else
        AMREX_IF_ON_DEVICE((
                if (msg) { AMREX_DEVICE_PRINTF("Abort %s\n", msg); }
                AMREX_DEVICE_ASSERT(0);
        ))
#endif
        AMREX_IF_ON_HOST((Error_host("Abort", msg);))
    }

    /**
     * \brief Host-side assertion helper used by BL_ASSERT().
     *
     * Writes a formatted assertion failure to stderr, adds any optional
     * message, and aborts the program.
     *
     * \param EX   Stringized expression that failed.
     * \param file Source file containing the assertion.
     * \param line Line number where the assertion triggered.
     * \param msg  Optional supplemental text (may be nullptr).
     * \param msg_size Length of \p msg when supplied.
     */
    void Assert_host (const char* EX, const char* file, int line, const char* msg,
                      std::size_t msg_size=0);

    /**
     * \brief Device/host assertion helper used by BL_ASSERT.
     *
     * On the host, always forwards to \c Assert_host.  On the device, prints
     * and aborts when assertions are enabled (debug builds or when
     * \c AMREX_USE_ASSERTION is defined); in release builds with
     * \c AMREX_USE_ASSERTION disabled, device-side arguments are ignored.
     *
     * \param EX   Stringized expression that failed.
     * \param file Source file containing the assertion.
     * \param line Line number where the assertion triggered.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Assert (const char* EX, const char* file, int line) {
#if defined(NDEBUG) && !defined(AMREX_USE_ASSERTION)
        AMREX_IF_ON_DEVICE((amrex::ignore_unused(EX,file,line);))
#else
        AMREX_IF_ON_DEVICE((
            AMREX_DEVICE_PRINTF("Assertion `%s' failed, file \"%s\", line %d",
                                    EX, file, line);
            AMREX_DEVICE_ASSERT(0);
        ))
#endif
        AMREX_IF_ON_HOST((Assert_host(EX,file,line,nullptr,0U);))
    }

    /**
     * \brief Assert helper that includes an additional message in the output.
     *
     * The same host/device dispatch described for the primary overload applies:
     * on the host, forwards to \c Assert_host; on the device, prints the
     * assertion and message and aborts when assertions are enabled.
     *
     * \param EX   Stringized expression that failed.
     * \param file Source file containing the assertion.
     * \param line Line number where the assertion triggered.
     * \param msg  Non-null text providing additional context.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Assert (const char* EX, const char* file, int line, const char* msg) {
#if defined(NDEBUG) && !defined(AMREX_USE_ASSERTION)
        AMREX_IF_ON_DEVICE((amrex::ignore_unused(EX,file,line,msg);))
#else
        AMREX_IF_ON_DEVICE((
            AMREX_DEVICE_PRINTF("Assertion `%s' failed, file \"%s\", line %d, Msg: %s",
                                    EX, file, line, msg);
            AMREX_DEVICE_ASSERT(0);
        ))
#endif
        AMREX_IF_ON_HOST((Assert_host(EX,file,line,msg,std::strlen(msg));))
    }

    //! \overload amrex::Assert(const char*, const char*, int, const char*)
    AMREX_FORCE_INLINE
    void Assert (const char* EX, const char* file, int line, const std::string& msg) {
        Assert_host(EX,file,line,msg.c_str(),msg.size());
    }

    /**
     * \brief Write directly to stderr after flushing standard I/O buffers.
     *
     * Used internally by amrex::Error(), amrex::Abort(), and amrex::Assert()
     * to emit the final abort message through C stdio.
     *
     * \param str Null-terminated string describing the failure (may be nullptr).
     */
    void write_to_stderr_without_buffering (const char* str);

    /**
     * \brief Install a custom error handler invoked by amrex::Error_host.
     * \param f Callback receiving the error message; \c nullptr selects the default path.
     */
    void SetErrorHandler (ErrorHandler f);

    //! Stream that receives informational output (defaults to std::cout).
    std::ostream& OutStream ();
    //! Stream that receives diagnostic output (defaults to std::cerr).
    std::ostream& ErrorStream ();

    //! Return the verbosity level configured via ParmParse or SetVerbose().
    [[nodiscard]] int Verbose () noexcept;
    /**
     * \brief Override the verbosity level used by amrex::Print().
     * \param v Non-negative verbosity level (0 = quiet).
     */
    void SetVerbose (int v) noexcept;

    //! Whether new allocations are filled with signaling NaNs when enabled.
    [[nodiscard]] bool InitSNaN () noexcept;
    /**
     * \brief Toggle whether new allocations are initialized with signaling NaNs.
     * \param v True to fill new allocations with SNaN values, false to skip.
     */
    void SetInitSNaN (bool v) noexcept;

    //! Get the entire command line including the executable path.
    [[nodiscard]] std::string get_command ();

    //! Get the number of command line arguments excluding argv[0].
    [[nodiscard]] int command_argument_count ();

    /**
     * \brief Fetch a copy of argv[number] (argv[0] is the executable path).
     * \param number Argument index to query.
     * \return Requested argument or an empty string when \p number is out of range.
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

    /**
     * \brief Runtime context stored for each active amrex::Initialize() call.
     *
     * Instances live on an internal stack so nested initialization sequences
     * can cooperate.
     */
    class AMReX
    {
    public:
        AMReX ();
        ~AMReX ();
        AMReX (AMReX const&) = delete;
        AMReX (AMReX &&) = delete;
        AMReX& operator= (AMReX const&) = delete;
        AMReX& operator= (AMReX &&) = delete;

        //! True if no AMReX instances are registered on the stack.
        static bool empty () noexcept { return m_instance.empty(); }

        //! Number of managed instances currently tracked (Initialize depth).
        static int size () noexcept { return static_cast<int>(m_instance.size()); }

        //! Pointer to the most recently pushed instance.
        static AMReX* top () noexcept { return m_instance.back().get(); }

        /**
         * \brief Transfer ownership of an existing AMReX* into the stack.
         *
         * If the pointer already lives in the stack, it is rotated to the top.
         *
         * \param pamrex Raw pointer previously allocated with `new AMReX`.
         */
        static void push (AMReX* pamrex);

        //! Push a freshly created AMReX instance stored in a unique_ptr.
        static void push (std::unique_ptr<AMReX> pamrex);

        //! Remove a specific instance pointer from the internal stack.
        static void erase (AMReX* pamrex);

        //! Geometry owned by this context (constructed with the \c AMReX object).
        [[nodiscard]] Geometry* getDefaultGeometry () noexcept { return m_geom; }

    private:

        static AMREX_EXPORT std::vector<std::unique_ptr<AMReX> > m_instance;

        Geometry* m_geom = nullptr;
    };

    enum struct FPExcept : std::uint8_t {
        none     = 0B0000,
        invalid  = 0B0001,
        zero     = 0B0010,
        overflow = 0B0100,
        all      = 0B0111
    };

    //! Return true when \p a is not \c FPExcept::none (i.e., the mask is non-empty).
    [[nodiscard]] inline bool any (FPExcept a) { return a != FPExcept::none; }

    //! Bitwise union of two floating-point exception masks.
    [[nodiscard]] inline FPExcept operator| (FPExcept a, FPExcept b)
    {
        using T = std::underlying_type_t<FPExcept>;
        // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
        return static_cast<FPExcept>(static_cast<T>(a) | static_cast<T>(b));
    }

    //! Bitwise intersection of two floating-point exception masks.
    [[nodiscard]] inline FPExcept operator& (FPExcept a, FPExcept b)
    {
        using T = std::underlying_type_t<FPExcept>;
        return static_cast<FPExcept>(static_cast<T>(a) & static_cast<T>(b));
    }

    //! Return currently enabled FP exceptions. Linux only.
    [[nodiscard]] FPExcept getFPExcept ();

    //! Set FP exception traps. Linux only. This enables set flags and
    //! DISABLES unset flags. This can be used to restore previous settings.
    //! \param excepts Bit mask describing which traps should be enabled.
    FPExcept setFPExcept (FPExcept excepts);

    /**
     * \brief Disable FP exceptions. Linux Only
     *
     * This function disables given exception traps and keeps the status of
     * the others. The example below disables FPE invalid and
     * divide-by-zero, and later restores the previous settings.
     *
     \code{.cpp}
         auto prev_excepts = disableFPExcept(FPExcept::invalid | FPExcept::zero);
         // ....
         setFPExcept(prev_excepts); // restore previous settings
     \endcode
     * \param excepts Bit mask describing which traps should be disabled.
     */
    [[nodiscard]] FPExcept disableFPExcept (FPExcept excepts);

    /**
     * \brief Enable FP exceptions. Linux Only
     *
     * This function enables given exception traps and keeps the status of
     * the others. The example below enables all FPE traps, and later
     * restores the previous settings.
     *
     \code{.cpp}
         auto prev_excepts = enableFPExcept(FPExcept::all);
         // ....
         setFPExcept(prev_excepts); // restore previous settings
     \endcode
     * \param excepts Bit mask describing which traps should be enabled.
     */
    [[nodiscard]] FPExcept enableFPExcept (FPExcept excepts);
}

#endif /*BL_AMREX_H*/
