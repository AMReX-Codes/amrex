#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Box.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BLFort.H>
#include <AMReX_Utility.H>
#include <AMReX_Random.H>
#include <AMReX_Print.H>
#include <AMReX_Arena.H>
#include <AMReX_BLBackTrace.H>
#include <AMReX_MemPool.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>

#ifdef AMREX_USE_FFT
#include <AMReX_FFT.H>
#endif

#ifdef AMREX_USE_HYPRE
#include <_hypre_utilities.h>
#ifdef AMREX_USE_CUDA
#include <_hypre_utilities.hpp>
#endif
#endif

#ifdef AMREX_USE_SUNDIALS
#include <AMReX_Sundials.H>
#endif

#ifdef AMREX_USE_CUPTI
#include <AMReX_CuptiTrace.H>
#endif

#include <AMReX_Machine.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#endif

#ifndef BL_AMRPROF
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_AsyncOut.H>
#endif

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#ifdef AMREX_USE_OMP
#include <AMReX_OpenMP.H>
#include <omp.h>
#endif

#if defined(__APPLE__) && defined(__x86_64__)
#include <xmmintrin.h>
#endif

#if !(defined(_MSC_VER) && defined(__CUDACC__))
//MSVC can't pre-processor cfenv with `Zc:preprocessor`
//https://developercommunity.visualstudio.com/content/problem/1271183/zcpreprocessor-e-crashes-when-given.html
#include <cfenv>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <csignal>
#include <iostream>
#include <iomanip>
#include <new>
#include <optional>
#include <stack>
#include <string>
#include <thread>
#include <limits>
#include <vector>
#include <algorithm>

#ifdef AMREX_USE_COVERITY
namespace {
    // coverity[+kill]
    void amrex_coverity_abort()
    {
        std::exit(EXIT_FAILURE);
    }
}
#endif

namespace amrex {

std::vector<std::unique_ptr<AMReX> > AMReX::m_instance;

namespace system
{
    std::string exename;
    int verbose = 1;
    bool signal_handling;
    bool handle_sigsegv;
    bool handle_sigterm;
    bool handle_sigint;
    bool handle_sigabrt;
    bool handle_sigfpe;
    bool handle_sigill;
    bool call_addr2line;
    bool throw_exception;
    bool regtest_reduction;
    bool abort_on_unused_inputs = false;
    std::ostream* osout = &std::cout;
    std::ostream* oserr = &std::cerr;
    ErrorHandler error_handler = nullptr;
#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
    bool init_snan = true;
#else
    bool init_snan = false;
#endif
}
}

namespace {
    std::string command_line;
    std::vector<std::string> command_arguments;
}

namespace {
    std::streamsize  prev_out_precision;
    std::streamsize  prev_err_precision;
    std::new_handler prev_new_handler;
    using SignalHandler = void (*)(int);
    SignalHandler prev_handler_sigsegv = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
    SignalHandler prev_handler_sigterm = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
    SignalHandler prev_handler_sigint  = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
    SignalHandler prev_handler_sigabrt = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
    SignalHandler prev_handler_sigfpe  = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
    SignalHandler prev_handler_sigill  = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
#if defined(__linux__)
    int           prev_fpe_excepts = 0;
    int           curr_fpe_excepts = 0;
#elif defined(__APPLE__) && defined(__x86_64__)
    unsigned int  prev_fpe_mask = 0u;
    unsigned int  curr_fpe_excepts = 0u;
#endif
}

#ifdef AMREX_USE_HYPRE
namespace {
    bool init_hypre = true;
#if defined(HYPRE_USING_CUDA) || defined(HYPRE_USING_HIP)
    bool hypre_spgemm_use_vendor = false;
    bool hypre_spmv_use_vendor = false;
    bool hypre_sptrans_use_vendor = false;
#endif
}
#endif

int amrex::Verbose () noexcept { return amrex::system::verbose; }

void amrex::SetVerbose (int v) noexcept { amrex::system::verbose = v; }

bool amrex::InitSNaN () noexcept { return amrex::system::init_snan; }

void amrex::SetInitSNaN (bool v) noexcept  { amrex::system::init_snan = v; }

void amrex::SetErrorHandler (amrex::ErrorHandler f) {
    amrex::system::error_handler = f;
}

//
// This is used by amrex::Error(), amrex::Abort(), and amrex::Assert()
// to ensure that when writing the message to stderr, that no additional
// heap-based memory is allocated.
//

void
amrex::write_to_stderr_without_buffering (const char* str)
{
    //
    // Flush all buffers.
    //
    std::fflush(nullptr);

    if (str)
    {
        std::ostringstream procall;
        procall << ParallelDescriptor::MyProc() << "::";
        auto tmp = procall.str();
        const char *cprocall = tmp.c_str();
        const char * const end = " !!!\n";
        std::fwrite(cprocall, strlen(cprocall), 1, stderr);
        std::fwrite(str, strlen(str), 1, stderr);
        std::fwrite(end, strlen(end), 1, stderr);
    }
}

namespace {
void
write_lib_id(const char* msg)
{
    std::fflush(nullptr);
    const char* const s = "amrex::";
    std::fwrite(s, strlen(s), 1, stderr);
    if ( msg )
    {
        std::fwrite(msg, strlen(msg), 1, stderr);
        std::fwrite("::", 2, 1, stderr);
    }
}
}

void
amrex::Error (const std::string& msg)
{
    Error(msg.c_str());
}

void
amrex::Abort (const std::string& msg)
{
    Abort(msg.c_str());
}

void
amrex::Warning (const std::string& msg)
{
    Warning(msg.c_str());
}

void
amrex::Error_host (const char* type, const char * msg)
{
    amrex::ignore_unused(type);
#ifdef AMREX_USE_COVERITY
    amrex_coverity_abort();
#else
    if (system::error_handler) {
        system::error_handler(msg);
    } else if (system::throw_exception) {
        throw RuntimeError(msg);
    } else {
        write_lib_id(type);
        write_to_stderr_without_buffering(msg);
#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_abort_omp_critical)
#endif
        ParallelDescriptor::Abort();
    }
#endif
}

void
amrex::Warning_host (const char * msg)
{
    if (msg) {
        amrex::Print(Print::AllProcs,amrex::ErrorStream()) << msg << '!' << '\n';
    }
}

void
amrex::Assert_host (const char* EX, const char* file, int line, const char* msg)
{
#ifdef AMREX_USE_COVERITY
    amrex_coverity_abort();
#else
    const int N = 512;

    char buf[N];

    if (msg) {
        snprintf(buf,
                 N,
                 "Assertion `%s' failed, file \"%s\", line %d, Msg: %s",
                 EX,
                 file,
                 line,
                 msg);
    } else {
        snprintf(buf,
                 N,
                 "Assertion `%s' failed, file \"%s\", line %d",
                 EX,
                 file,
                 line);
    }

    if (system::error_handler) {
        system::error_handler(buf);
    } else if (system::throw_exception) {
        throw RuntimeError(buf);
    } else {
       write_to_stderr_without_buffering(buf);
#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_abort_omp_critical)
#endif
       ParallelDescriptor::Abort();
   }
#endif
}

namespace
{
    std::stack<std::function<void()>> The_Finalize_Function_Stack;
    std::stack<std::function<void()>> The_Initialize_Function_Stack;
}

void
amrex::ExecOnFinalize (std::function<void()> f)
{
    The_Finalize_Function_Stack.push(std::move(f));
}

void
amrex::ExecOnInitialize (std::function<void()> f)
{
    The_Initialize_Function_Stack.push(std::move(f));
}

amrex::AMReX*
amrex::Initialize (MPI_Comm mpi_comm,
                   std::ostream& a_osout, std::ostream& a_oserr,
                   ErrorHandler a_errhandler)
{
    int argc = 0;
    char** argv = nullptr;
    return Initialize(argc, argv, false, mpi_comm, {}, a_osout, a_oserr, a_errhandler);
}

amrex::AMReX*
amrex::Initialize (int& argc, char**& argv, bool build_parm_parse,
                   MPI_Comm mpi_comm, const std::function<void()>& func_parm_parse,
                   std::ostream& a_osout, std::ostream& a_oserr,
                   ErrorHandler a_errhandler)
{
    system::exename.clear();
//    system::verbose = 0;
    system::regtest_reduction = false;
    system::signal_handling = true;
    system::handle_sigsegv = true;
    system::handle_sigterm = false;
    system::handle_sigint  = true;
    system::handle_sigabrt = true;
    system::handle_sigfpe  = true;
    system::handle_sigill  = true;
    system::call_addr2line = true;
    system::throw_exception = false;
    system::osout = &a_osout;
    system::oserr = &a_oserr;
    system::error_handler = a_errhandler;

    ParallelDescriptor::StartParallel(&argc, &argv, mpi_comm);

    prev_out_precision = system::osout->precision(10);
    prev_err_precision = system::oserr->precision(10);

#ifdef AMREX_PMI
    ParallelDescriptor::PMI_Initialize();
#endif

    //
    // Make sure to catch new failures.
    //
    prev_new_handler = std::set_new_handler(amrex::OutOfMemory);

    command_line.clear();
    command_arguments.clear();

    if (argc > 0)
    {
        if (argv[0][0] != '/') {
            system::exename = FileSystem::CurrentPath();
            system::exename += "/";
        }
        system::exename += argv[0];

        for (int i = 0; i < argc; ++i) {
            if (i != 0) { command_line.append(" "); }
            command_line.append(argv[i]);
            command_arguments.emplace_back(argv[i]);
        }
    }

#if defined(AMREX_USE_UPCXX)
    upcxx::init();
#endif

    while ( ! The_Initialize_Function_Stack.empty())
    {
        //
        // Call the registered function.
        //
        The_Initialize_Function_Stack.top()();
        //
        // And then remove it from the stack.
        //
        The_Initialize_Function_Stack.pop();
    }

    BL_PROFILE_INITIALIZE();

#ifndef BL_AMRPROF
    if (build_parm_parse)
    {
        if (argc == 1)
        {
            ParmParse::Initialize(0,nullptr,nullptr);
        }
        else if (argc > 1)
        {
            if (argv[1][0] == '-')
            {
                // If arguments list starts with "-", do not use ParmParse.
                // Application code can then parse the command line. This will
                // prevent "-h" or "--help" from creating errors in ParmParse,
                // but only if it's the first argument after the executable.
                ParmParse::Initialize(0,nullptr,nullptr);
            }
            else
            {
                // This counts command line arguments before a "--"
                // and only sends the preceding arguments to ParmParse;
                // the rest get ignored.
                int ppargc = 1;
                for (; ppargc < argc; ++ppargc) {
                    if (std::strcmp(argv[ppargc], "--") == 0) { break; }
                }
                if (ppargc > 1)
                {
                    if (std::strchr(argv[1],'=') || (argc > 2 ? argv[2][0] == '=' : false) )
                    {
                        // No inputs file to parse
                        ParmParse::Initialize(ppargc-1,argv+1,nullptr);
                    }
                    else
                    {
                        // argv[1] is an inputs file
                        ParmParse::Initialize(ppargc-2,argv+2,argv[1]);
                    }
                }
            }
        }
    } else {
        ParmParse::Initialize(0,nullptr,nullptr);
    }

    if (func_parm_parse) {
        func_parm_parse();
    }

    {
        ParmParse pp("amrex");
        if (! pp.query("verbose", "v", system::verbose)) {
            pp.add("verbose", system::verbose);
        }
        pp.queryAdd("init_snan", system::init_snan);
    }

    if (system::verbose > 0) {
        amrex::Print() << "Initializing AMReX (" << amrex::Version() << ")...\n";
    }

#ifdef AMREX_USE_MPI
    if (system::verbose > 0) {
        amrex::Print() << "MPI initialized with "
                       << ParallelDescriptor::NProcs()
                       << " MPI processes\n";
        int provided = -1;
        MPI_Query_thread(&provided);
        amrex::Print() << "MPI initialized with thread support level " << provided << '\n';
    }
#endif

#ifdef AMREX_USE_OMP
    amrex::OpenMP::Initialize();

    // status output
    if (system::verbose > 0) {
//    static_assert(_OPENMP >= 201107, "OpenMP >= 3.1 is required.");
        amrex::Print() << "OMP initialized with "
                       << omp_get_max_threads()
                       << " OMP threads\n";
    }

    // warn if over-subscription is detected
    if (system::verbose > 0) {
        auto ncores = int(std::thread::hardware_concurrency());
        if (ncores != 0 && // It might be zero according to the C++ standard.
            ncores < omp_get_max_threads() * ParallelDescriptor::NProcsPerNode())
        {
            amrex::Print(amrex::ErrorStream())
                << "AMReX Warning: You might be oversubscribing CPU cores with OMP threads.\n"
                << "               There are " << ncores << " cores per node.\n"
#if defined(AMREX_USE_MPI)
                << "               There are " << ParallelDescriptor::NProcsPerNode() << " MPI ranks (processes) per node.\n"
#endif
                << "               But OMP is initialized with " << omp_get_max_threads() << " threads per process.\n"
                << "               You should consider setting OMP_NUM_THREADS="
                << ncores/ParallelDescriptor::NProcsPerNode() << " or less in the environment.\n";
        }
    }
#endif

    Machine::Initialize();

#ifdef AMREX_USE_GPU
    // Initialize after ParmParse so that we can read inputs.
    Gpu::Device::Initialize();
#ifdef AMREX_USE_CUPTI
    CuptiInitialize();
#endif
#endif

    {
        ParmParse pp("amrex");
        pp.query("regtest_reduction", system::regtest_reduction);
        pp.queryAdd("signal_handling", system::signal_handling);
        pp.queryAdd("throw_exception", system::throw_exception);
        pp.query("call_addr2line", system::call_addr2line);
        pp.queryAdd("abort_on_unused_inputs", system::abort_on_unused_inputs);

#ifdef AMREX_USE_SYCL
        // Disable SIGSEGV handling by default for Intel GPUs, because it is
        // currently used by their managed memory implementation with discrete
        // GPUs
        if (Gpu::Device::deviceVendor().find("Intel") != std::string::npos) {
            system::handle_sigsegv = 0;
        }
#endif

        if (system::signal_handling)
        {
            pp.queryAdd("handle_sigsegv", system::handle_sigsegv);
            pp.queryAdd("handle_sigterm", system::handle_sigterm);
            pp.queryAdd("handle_sigint" , system::handle_sigint );
            pp.queryAdd("handle_sigabrt", system::handle_sigabrt);
            pp.queryAdd("handle_sigfpe" , system::handle_sigfpe );
            pp.queryAdd("handle_sigill" , system::handle_sigill );

            // We save the signal handlers and restore them in Finalize.
            if (system::handle_sigsegv) {
                prev_handler_sigsegv = std::signal(SIGSEGV, BLBackTrace::handler);
            } else {
                prev_handler_sigsegv = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
            }

            if (system::handle_sigterm) {
                prev_handler_sigterm = std::signal(SIGTERM,  BLBackTrace::handler);
            } else {
                prev_handler_sigterm = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
            }

            if (system::handle_sigint) {
                prev_handler_sigint = std::signal(SIGINT,  BLBackTrace::handler);
            } else {
                prev_handler_sigint = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
            }

            if (system::handle_sigabrt) {
                prev_handler_sigabrt = std::signal(SIGABRT, BLBackTrace::handler);
            } else {
                prev_handler_sigabrt = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
            }

            if (system::handle_sigfpe) {
                prev_handler_sigfpe = std::signal(SIGFPE,  BLBackTrace::handler);
            } else {
                prev_handler_sigfpe = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
            }

            if (system::handle_sigill) {
                prev_handler_sigill = std::signal(SIGILL,  BLBackTrace::handler);
            } else {
                prev_handler_sigill = SIG_ERR; // NOLINT(performance-no-int-to-ptr)
            }

            if (system::handle_sigfpe)
            {
                bool invalid = false, divbyzero=false, overflow=false;
                pp.queryAdd("fpe_trap_invalid", invalid);
                pp.queryAdd("fpe_trap_zero", divbyzero);
                pp.queryAdd("fpe_trap_overflow", overflow);

#if defined(__linux__)
                curr_fpe_excepts = 0;
                if (invalid)   { curr_fpe_excepts |= FE_INVALID;   }
                if (divbyzero) { curr_fpe_excepts |= FE_DIVBYZERO; }
                if (overflow)  { curr_fpe_excepts |= FE_OVERFLOW;  }
                prev_fpe_excepts = fegetexcept();
                if (curr_fpe_excepts != 0) {
                    feenableexcept(curr_fpe_excepts);  // trap floating point exceptions
                }

#elif defined(__APPLE__) && defined(__x86_64__)
                prev_fpe_mask = _MM_GET_EXCEPTION_MASK();
                curr_fpe_excepts = 0u;
                if (invalid)   { curr_fpe_excepts |= _MM_MASK_INVALID;  }
                if (divbyzero) { curr_fpe_excepts |= _MM_MASK_DIV_ZERO; }
                if (overflow)  { curr_fpe_excepts |= _MM_MASK_OVERFLOW; }
                if (curr_fpe_excepts != 0u) {
                    _MM_SET_EXCEPTION_MASK(prev_fpe_mask & ~curr_fpe_excepts);
                }
#endif
            }

#if defined(__APPLE__) && defined(__aarch64__)
            if (system::handle_sigill)
            {
                bool invalid = false, divbyzero=false, overflow=false;
                pp.queryAdd("fpe_trap_invalid", invalid);
                pp.queryAdd("fpe_trap_zero", divbyzero);
                pp.queryAdd("fpe_trap_overflow", overflow);

                fenv_t env;
                fegetenv(&env);
                if (invalid)   { env.__fpcr |= __fpcr_trap_invalid;   }
                if (divbyzero) { env.__fpcr |= __fpcr_trap_divbyzero; }
                if (overflow)  { env.__fpcr |= __fpcr_trap_overflow;  }
                fesetenv(&env);
                // SIGILL ref: https://developer.apple.com/forums/thread/689159
            }
#endif
        }

#ifdef AMREX_USE_HYPRE
        pp.queryAdd("init_hypre", init_hypre);
#if defined(HYPRE_USING_CUDA) || defined(HYPRE_USING_HIP)
        pp.queryAdd("hypre_spgemm_use_vendor", hypre_spgemm_use_vendor);
        pp.queryAdd("hypre_spmv_use_vendor", hypre_spmv_use_vendor);
        pp.queryAdd("hypre_sptrans_use_vendor", hypre_sptrans_use_vendor);
#endif
#endif
    }

    ParallelDescriptor::Initialize();

    BL_TINY_PROFILE_MEMORYINITIALIZE();
    Arena::Initialize();
    amrex_mempool_init();

    //
    // Initialize random seed after we're running in parallel.
    //
    amrex::InitRandom(ParallelDescriptor::MyProc()+1, ParallelDescriptor::NProcs());

    // For thread safety, we should do these initializations here.
    BaseFab_Initialize();
    BoxArray::Initialize();
    DistributionMapping::Initialize();
    FArrayBox::Initialize();
    IArrayBox::Initialize();
    FabArrayBase::Initialize();
    MultiFab::Initialize();
    iMultiFab::Initialize();
    VisMF::Initialize();
    AsyncOut::Initialize();
    VectorGrowthStrategy::Initialize();

#ifdef AMREX_USE_FFT
    FFT::Initialize();
#endif

#ifdef AMREX_USE_EB
    EB2::Initialize();
#endif

    BL_PROFILE_INITPARAMS();
#endif // ifndef BL_AMRPROF

#ifdef AMREX_USE_HYPRE
    if (init_hypre) {
        HYPRE_Init();
#if defined(HYPRE_USING_CUDA) || defined(HYPRE_USING_HIP)

#if defined(HYPRE_RELEASE_NUMBER) && (HYPRE_RELEASE_NUMBER >= 22400)

#ifdef HYPRE_USING_DEVICE_POOL
        /* device pool allocator */
        hypre_uint mempool_bin_growth   = 8,
            mempool_min_bin      = 3,
            mempool_max_bin      = 9;
        size_t mempool_max_cached_bytes = 2000LL * 1024 * 1024;

        /* To be effective, hypre_SetCubMemPoolSize must immediately follow HYPRE_Init */
        HYPRE_SetGPUMemoryPoolSize( mempool_bin_growth, mempool_min_bin,
                                    mempool_max_bin, mempool_max_cached_bytes );
#endif
#if (HYPRE_RELEASE_NUMBER >= 22500)
        HYPRE_SetSpGemmUseVendor(hypre_spgemm_use_vendor);
        HYPRE_SetSpMVUseVendor(hypre_spmv_use_vendor);
        HYPRE_SetSpTransUseVendor(hypre_sptrans_use_vendor);
#elif (HYPRE_USING_CUDA)
        HYPRE_SetSpGemmUseCusparse(hypre_spgemm_use_vendor);
#endif
        HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
        HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
        HYPRE_SetUseGpuRand(true);
#else
        hypre_HandleDefaultExecPolicy(hypre_handle()) = HYPRE_EXEC_DEVICE;
        hypre_HandleSpgemmUseCusparse(hypre_handle()) = 0;
#endif
#endif

    }
#endif

#ifdef AMREX_USE_SUNDIALS
    sundials::Initialize(amrex::OpenMP::get_max_threads());
#endif

    if (system::verbose > 0) {
        amrex::Print() << "AMReX (" << amrex::Version() << ") initialized" << '\n';
    }

    BL_TINY_PROFILE_INITIALIZE();

    AMReX::push(new AMReX()); // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)
    return AMReX::top(); // NOLINT
}

bool
amrex::Initialized ()
{
    return !amrex::AMReX::empty();
}

void
amrex::Finalize ()
{
    amrex::Finalize(AMReX::top());
}

void
amrex::Finalize (amrex::AMReX* pamrex)
{
#ifdef AMREX_USE_GPU
    Gpu::streamSynchronizeAll();
#endif

    AMReX::erase(pamrex);

#ifdef AMREX_USE_HYPRE
    if (init_hypre) { HYPRE_Finalize(); }
#endif

    BL_TINY_PROFILE_FINALIZE();
    BL_PROFILE_FINALIZE();

#ifdef AMREX_USE_GPU
    amrex::DeallocateRandomSeedDevArray();
#endif

#ifdef BL_LAZY
    Lazy::Finalize();
#endif

    while (!The_Finalize_Function_Stack.empty())
    {
        //
        // Call the registered function.
        //
        The_Finalize_Function_Stack.top()();
        //
        // And then remove it from the stack.
        //
        The_Finalize_Function_Stack.pop();
    }

    // The MemPool stuff is not using The_Finalize_Function_Stack so that
    // it can be used in Fortran BoxLib.
#ifndef BL_AMRPROF
    if (amrex::system::verbose > 1)
    {
        int mp_min, mp_max, mp_tot;
        amrex_mempool_get_stats(mp_min, mp_max, mp_tot);  // in MB
        if (ParallelDescriptor::NProcs() == 1) {
            if (mp_tot > 0) {
                amrex::Print() << "MemPool: "
#ifdef AMREX_USE_OMP
                               << "min used in a thread: " << mp_min << " MB, "
                               << "max used in a thread: " << mp_max << " MB, "
#endif
                               << "tot used: " << mp_tot << " MB." << '\n';
            }
        } else {
            int global_max = mp_tot;
            int global_min = mp_tot;
            ParallelDescriptor::ReduceIntMax(global_max);
            if (global_max > 0) {
                ParallelDescriptor::ReduceIntMin(global_min);
                amrex::Print() << "MemPool: "
                               << "min used in a rank: " << global_min << " MB, "
                               << "max used in a rank: " << global_max << " MB.\n";
            }
        }
    }
#endif

#ifdef AMREX_MEM_PROFILING
    MemProfiler::report("Final");
    MemProfiler::Finalize();
#endif

#ifdef AMREX_USE_SUNDIALS
    sundials::Finalize();
#endif

    BL_TINY_PROFILE_MEMORYFINALIZE();

    amrex_mempool_finalize();
    Arena::Finalize();

#ifndef BL_AMRPROF
    if (system::signal_handling)
    {
        if (prev_handler_sigsegv != SIG_ERR) { std::signal(SIGSEGV, prev_handler_sigsegv); } // NOLINT(performance-no-int-to-ptr)
        if (prev_handler_sigterm != SIG_ERR) { std::signal(SIGTERM, prev_handler_sigterm); } // NOLINT(performance-no-int-to-ptr)
        if (prev_handler_sigint  != SIG_ERR) { std::signal(SIGINT , prev_handler_sigint);  } // NOLINT(performance-no-int-to-ptr)
        if (prev_handler_sigabrt != SIG_ERR) { std::signal(SIGABRT, prev_handler_sigabrt); } // NOLINT(performance-no-int-to-ptr)
        if (prev_handler_sigfpe  != SIG_ERR) { std::signal(SIGFPE , prev_handler_sigfpe);  } // NOLINT(performance-no-int-to-ptr)
        if (prev_handler_sigill  != SIG_ERR) { std::signal(SIGILL , prev_handler_sigill);  } // NOLINT(performance-no-int-to-ptr)
#if defined(__linux__)
#if !defined(__PGI) || (__PGIC__ >= 16)
        if (curr_fpe_excepts != 0) {
            fedisableexcept(curr_fpe_excepts);
            feenableexcept(prev_fpe_excepts);
        }
#endif
#elif defined(__APPLE__) && defined(__x86_64__)
        if (curr_fpe_excepts != 0u) {
            _MM_SET_EXCEPTION_MASK(prev_fpe_mask);
        }
#endif
    }
#endif

#ifdef AMREX_USE_GPU
    Gpu::Device::Finalize();
#endif

#ifdef AMREX_USE_OMP
    amrex::OpenMP::Finalize();
#endif

#if defined(AMREX_USE_UPCXX)
    upcxx::finalize();
#endif

    std::set_new_handler(prev_new_handler);

    amrex::OutStream().precision(prev_out_precision);
    amrex::ErrorStream().precision(prev_err_precision);

    bool is_ioproc = ParallelDescriptor::IOProcessor();

    /* Don't shut down MPI if GASNet is still using MPI */
#ifndef GASNET_CONDUIT_MPI
    ParallelDescriptor::EndParallel();
#endif

    if (amrex::system::verbose > 0 && is_ioproc) {
        amrex::OutStream() << "AMReX (" << amrex::Version() << ") finalized" << '\n';
    }
}

std::ostream&
amrex::OutStream ()
{
    return *system::osout;
}

std::ostream&
amrex::ErrorStream ()
{
    return *system::oserr;
}

std::string
amrex::get_command ()
{
    return command_line;
}

int
amrex::command_argument_count ()
{
    return static_cast<int>(command_arguments.size())-1;
}

std::string
amrex::get_command_argument (int number)
{
    if (number < static_cast<int>(command_arguments.size())) {
        return command_arguments[number];
    } else {
        return std::string();
    }
}

namespace amrex
{

AMReX::AMReX ()
    : m_geom(new Geometry())
{
}

AMReX::~AMReX ()
{
    delete m_geom;
}

void
AMReX::push (AMReX* pamrex)
{
    auto r = std::find_if(m_instance.begin(), m_instance.end(),
                          [=] (const std::unique_ptr<AMReX>& x) -> bool
                          { return x.get() == pamrex; });
    if (r == m_instance.end()) {
        m_instance.emplace_back(pamrex);
    } else if (r+1 != m_instance.end()) {
        std::rotate(r, r+1, m_instance.end());
    }
}

void
AMReX::erase (AMReX* pamrex)
{
    auto r = std::find_if(m_instance.begin(), m_instance.end(),
                          [=] (const std::unique_ptr<AMReX>& x) -> bool
                          { return x.get() == pamrex; });
    if (r != m_instance.end()) {
        m_instance.erase(r);
    }
}

FPExcept getFPExcept ()
{
    auto r = FPExcept::none;
#if defined(__linux__)
    auto excepts = fegetexcept();
    if (excepts & FE_INVALID  ) { r = r | FPExcept::invalid ; }
    if (excepts & FE_DIVBYZERO) { r = r | FPExcept::zero    ; }
    if (excepts & FE_OVERFLOW ) { r = r | FPExcept::overflow; }
#endif
    return r;
}

FPExcept setFPExcept (FPExcept excepts)
{
    auto prev = getFPExcept();
#if defined(__linux__)
    int flags = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW;
    fedisableexcept(flags);
    flags = 0;
    if (any(excepts & FPExcept::invalid )) { flags |= FE_INVALID  ; }
    if (any(excepts & FPExcept::zero    )) { flags |= FE_DIVBYZERO; }
    if (any(excepts & FPExcept::overflow)) { flags |= FE_OVERFLOW ; }
    feenableexcept(flags);
#else
    amrex::ignore_unused(excepts);
#endif
    return prev;
}

FPExcept disableFPExcept (FPExcept excepts)
{
    auto prev = getFPExcept();
#if defined(__linux__)
    int flags = 0;
    if (any(excepts & FPExcept::invalid )) { flags |= FE_INVALID  ; }
    if (any(excepts & FPExcept::zero    )) { flags |= FE_DIVBYZERO; }
    if (any(excepts & FPExcept::overflow)) { flags |= FE_OVERFLOW ; }
    fedisableexcept(flags);
#else
    amrex::ignore_unused(excepts);
#endif
    return prev;
}

FPExcept enableFPExcept (FPExcept excepts)
{
    auto prev = getFPExcept();
#if defined(__linux__)
    int flags = 0;
    if (any(excepts & FPExcept::invalid )) { flags |= FE_INVALID  ; }
    if (any(excepts & FPExcept::zero    )) { flags |= FE_DIVBYZERO; }
    if (any(excepts & FPExcept::overflow)) { flags |= FE_OVERFLOW ; }
    feenableexcept(flags);
#else
    amrex::ignore_unused(excepts);
#endif
    return prev;
}

}
