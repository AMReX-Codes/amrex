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
#include <stack>
#include <limits>
#include <vector>
#include <algorithm>

namespace amrex {

std::vector<std::unique_ptr<AMReX> > AMReX::m_instance;

namespace system
{
    std::string exename;
    int verbose = 1;
    int signal_handling;
    int call_addr2line;
    int throw_exception;
    int regtest_reduction;
    int abort_on_unused_inputs = 0;
    std::ostream* osout = &std::cout;
    std::ostream* oserr = &std::cerr;
    ErrorHandler error_handler = nullptr;
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
    typedef void (*SignalHandler)(int);
    SignalHandler prev_handler_sigsegv;
    SignalHandler prev_handler_sigterm;
    SignalHandler prev_handler_sigint;
    SignalHandler prev_handler_sigabrt;
    SignalHandler prev_handler_sigfpe;
#if defined(__linux__)
    int           prev_fpe_excepts;
    int           curr_fpe_excepts;
#elif defined(__APPLE__) && defined(__x86_64__)
    unsigned int  prev_fpe_mask;
    unsigned int  curr_fpe_excepts;
#endif
}

#ifdef AMREX_USE_HYPRE
namespace {
    int init_hypre = 1;
#if defined(HYPRE_USING_CUDA) || defined(HYPRE_USING_HIP)
    int hypre_spgemm_use_vendor = 0;
    int hypre_spmv_use_vendor = 0;
    int hypre_sptrans_use_vendor = 0;
#endif
}
#endif

int amrex::Verbose () noexcept { return amrex::system::verbose; }

void amrex::SetVerbose (int v) noexcept { amrex::system::verbose = v; }

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
    fflush(NULL);

    if (str)
    {
        std::ostringstream procall;
        procall << ParallelDescriptor::MyProc() << "::";
        auto tmp = procall.str();
        const char *cprocall = tmp.c_str();
        const char * const end = " !!!\n";
        fwrite(cprocall, strlen(cprocall), 1, stderr);
        fwrite(str, strlen(str), 1, stderr);
        fwrite(end, strlen(end), 1, stderr);
    }
}

namespace {
void
write_lib_id(const char* msg)
{
    fflush(0);
    const char* const s = "amrex::";
    fwrite(s, strlen(s), 1, stderr);
    if ( msg )
    {
        fwrite(msg, strlen(msg), 1, stderr);
        fwrite("::", 2, 1, stderr);
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
amrex::Error_host (const char * msg)
{
    if (system::error_handler) {
        system::error_handler(msg);
    } else if (system::throw_exception) {
        throw RuntimeError(msg);
    } else {
        write_lib_id("Error");
        write_to_stderr_without_buffering(msg);
#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_abort_omp_critical)
#endif
        ParallelDescriptor::Abort();
    }
}

void
amrex::Warning_host (const char * msg)
{
    if (msg) {
        amrex::Print(Print::AllProcs,amrex::ErrorStream()) << msg << '!' << '\n';
    }
}

void
amrex::Abort_host (const char * msg)
{
    if (system::error_handler) {
        system::error_handler(msg);
    } else if (system::throw_exception) {
        throw RuntimeError(msg);
    } else {
       write_lib_id("Abort");
       write_to_stderr_without_buffering(msg);
#ifdef AMREX_USE_OMP
#pragma omp critical (amrex_abort_omp_critical)
#endif
       ParallelDescriptor::Abort();
   }
}

void
amrex::Assert_host (const char* EX, const char* file, int line, const char* msg)
{
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
}

namespace
{
    std::stack<amrex::PTR_TO_VOID_FUNC> The_Finalize_Function_Stack;
    std::stack<amrex::PTR_TO_VOID_FUNC> The_Initialize_Function_Stack;
}

void
amrex::ExecOnFinalize (PTR_TO_VOID_FUNC fp)
{
    The_Finalize_Function_Stack.push(fp);
}

void
amrex::ExecOnInitialize (PTR_TO_VOID_FUNC fp)
{
    The_Initialize_Function_Stack.push(fp);
}

amrex::AMReX*
amrex::Initialize (MPI_Comm mpi_comm,
                   std::ostream& a_osout, std::ostream& a_oserr,
                   ErrorHandler a_errhandler)
{
    int argc = 0;
    char** argv = 0;
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
    system::regtest_reduction = 0;
    system::signal_handling = 1;
    system::call_addr2line = 1;
    system::throw_exception = 0;
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
            if (i != 0) command_line.append(" ");
            command_line.append(argv[i]);
            command_arguments.push_back(std::string(argv[i]));
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
        (*The_Initialize_Function_Stack.top())();
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
            ParmParse::Initialize(0,0,0);
        }
        else if (argc > 1)
        {
            if (argv[1][0] == '-')
            {
                // If arguments list starts with "-", do not use ParmParse.
                // Application code can then parse the command line. This will
                // prevent "-h" or "--help" from creating errors in ParmParse,
                // but only if it's the first argument after the executable.
                ParmParse::Initialize(0,0,0);
            }
            else
            {
                // This counts command line arguments before a "--"
                // and only sends the preceeding arguments to ParmParse;
                // the rest get ingored.
                int ppargc = 1;
                for (; ppargc < argc; ++ppargc) {
                    if (strcmp(argv[ppargc], "--") == 0) break;
                }
                if (ppargc > 1)
                {
                    if (strchr(argv[1],'=') || (argc > 2 ? argv[2][0] == '=' : false) )
                    {
                        // No inputs file to parse
                        ParmParse::Initialize(ppargc-1,argv+1,0);
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
        ParmParse::Initialize(0,0,0);
    }

    if (func_parm_parse) {
        func_parm_parse();
    }

    {
        ParmParse pp("amrex");
        pp.queryAdd("v", system::verbose);
        pp.queryAdd("verbose", system::verbose);
    }

#ifdef AMREX_USE_GPU
    // Initialize after ParmParse so that we can read inputs.
    Gpu::Device::Initialize();
#ifdef AMREX_USE_CUPTI
    CuptiInitialize();
#endif
#endif

    {
        ParmParse pp("amrex");
        pp.queryAdd("regtest_reduction", system::regtest_reduction);
        pp.queryAdd("signal_handling", system::signal_handling);
        pp.queryAdd("throw_exception", system::throw_exception);
        pp.queryAdd("call_addr2line", system::call_addr2line);
        pp.queryAdd("abort_on_unused_inputs", system::abort_on_unused_inputs);

        if (system::signal_handling)
        {
            // We could save the singal handlers and restore them in Finalize.
            prev_handler_sigsegv = signal(SIGSEGV, BLBackTrace::handler); // catch seg fault
            prev_handler_sigint = signal(SIGINT,  BLBackTrace::handler);
            prev_handler_sigabrt = signal(SIGABRT, BLBackTrace::handler);

            int term = 0;
            pp.queryAdd("handle_sigterm", term);
            if (term) {
                prev_handler_sigterm = signal(SIGTERM,  BLBackTrace::handler);
            } else {
                prev_handler_sigterm = SIG_ERR;
            }

            prev_handler_sigfpe = SIG_ERR;

            int invalid = 0, divbyzero=0, overflow=0;
            pp.queryAdd("fpe_trap_invalid", invalid);
            pp.queryAdd("fpe_trap_zero", divbyzero);
            pp.queryAdd("fpe_trap_overflow", overflow);

#if defined(__linux__)
            curr_fpe_excepts = 0;
            if (invalid)   curr_fpe_excepts |= FE_INVALID;
            if (divbyzero) curr_fpe_excepts |= FE_DIVBYZERO;
            if (overflow)  curr_fpe_excepts |= FE_OVERFLOW;
#if !defined(AMREX_USE_DPCPP) && (!defined(__PGI) || (__PGIC__ >= 16))
            // xxxxx DPCPP todo: fpe trap
            prev_fpe_excepts = fegetexcept();
            if (curr_fpe_excepts != 0) {
                feenableexcept(curr_fpe_excepts);  // trap floating point exceptions
                prev_handler_sigfpe = signal(SIGFPE,  BLBackTrace::handler);
            }
#endif

#elif defined(__APPLE__) && defined(__x86_64__)
            prev_fpe_mask = _MM_GET_EXCEPTION_MASK();
            curr_fpe_excepts = 0u;
            if (invalid)   curr_fpe_excepts |= _MM_MASK_INVALID;
            if (divbyzero) curr_fpe_excepts |= _MM_MASK_DIV_ZERO;
            if (overflow)  curr_fpe_excepts |= _MM_MASK_OVERFLOW;
            if (curr_fpe_excepts != 0u) {
                _MM_SET_EXCEPTION_MASK(prev_fpe_mask & ~curr_fpe_excepts);
                prev_handler_sigfpe = signal(SIGFPE,  BLBackTrace::handler);
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

#ifdef AMREX_USE_EB
    EB2::Initialize();
#endif

    BL_PROFILE_INITPARAMS();
#endif // ifndef BL_AMRPROF

    machine::Initialize();

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

    if (system::verbose > 0)
    {
#ifdef BL_USE_MPI

        amrex::Print() << "MPI initialized with "
                       << ParallelDescriptor::NProcs()
                       << " MPI processes\n";

        int provided = -1;
        MPI_Query_thread(&provided);
        amrex::Print() << "MPI initialized with thread support level " << provided << std::endl;
#endif

#ifdef AMREX_USE_OMP
//    static_assert(_OPENMP >= 201107, "OpenMP >= 3.1 is required.");
        amrex::Print() << "OMP initialized with "
                       << omp_get_max_threads()
                       << " OMP threads\n";
#endif

        amrex::Print() << "AMReX (" << amrex::Version() << ") initialized" << std::endl;
    }

    BL_TINY_PROFILE_INITIALIZE();

    AMReX::push(new AMReX());
    return AMReX::top();
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
    if (init_hypre) HYPRE_Finalize();
#endif

    BL_TINY_PROFILE_FINALIZE();
    BL_PROFILE_FINALIZE();

#ifdef AMREX_USE_CUDA
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
        (*The_Finalize_Function_Stack.top())();
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
                               << "tot used: " << mp_tot << " MB." << std::endl;
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

    amrex_mempool_finalize();
    Arena::Finalize();

#ifndef BL_AMRPROF
    if (system::signal_handling)
    {
        if (prev_handler_sigsegv != SIG_ERR) signal(SIGSEGV, prev_handler_sigsegv);
        if (prev_handler_sigterm != SIG_ERR) signal(SIGTERM, prev_handler_sigterm);
        if (prev_handler_sigint != SIG_ERR) signal(SIGINT, prev_handler_sigint);
        if (prev_handler_sigabrt != SIG_ERR) signal(SIGABRT, prev_handler_sigabrt);
        if (prev_handler_sigfpe != SIG_ERR) signal(SIGFPE, prev_handler_sigfpe);
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
        amrex::OutStream() << "AMReX (" << amrex::Version() << ") finalized" << std::endl;
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
    return command_arguments.size()-1;
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

}
