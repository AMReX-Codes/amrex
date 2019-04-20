#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <new>
#include <stack>
#include <limits>
#include <vector>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Box.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BLFort.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Arena.H>

#include <AMReX_Gpu.H>

#include <AMReX_Machine.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#endif

#ifndef BL_AMRPROF
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_VisMF.H>
#endif

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_BLBackTrace.H>
#include <AMReX_MemPool.H>

namespace amrex {
namespace system
{
    std::string exename;
    int verbose = 1;
    int signal_handling;
    int call_addr2line;
    int throw_exception;
    int regtest_reduction;
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
    int           prev_fpe_excepts;
    int           curr_fpe_excepts;
}

std::string amrex::Version ()
{
#ifdef AMREX_GIT_VERSION
    return std::string(AMREX_GIT_VERSION);
#else
    return std::string("Unknown");
#endif
}

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
	const char *cprocall = procall.str().c_str();
        const char * const end = " !!!\n";
	fwrite(cprocall, strlen(cprocall), 1, stderr);
        fwrite(str, strlen(str), 1, stderr);
        fwrite(end, strlen(end), 1, stderr);
    }
}

static
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

void
amrex::Error (const char* msg)
{
#if defined(__CUDA_ARCH__)
#if !defined(__APPLE__)
    if (msg) printf("%s\n", msg);
#endif
#else
    if (system::error_handler) {
        system::error_handler(msg);
    } else if (system::throw_exception) {
        throw RuntimeError(msg);
    } else {
        write_lib_id("Error");
        write_to_stderr_without_buffering(msg);
        ParallelDescriptor::Abort();
    }
#endif
}

void
amrex::Error (const std::string& msg)
{
    Error(msg.c_str());
}

void
amrex::Abort (const char* msg)
{
#if defined(__CUDA_ARCH__)
#if !defined(__APPLE__)
    if (msg) printf("Abort %s\n", msg);
    assert(0);
#endif
#else
    if (system::error_handler) {
        system::error_handler(msg);
    } else if (system::throw_exception) {
        throw RuntimeError(msg);
    } else {
       write_lib_id("Abort");
       write_to_stderr_without_buffering(msg);
#ifdef _OPENMP
#pragma omp critical (amrex_abort_omp_critical)
#endif
       ParallelDescriptor::Abort();
   }
#endif
}

void
amrex::Abort (const std::string& msg)
{
    Abort(msg.c_str());
}

void
amrex::Warning (const char* msg)
{
#if defined(__CUDA_ARCH__)
#if !defined(__APPLE__)
    if (msg) printf("%s\n", msg);
#endif
#else
    if (msg)
    {
	amrex::Print(Print::AllProcs,amrex::ErrorStream()) << msg << '!' << '\n';
    }
#endif
}

void
amrex::Warning (const std::string& msg)
{
    Warning(msg.c_str());
}

void
amrex::Assert (const char* EX,
               const char* file,
               int         line,
               const char* msg)
{
#if defined(__CUDA_ARCH__)
#if !defined(__APPLE__)
    if (msg) {
        printf("Assertion `%s' failed, file \"%s\", line %d, Msg: %s",
               EX, file, line, msg);
    } else {
        printf("Assertion `%s' failed, file \"%s\", line %d",
               EX, file, line);
    }

    assert(0);
#endif
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
       ParallelDescriptor::Abort();
   }
#endif
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

void
amrex::Initialize (MPI_Comm mpi_comm,
                   std::ostream& a_osout, std::ostream& a_oserr,
                   ErrorHandler a_errhandler)
{
    int argc = 0;
    char** argv = 0;
    Initialize(argc, argv, false, mpi_comm, {}, a_osout, a_oserr, a_errhandler);
}

void
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
            constexpr int bufSize = 1024;
            char temp[bufSize];
            char *rCheck = getcwd(temp, bufSize);
            if(rCheck == 0) {
                amrex::Abort("**** Error:  getcwd buffer too small.");
            }
            system::exename = temp;
            system::exename += "/";
        }
        system::exename += argv[0];

        for (int i = 0; i < argc; ++i) {
            if (i != 0) command_line.append(" ");
            command_line.append(argv[i]);
            command_arguments.push_back(std::string(argv[i]));
        }
    }

#if defined(PERILLA_USE_UPCXX) || defined(AMREX_USE_UPCXX)
    upcxx::init();
#endif

#ifdef BL_USE_MPI3
    BL_MPI_REQUIRE( MPI_Win_create_dynamic(MPI_INFO_NULL, ParallelDescriptor::Communicator(),
                                           &ParallelDescriptor::cp_win) );
    BL_MPI_REQUIRE( MPI_Win_create_dynamic(MPI_INFO_NULL, ParallelDescriptor::Communicator(),
                                           &ParallelDescriptor::fb_win) );
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
            if (strchr(argv[1],'='))
            {
                ParmParse::Initialize(argc-1,argv+1,0);
            }
            else
            {
                ParmParse::Initialize(argc-2,argv+2,argv[1]);
            }
        }
    } else {
        ParmParse::Initialize(0,0,0);
    }

    if (func_parm_parse) {
        func_parm_parse();
    }

#ifdef AMREX_USE_GPU
    // Initialize after ParmParse so that we can read inputs.
    Gpu::Device::Initialize();
#endif

    {
	ParmParse pp("amrex");
	pp.query("v", system::verbose);
	pp.query("verbose", system::verbose);
        pp.query("regtest_reduction", system::regtest_reduction);
        pp.query("signal_handling", system::signal_handling);
        pp.query("throw_exception", system::throw_exception);
        pp.query("call_addr2line", system::call_addr2line);

        if (system::signal_handling)
        {
            // We could save the singal handlers and restore them in Finalize.
            prev_handler_sigsegv = signal(SIGSEGV, BLBackTrace::handler); // catch seg fault
            prev_handler_sigint = signal(SIGINT,  BLBackTrace::handler);
            prev_handler_sigabrt = signal(SIGABRT, BLBackTrace::handler);

            int term = 0;
            pp.query("handle_sigterm", term);
            if (term) {
                prev_handler_sigterm = signal(SIGTERM,  BLBackTrace::handler);
            } else {
                prev_handler_sigterm = SIG_ERR;
            }

            prev_handler_sigfpe = SIG_ERR;

            int invalid = 0, divbyzero=0, overflow=0;
            pp.query("fpe_trap_invalid", invalid);
            pp.query("fpe_trap_zero", divbyzero);
            pp.query("fpe_trap_overflow", overflow);
            curr_fpe_excepts = 0;
            if (invalid)   curr_fpe_excepts |= FE_INVALID;
            if (divbyzero) curr_fpe_excepts |= FE_DIVBYZERO;
            if (overflow)  curr_fpe_excepts |= FE_OVERFLOW;
#if defined(__linux__) && !defined(__NEC__)
#if !defined(__PGI) || (__PGIC__ >= 16)
            prev_fpe_excepts = fegetexcept();
            if (curr_fpe_excepts != 0) {
                feenableexcept(curr_fpe_excepts);  // trap floating point exceptions
                prev_handler_sigfpe = signal(SIGFPE,  BLBackTrace::handler);
            }
#endif
#endif
        }
    }

    ParallelDescriptor::Initialize();

    //
    // Initialize random seed after we're running in parallel.
    //
    amrex::InitRandom(ParallelDescriptor::MyProc()+1, ParallelDescriptor::NProcs());

    Arena::Initialize();
    amrex_mempool_init();

    // For thread safety, we should do these initializations here.
    BaseFab_Initialize();
    BoxArray::Initialize();
    DistributionMapping::Initialize();
    FArrayBox::Initialize();
    IArrayBox::Initialize();
    Gpu::AsyncFab::Initialize();
    FabArrayBase::Initialize();
    MultiFab::Initialize();
    iMultiFab::Initialize();
    VisMF::Initialize();
#ifdef AMREX_USE_EB
    EB2::Initialize();
#endif
    BL_PROFILE_INITPARAMS();
#endif
    machine::Initialize();

    if (double(std::numeric_limits<long>::max()) < 9.e18)
    {
        if (system::verbose) {
            amrex::Print() << "!\n! WARNING: Maximum of long int, "
                           << std::numeric_limits<long>::max() 
                           << ", might be too small for big runs.\n!\n";
        }
    }

    if (system::verbose > 0)
    {
#ifdef BL_USE_MPI
        amrex::Print() << "MPI initialized with "
                       << ParallelDescriptor::NProcs()
                       << " MPI processes\n";
#endif
        
#ifdef _OPENMP
//    static_assert(_OPENMP >= 201107, "OpenMP >= 3.1 is required.");
        amrex::Print() << "OMP initialized with "
                       << omp_get_max_threads()
                       << " OMP threads\n";
#endif

        amrex::Print() << "AMReX (" << amrex::Version() << ") initialized" << std::endl;
    }

    BL_TINY_PROFILE_INITIALIZE();
}

void
amrex::Finalize (bool finalize_parallel)
{
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
#ifdef _OPENMP
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

#ifdef BL_MEM_PROFILING
    MemProfiler::report("Final");
    MemProfiler::Finalize();
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
#if defined(__linux__) && !defined(__NEC__)
#if !defined(__PGI) || (__PGIC__ >= 16)
        if (curr_fpe_excepts != 0) {
            fedisableexcept(curr_fpe_excepts);
            feenableexcept(prev_fpe_excepts);
        }
#endif
#endif
    }
#endif

#ifdef AMREX_USE_GPU
    Gpu::Device::Finalize();
#endif

#if defined(PERILLA_USE_UPCXX) || defined(AMREX_USE_UPCXX)
    upcxx::finalize();
#endif

    std::set_new_handler(prev_new_handler);

    amrex::OutStream().precision(prev_out_precision);
    amrex::ErrorStream().precision(prev_err_precision);

    bool is_ioproc = ParallelDescriptor::IOProcessor();

    if (finalize_parallel) {
    /* Don't shut down MPI if GASNet is still using MPI */
#ifndef GASNET_CONDUIT_MPI
        ParallelDescriptor::EndParallel();
#endif
    }

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
