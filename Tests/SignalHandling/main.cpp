#include <AMReX.H>
#include <AMReX_BLBackTrace.H>

#include <csignal>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#if defined(__linux__) || defined(__APPLE__)
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

// AMReX writes Backtrace files on Linux with glibc and on macOS.
#if ((defined(__linux__) && defined(__GLIBC__)) || defined(__APPLE__)) && \
    !defined(AMREX_USE_MPI) && !defined(AMREX_USE_GPU)
#define CHECK_BACKTRACE_FILE 1
#endif
#if defined(CHECK_BACKTRACE_FILE) && defined(__APPLE__)
extern char** environ; // glibc declares it in unistd.h
#endif

// amrex::Initialize installs the signal handlers of AMReX, and
// amrex::Finalize restores the previous ones with their flags and gives back
// the alternate signal stack, also after nested calls. A crash writes a
// Backtrace file.

namespace {

int failures = 0;

#if defined(__linux__) || defined(__APPLE__)
void check (bool ok, char const* what)
{
    if (!ok) {
        std::fprintf(stderr, "FAILED: %s\n", what);
        ++failures;
    }
}

// A fault after it was called ends the process with the default action.
void previous_handler (int sig, siginfo_t*, void*) { std::signal(sig, SIG_DFL); }

bool is_previous_handler (int sig)
{
    struct sigaction cur{};
    sigaction(sig, nullptr, &cur);
    return (cur.sa_flags & SA_SIGINFO) && cur.sa_sigaction == previous_handler;
}

bool has_alt_stack ()
{
    stack_t st{};
    return sigaltstack(nullptr, &st) == 0 && (st.ss_flags & SS_DISABLE) == 0;
}
#endif

#ifdef CHECK_BACKTRACE_FILE
char const* crash_arg = "signal_handling_test.crash=1";

// In a new process of this program, a segfault on two threads at once after
// amrex::Initialize writes a Backtrace file and ends the process with an
// error. On macOS, the frames are of the code that raised the signal, not of
// the handler.
void check_backtrace_file (int argc, char* argv[])
{
#ifdef AMREX_USE_OMP
    const std::string file = "Backtrace.0.0";
#else
    const std::string file = "Backtrace.0";
#endif
    std::remove(file.c_str());
    std::vector<char*> args(argv, argv + argc);
    args.push_back(const_cast<char*>(crash_arg));
    args.push_back(nullptr);
    std::cout.flush(); // before the output of the new process
    std::fflush(nullptr);
    pid_t pid = 0;
    int status = 0;
    const bool ran = posix_spawnp(&pid, argv[0], nullptr, nullptr, args.data(), environ) == 0;
    check(ran, "run this program in a new process");
    if (!ran) { return; }
    // at most five minutes
    pid_t r = 0;
    for (int i = 0; i < 3000 && (r = waitpid(pid, &status, WNOHANG)) == 0; ++i) {
        usleep(100000);
    }
    if (r == 0) {
        kill(pid, SIGKILL);
        waitpid(pid, &status, 0);
    }
    check(r == pid, "the segfault ends the process within five minutes");
    if (r != pid) { return; }
    check(!(WIFEXITED(status) && WEXITSTATUS(status) == 0), "the segfault ends the process with an error");
    std::ifstream ifs(file);
    std::stringstream ss;
    ss << ifs.rdbuf();
    const std::string bt = ss.str();
    check(!bt.empty(), "the segfault writes a Backtrace file");
#ifdef __APPLE__
    check(bt.find(" 1: [") != std::string::npos && bt.find("handle_signal") == std::string::npos,
          "the Backtrace file has the frames of the code that raised the signal");
#endif
    std::remove(file.c_str());
}
#endif

}

int main (int argc, char* argv[])
{
#ifdef CHECK_BACKTRACE_FILE
    if (argc > 1 && std::string(argv[argc-1]) == crash_arg) {
        amrex::Initialize(argc, argv);
        std::thread t([] { std::raise(SIGSEGV); });
        std::raise(SIGSEGV);
        t.join();
        return 0;
    }
#endif

#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

#if defined(__linux__) || defined(__APPLE__)
    int const signals[] = {SIGSEGV, SIGBUS, SIGFPE, SIGILL, SIGTRAP, SIGABRT, SIGINT, SIGTERM};
    struct sigaction sa{};
    sa.sa_sigaction = previous_handler;
    sa.sa_flags = SA_SIGINFO;
    sigemptyset(&sa.sa_mask);
    for (int sig : signals) { sigaction(sig, &sa, nullptr); }
    const bool had_alt_stack = has_alt_stack();
#endif
    bool handled = false;

    for (int cycle = 0; cycle < 2; ++cycle) {
        amrex::Initialize(argc, argv);
        // Not, e.g., for Intel GPUs, whose runtime handles SIGSEGV.
        handled = amrex::system::signal_handling && amrex::system::handle_sigsegv;
#if defined(__linux__) || defined(__APPLE__)
        if (handled) {
            check(!is_previous_handler(SIGSEGV), "Initialize installs the SIGSEGV handler of AMReX");
            amrex::BLBackTrace::install_signal_handlers();
            amrex::BLBackTrace::restore_signal_handlers();
            check(!is_previous_handler(SIGSEGV), "a nested restore keeps the handler of AMReX");
        }
#endif
        amrex::Finalize();
#if defined(__linux__) || defined(__APPLE__)
        if (handled) {
            bool restored = true;
            for (int sig : signals) { restored = restored && is_previous_handler(sig); }
            check(restored, "Finalize restores the previous signal handlers with SA_SIGINFO");
            check(has_alt_stack() == had_alt_stack, "Finalize gives back the alternate signal stack");
        }
#endif
    }

#ifdef CHECK_BACKTRACE_FILE
    if (handled) { check_backtrace_file(argc, argv); }
#endif

    amrex::ignore_unused(handled);

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif

    if (failures == 0) { std::printf("SignalHandling test passed\n"); }
    return failures == 0 ? 0 : 1;
}
