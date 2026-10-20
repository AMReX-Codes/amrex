#include <AMReX_BLBackTrace.H>
#include <AMReX_CrtReport.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_AsyncOut.H>
#include <AMReX.H>
#include <AMReX_Utility.H>
#ifdef AMREX_USE_MPI
#include <AMReX_MPMD.H>
#endif

#ifdef AMREX_TINY_PROFILING
#include <AMReX_TinyProfiler.H>
#endif

#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#if !(defined(_MSC_VER) && defined(__CUDACC__))
//MSVC can't pre-processor cfenv with `Zc:preprocessor`
//https://developercommunity.visualstudio.com/content/problem/1271183/zcpreprocessor-e-crashes-when-given.html
#include <cfenv>
#endif

#if defined(AMREX_EXPORT_DYNAMIC) && defined(__APPLE__)
#include <cxxabi.h>
#include <dlfcn.h>
#define AMREX_BACKTRACE_SUPPORTED 1
#elif defined(__linux__) && defined(__GLIBC__)
#include <link.h>
#define AMREX_BACKTRACE_SUPPORTED 1
#endif

#ifndef _WIN32
#  if defined(__GLIBC__) || defined(__APPLE__)
#    include <execinfo.h>
#  endif
#include <sys/mman.h>
#include <unistd.h>
#endif

namespace amrex {

std::stack<std::pair<std::string, std::string> >  BLBackTrace::bt_stack;

namespace {
// A signal handler that was in place before install_signal_handlers.
struct PrevSignalHandler
{
    bool saved = false;
#ifdef _WIN32
    void (*handler)(int) = SIG_DFL;
#else
    struct sigaction action{};
#endif
};

// The signals that AMReX handles if their flag is true.
struct HandledSignal
{
    int sig;
    bool* handle;
    char const* message;
    PrevSignalHandler prev;
};

HandledSignal g_signals[] = {
    {.sig = SIGSEGV, .handle = &amrex::system::handle_sigsegv, .message = "Segfault", .prev = {}},
    {.sig = SIGFPE , .handle = &amrex::system::handle_sigfpe , .message = "Erroneous arithmetic operation", .prev = {}},
    {.sig = SIGILL , .handle = &amrex::system::handle_sigill , .message = "SIGILL Invalid, privileged, or ill-formed instruction", .prev = {}},
    {.sig = SIGTERM, .handle = &amrex::system::handle_sigterm, .message = "SIGTERM", .prev = {}},
    {.sig = SIGINT , .handle = &amrex::system::handle_sigint , .message = "SIGINT", .prev = {}},
    {.sig = SIGABRT, .handle = &amrex::system::handle_sigabrt, .message = "SIGABRT", .prev = {}},
#if defined(SIGBUS) && !defined(_WIN32)
    {.sig = SIGBUS , .handle = &amrex::system::handle_sigbus , .message = "SIGBUS Bus error", .prev = {}},
#endif
#if defined(SIGTRAP) && !defined(_WIN32)
    {.sig = SIGTRAP, .handle = &amrex::system::handle_sigtrap, .message = "SIGTRAP Trace or breakpoint trap", .prev = {}},
#endif
};

// How often install_signal_handlers was called without
// restore_signal_handlers (e.g., by nested amrex::Initialize), and at which
// of these depths it installed the handlers, or 0.
int g_install_depth = 0;
int g_installed_depth = 0;

#ifndef _WIN32
constexpr std::size_t g_alt_stack_size = std::size_t(1) << 20;

// Disables the alternate stack sp of this thread unless it is in use.
// Returns whether this thread does not use sp (anymore).
bool give_back_alt_stack (void* sp)
{
    stack_t cur{};
    if (sigaltstack(nullptr, &cur) != 0) { return false; }
    if ((cur.ss_flags & SS_DISABLE) != 0 || cur.ss_sp != sp) { return true; }
    if ((cur.ss_flags & SS_ONSTACK) != 0) { return false; }
    // macOS rejects SS_DISABLE with a size below MINSIGSTKSZ.
    stack_t ss{};
    ss.ss_sp = sp;
    ss.ss_size = g_alt_stack_size;
    ss.ss_flags = SS_DISABLE;
    return sigaltstack(&ss, nullptr) == 0;
}

// An alternate stack for the signal handler on a thread that installs it and
// has none, so that the handler can run after a stack overflow there. It is
// mapped below a guard page, and unmapped when the thread exits.
struct AltStack
{
    void* sp = nullptr;
    AltStack () = default;
    AltStack (AltStack const&) = delete;
    AltStack (AltStack&&) = delete;
    AltStack& operator= (AltStack const&) = delete;
    AltStack& operator= (AltStack&&) = delete;
    ~AltStack ()
    {
        if (sp != nullptr && give_back_alt_stack(sp)) {
            const auto page = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
            munmap(static_cast<char*>(sp) - page, g_alt_stack_size + page);
        }
    }
};
thread_local AltStack t_alt_stack;

// The stack left if this runs on an alternate signal stack, which can be
// small (e.g., one of another component on its own thread), else SIZE_MAX.
std::size_t alt_stack_left ()
{
    stack_t st{};
    if (sigaltstack(nullptr, &st) != 0 || (st.ss_flags & SS_ONSTACK) == 0) { return SIZE_MAX; }
    const char here = 0;
    const auto sp = reinterpret_cast<std::uintptr_t>(&here);
    const auto lo = reinterpret_cast<std::uintptr_t>(st.ss_sp);
    return (sp > lo) ? sp - lo : 0;
}
#endif

HandledSignal const* find_signal (int s)
{
    for (auto const& h : g_signals) {
        if (h.sig == s) { return &h; }
    }
    return nullptr;
}

void
handle_signal (int s)
{
    std::signal(s, SIG_DFL);

#ifndef _WIN32
    // Writing the Backtrace file could overflow a small alternate stack. It
    // needs about 12 KiB on macOS, more on Linux (e.g., for backtrace()).
#ifdef __APPLE__
    constexpr std::size_t stack_needed = std::size_t(1) << 14;
#else
    constexpr std::size_t stack_needed = std::size_t(1) << 16;
#endif
    if (alt_stack_left() < stack_needed) {
        std::raise(s); // with the default action once this handler returns
        return;
    }
#endif

    AsyncOut::Finalize();

    if (HandledSignal const* h = find_signal(s)) {
        amrex::ErrorStream() << h->message << '\n';
    }

#if defined(AMREX_BACKTRACE_SUPPORTED) || defined(AMREX_TINY_PROFILING)

    auto& bt_stack = BLBackTrace::bt_stack;

    std::string errfilename;
    {
        std::ostringstream ss;
#ifdef AMREX_USE_MPI
        if (MPMD::Initialized()) {
            ss << "Backtrace.prog" << MPMD::MyProgId() << ".";
        } else
#endif
        {
            ss << "Backtrace.";
        }
        ss << ParallelDescriptor::MyProc();
#ifdef AMREX_USE_OMP
        ss << "." << omp_get_thread_num();
#endif
        errfilename = ss.str();
    }

    if (FILE* p = fopen(errfilename.c_str(), "w")) {
#if defined(AMREX_BACKTRACE_SUPPORTED)
        BLBackTrace::print_backtrace_info(p);
#endif
        fclose(p);
    }

    amrex::ErrorStream() << "See " << errfilename << " file for details" << '\n';

    if (!bt_stack.empty()) {
        std::ofstream errfile;
        errfile.open(errfilename.c_str(), std::ofstream::out | std::ofstream::app);
        if (errfile.is_open()) {
            errfile << '\n';
            while (!bt_stack.empty()) {
                errfile << "== BACKTRACE == " << bt_stack.top().first
                        <<", " << bt_stack.top().second << "\n";
                bt_stack.pop();
            }
            errfile << '\n';
        }
    }

#ifdef AMREX_TINY_PROFILING
    {
        std::ofstream errfile;
        errfile.open(errfilename.c_str(), std::ofstream::out | std::ofstream::app);
        if (errfile.is_open()) {
            errfile << '\n';
            TinyProfiler::PrintCallStack(errfile);
            errfile << '\n';
        }
    }
#endif

    if (ParallelDescriptor::NProcs() > 1) {
        amrex::Sleep(3);
    }

#endif

    ParallelDescriptor::Abort(s, false);
}

#ifdef _WIN32
void signal_handler (int s)
{
    handle_signal(s);
}
#else
void signal_handler (int s, siginfo_t* /*info*/, void* /*uctx*/)
{
    handle_signal(s);
}
#endif
}

void
BLBackTrace::handler (int s)
{
    handle_signal(s);
}

void
BLBackTrace::install_signal_handlers (bool install)
{
    ++g_install_depth;
    if (!install || g_installed_depth != 0) { return; }
    g_installed_depth = g_install_depth;

#if defined(__linux__) && defined(AMREX_BACKTRACE_SUPPORTED)
    // glibc loads libgcc_s on the first call of backtrace(), which should not
    // happen in the signal handler, e.g., with little stack left.
    {
        void* frame = nullptr;
        backtrace(&frame, 1);
    }
#endif

#ifndef _WIN32
    // Use the alternate stack only if the thread has none of its own (e.g.,
    // of Python's faulthandler, or a sanitizer), which may be too small, and
    // only if SIGSEGV or SIGBUS is handled, which a stack overflow raises.
    // These then always run on it on this thread. Other threads (e.g., of
    // OpenMP) do not get one.
    bool use_alt_stack = false;
    if (amrex::system::handle_sigsegv || amrex::system::handle_sigbus)
    {
        stack_t cur{};
        if (sigaltstack(nullptr, &cur) == 0) {
            auto& alt = t_alt_stack.sp;
            if ((cur.ss_flags & SS_DISABLE) == 0) {
                use_alt_stack = (alt != nullptr && cur.ss_sp == alt);
            } else {
                if (alt == nullptr) {
                    const auto page = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
                    void* p = mmap(nullptr, g_alt_stack_size + page, PROT_READ | PROT_WRITE,
                                   MAP_PRIVATE | MAP_ANON, -1, 0);
                    if (p != MAP_FAILED) {
                        mprotect(p, page, PROT_NONE); // guard page below the stack
                        alt = static_cast<char*>(p) + page;
                    }
                }
                if (alt != nullptr) {
                    stack_t ss{};
                    ss.ss_sp = alt;
                    ss.ss_size = g_alt_stack_size;
                    use_alt_stack = (sigaltstack(&ss, nullptr) == 0);
                }
            }
        }
    }
#endif

    for (auto& h : g_signals) {
        if (!*h.handle) { continue; }
#ifdef _WIN32
        auto* prev = std::signal(h.sig, signal_handler);
        h.prev.saved = (prev != SIG_ERR); // NOLINT(performance-no-int-to-ptr)
        h.prev.handler = prev;
#else
        struct sigaction sa{};
        sa.sa_sigaction = signal_handler;
        sa.sa_flags = SA_SIGINFO | SA_RESTART;
#ifdef SIGBUS
        // A stack overflow raises SIGSEGV or SIGBUS. Other signals do not use
        // the alternate stack. On other threads, these use the alternate stack
        // of the thread if it has one, which can be small (e.g., of another
        // component), so the handler ends the process right away if little
        // of it is left.
        if (use_alt_stack && (h.sig == SIGSEGV || h.sig == SIGBUS)) {
            sa.sa_flags |= SA_ONSTACK;
        }
#endif
        sigemptyset(&sa.sa_mask);
        h.prev.saved = (sigaction(h.sig, &sa, &h.prev.action) == 0);
#endif
        if (h.sig == SIGABRT && h.prev.saved && amrex::system::handle_crt_reports) {
            // abort() raises SIGABRT for this handler without a CRT report
            detail::CrtReportSkipAbortReport(true);
        }
    }
}

void
BLBackTrace::restore_signal_handlers ()
{
    if (g_install_depth == 0) { return; }
    const bool installed = (g_installed_depth == g_install_depth);
    --g_install_depth;
    if (!installed) { return; }
    g_installed_depth = 0;

    for (auto& h : g_signals) {
        if (!h.prev.saved) { continue; }
        if (h.sig == SIGABRT) {
            // restore the CRT report of abort() together with the previous handler
            detail::CrtReportSkipAbortReport(false);
        }
#ifdef _WIN32
        std::signal(h.sig, h.prev.handler);
#else
        sigaction(h.sig, &h.prev.action, nullptr);
#endif
        h.prev.saved = false;
    }

#ifndef _WIN32
    // Other threads keep theirs until they exit.
    if (t_alt_stack.sp != nullptr) { give_back_alt_stack(t_alt_stack.sp); }
#endif
}

void
BLBackTrace::print_backtrace_info (const std::string& filename)
{
    if (FILE* p = fopen(filename.c_str(), "w"))
    {
      BLBackTrace::print_backtrace_info(p);
      fclose(p);
    }
    else
    {
        amrex::Print() << "Warning @ BLBackTrace::print_backtrace_info: "
                       << filename << " is not a valid output file."
                       << '\n';
    }
}

#ifdef AMREX_BACKTRACE_SUPPORTED
namespace {
    std::string run_command (std::string const& cmd)
    {
        std::string r;
        if (FILE * ps = popen(cmd.c_str(), "r")) {
            char print_buff[512];
            while (fgets(print_buff, sizeof(print_buff), ps)) {
                r += print_buff;
            }
            pclose(ps);
        }
        return r;
    }

    int file_exists (const char* file)
    {
        int r = 0;
        if (FILE *fp = fopen(file, "r")) {
            fclose(fp);
            r = 1;
        }
        return r;
    }

#ifdef __linux__
    bool command_exists(std::string const &cmd)
    {
        // command -v is part of POSIX so should be available
        std::string check_command = "command -v " + cmd + " > /dev/null 2>&1";
        int r = std::system(check_command.c_str());
        // return value of std::system is implementation defined and can be
        // decoded using WEXITSTATUS but it should be 0 on success
        return r == 0;
    }

    // Run addr2line once for all addresses. cmd must contain -a so that
    // each entry starts with "0x...: ". Results are keyed by address value.
    std::map<unsigned long long, std::string>
    batch_addr2line (std::string const& cmd, Vector<std::string> const& addrs)
    {
        std::map<unsigned long long, std::string> r;
        if (addrs.empty()) { return r; }
        std::string full_cmd = cmd;
        for (auto const& a : addrs) {
            full_cmd.append(" ").append(a);
        }
        std::istringstream iss(run_command(full_cmd));
        std::string line;
        std::string* cur = nullptr;
        while (std::getline(iss, line)) {
            if (line.starts_with("0x")) {
                auto pos = line.find(": ");
                if (pos != std::string::npos) {
                    cur = &r[std::strtoull(line.c_str(), nullptr, 16)];
                    *cur = line.substr(pos+2) + '\n';
                    continue;
                }
            }
            if (cur) { *cur += line + '\n'; }
        }
        return r;
    }

    // If addr is in the executable, return true and set offset to the address
    // that addr2line -e exe expects. Return false for shared libraries.
    bool exe_offset (std::uintptr_t addr, std::uintptr_t& offset)
    {
        std::uintptr_t a[2] = {addr, 0};
        // The first object dl_iterate_phdr visits is the main program.
        auto f = [] (dl_phdr_info* info, std::size_t, void* data) -> int {
            auto* q = static_cast<std::uintptr_t*>(data);
            for (decltype(info->dlpi_phnum) j = 0; j < info->dlpi_phnum; ++j) {
                auto const& ph = info->dlpi_phdr[j];
                const std::uintptr_t lo = info->dlpi_addr + ph.p_vaddr;
                if (ph.p_type == PT_LOAD && q[0] >= lo && q[0] < lo + ph.p_memsz) {
                    q[1] = q[0] - info->dlpi_addr;
                    return 1;
                }
            }
            return -1;
        };
        if (addr != 0 && dl_iterate_phdr(f, a) == 1) {
            offset = a[1];
            return true;
        }
        return false;
    }

    // Path of the running executable. amrex::system::exename is built from
    // argv[0] and is wrong if the program was found through PATH.
    std::string exe_path ()
    {
        char buf[4096];
        const ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf)-1);
        if (len > 0) {
            buf[len] = '\0';
            return std::string(buf);
        }
        return amrex::system::exename;
    }
#endif
}
#endif

void
BLBackTrace::print_backtrace_info (FILE* f)
{
#ifdef AMREX_BACKTRACE_SUPPORTED

    const int nbuf = 64;
    void *bt_buffer[nbuf];
    int nentries = backtrace(bt_buffer, nbuf);

#ifdef __linux__

    {
        constexpr std::size_t len = 64;
        char host_name[len];
        host_name[len-1] = '\0';
        // The returned buffer may not be null-terminated when truncation occurs.
        // So we insert a null terminator just in case.
        const int ret = gethostname(host_name, len-1);
        if (ret == 0) {
            std::fprintf(f, "Host Name: %s\n", host_name);
        }
    }

    char **strings = backtrace_symbols(bt_buffer, nentries);
    if (strings != nullptr) {
        int have_eu_addr2line = 0;
        int have_addr2line = 0;
        std::string eu_cmd;
        {
            if (command_exists("eu-addr2line")) {
                have_eu_addr2line = 1;
                eu_cmd = "eu-addr2line";
            } else {
                std::string eu_fallback_path = "/usr/bin/eu-addr2line";
                have_eu_addr2line = file_exists(eu_fallback_path.c_str());
                eu_cmd = std::move(eu_fallback_path);
            }
            if (have_eu_addr2line) {
                const pid_t pid = getpid();
                // cmd = "/usr/bin/eu-addr2line -C -f -i --pretty-print -p "
                eu_cmd += " -C -f -i -p " + std::to_string(pid);
            }
        }
        const std::string exe_file = exe_path();
        std::string cmd;
        {
            // GNU addr2line misreads clang's line-tables-only DWARF 5 output.
            // llvm-addr2line handles both gcc and clang output.
            const char* tools[] = {"llvm-addr2line", "addr2line"};
            for (auto const* tool : tools) {
                if (command_exists(tool)) {
                    have_addr2line = 1;
                    cmd = tool;
                    break;
                }
                std::string fallback_path = std::string("/usr/bin/") + tool;
                if (file_exists(fallback_path.c_str())) {
                    have_addr2line = 1;
                    cmd = std::move(fallback_path);
                    break;
                }
            }
            if (have_addr2line) {
                cmd += " -aCpfie '" + exe_file + "'";
            }
        }

        std::fprintf(f, "=== If no file names and line numbers are shown below, one can run\n");
        std::fprintf(f, "            llvm-addr2line -Cpfie my_exefile my_line_address\n");
        std::fprintf(f, "    (or addr2line, which may fail for executables built with clang)\n");
        std::fprintf(f, "    to convert `my_line_address` (e.g., 0x4a6b) into file name and line number.\n");
        std::fprintf(f, "    Or one can use amrex/Tools/Backtrace/parse_bt.py.\n\n");

        std::fprintf(f, "=== Please note that the line number reported by addr2line may not be accurate.\n");
        std::fprintf(f, "    One can use\n");
        std::fprintf(f, "            readelf -wl my_exefile | grep my_line_address'\n");
        std::fprintf(f, "    to find out the offset for that line.\n\n");

        Vector<std::string> results(nentries);

#if !defined(AMREX_USE_OMP) || !defined(__INTEL_COMPILER)
        if (amrex::system::call_addr2line) {
            // A return address can be past the end of the calling function
            // (e.g., after a noreturn call), so look up the address before it.
            auto call_site = [&] (int i) -> std::uintptr_t {
                auto a = reinterpret_cast<std::uintptr_t>(bt_buffer[i]);
                return (a != 0) ? a-1 : 0;
            };
            auto hex = [] (std::uintptr_t a) {
                char print_buff[32];
                std::snprintf(print_buff, sizeof(print_buff), "0x%llx",
                              static_cast<unsigned long long>(a));
                return std::string(print_buff);
            };

            // 0: nothing useful, 1: function only (has '?'), 2: complete
            auto rank = [] (std::string const& s) {
                if (s.empty() || s.starts_with("??")) { return 0; }
                return (s.find('?') == std::string::npos) ? 2 : 1;
            };

            // addr2line -e exe is used only for frames in the executable,
            // looked up by their offset there. Frames in shared libraries
            // are left to eu-addr2line.
            Vector<int> frames;
            Vector<std::string> addrs;
            for (int i = 0; i < nentries; ++i) {
                if (have_eu_addr2line && bt_buffer[i] != nullptr) {
                    results[i] = run_command(eu_cmd + " " + hex(call_site(i)));
                }
                std::uintptr_t offset = 0;
                if (have_addr2line && !exe_file.empty() && rank(results[i]) < 2
                    && exe_offset(call_site(i), offset))
                {
                    frames.push_back(i);
                    addrs.push_back(hex(offset));
                }
            }
            auto res = batch_addr2line(cmd, addrs);
            for (Long k = 0; k < frames.size(); ++k) {
                auto it = res.find(std::strtoull(addrs[k].c_str(), nullptr, 16));
                auto& r = results[frames[k]];
                if (it != res.end() && rank(it->second) > rank(r)) {
                    r = it->second;
                }
            }
        }
#endif

        for (int i = 0; i < nentries; ++i)
        {
            std::fprintf(f, "%2d: %s\n", i, strings[i]);
            if (!results[i].empty()) {
                std::fprintf(f, "    %s", results[i].c_str());
            }
            std::fprintf(f, "\n");
        }
        std::free(strings);
    }

#elif defined(AMREX_BACKTRACE_SUPPORTED) && defined(__APPLE__)

    const pid_t pid = getpid();
    const std::string cmd = "/usr/bin/atos -p " + std::to_string(pid);
    const int have_atos = file_exists("/usr/bin/atos");

    for (int i = 0; i < nentries; ++i) {
        Dl_info info;
        if (bt_buffer[i] != nullptr)
        {
            if (dladdr(bt_buffer[i], &info))
            {
                std::string line;
                if (amrex::system::call_addr2line && have_atos) {
                    char print_buff[32];
                    std::snprintf(print_buff,sizeof(print_buff),"%p",bt_buffer[i]);
                    const std::string full_cmd = cmd + " " + print_buff;
                    line = run_command(full_cmd);
                }
                if (line.empty()) {
                    int status;
                    char * demangled_name = abi::__cxa_demangle(info.dli_sname, nullptr, 0, &status);
                    if (status == 0) {
                        line += demangled_name;
                    } else {
                        line += info.dli_fname;
                    }
                    line += '\n';
                    std::free(demangled_name);
                }
                std::fprintf(f, "%2d: %s\n", i, line.c_str());
            }
        }
    }

#endif

#else
    amrex::ignore_unused(f);
#endif
}

BLBTer::BLBTer(const std::string& s, const char* file, int line)
{
    std::ostringstream ss;
    ss << "Line " << line << ", File " << file;
    line_file = ss.str();

#ifdef AMREX_USE_OMP
    if (omp_in_parallel()) {
        std::ostringstream ss0;
        ss0 << "Proc. " << ParallelDescriptor::MyProc()
            << ", Thread " << omp_get_thread_num()
            << ": \"" << s << "\"";
#ifndef AMREX_BT_STACK_THREADPRIVATE
        #pragma omp critical (amrex_bt_stack)
#endif
        BLBackTrace::bt_stack.emplace(ss0.str(), line_file);
    }
    else {
#ifdef AMREX_BT_STACK_THREADPRIVATE
        #pragma omp parallel
#endif
        {
            std::ostringstream ss0;
            ss0 << "Proc. " << ParallelDescriptor::MyProc()
                << ", Master Thread"
                << ": \"" << s << "\"";
            BLBackTrace::bt_stack.emplace(ss0.str(), line_file);
        }
    }
#else
    std::ostringstream ss0;
    ss0 << "Proc. " << ParallelDescriptor::MyProc()
        << ": \"" << s << "\"";
    BLBackTrace::bt_stack.emplace(ss0.str(), line_file);
#endif
}

BLBTer::~BLBTer()
{
#ifdef AMREX_USE_OMP
    if (omp_in_parallel()) {
#ifndef AMREX_BT_STACK_THREADPRIVATE
        #pragma omp critical (amrex_bt_stack)
#endif
        pop_bt_stack();
    }
    else {
#ifdef AMREX_BT_STACK_THREADPRIVATE
        #pragma omp parallel
#endif
        {
            pop_bt_stack();
        }
    }
#else
    pop_bt_stack();
#endif
}

void
BLBTer::pop_bt_stack()
{
    if (!BLBackTrace::bt_stack.empty()) {
        if (BLBackTrace::bt_stack.top().second == line_file) {
            BLBackTrace::bt_stack.pop();
        }
    }
}

}
