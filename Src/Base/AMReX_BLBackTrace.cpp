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

#if defined(__APPLE__)
#include <cxxabi.h>
#include <dlfcn.h>
#include <mach-o/loader.h>
#include <pthread.h>
#include <sys/ucontext.h>
#define AMREX_BACKTRACE_SUPPORTED 1
#elif defined(__linux__) && defined(__GLIBC__)
#include <link.h>
#include <sys/ucontext.h>
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
// info and uctx are the signal information and context in the signal
// handler, otherwise nullptr. in_handler is true for the Backtrace file of
// the signal handler.
void print_backtrace_info_impl (FILE* f, void* info, void* uctx, bool in_handler);

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
handle_signal (int s, void* info, void* uctx)
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
        print_backtrace_info_impl(p, info, uctx, true);
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

#else
    amrex::ignore_unused(info, uctx);
#endif

    ParallelDescriptor::Abort(s, false);
}

#ifdef _WIN32
void signal_handler (int s)
{
    handle_signal(s, nullptr, nullptr);
}
#else
void signal_handler (int s, siginfo_t* info, void* uctx)
{
    handle_signal(s, info, uctx);
}
#endif
}

void
BLBackTrace::handler (int s)
{
    handle_signal(s, nullptr, nullptr);
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

    std::string hex (std::uintptr_t a)
    {
        char print_buff[32];
        std::snprintf(print_buff, sizeof(print_buff), "0x%llx",
                      static_cast<unsigned long long>(a));
        return std::string(print_buff);
    }

    // A return address can be past the end of the calling function (e.g.,
    // after a noreturn call), so look up the address before it.
    std::uintptr_t call_site (void const* return_address)
    {
        auto a = reinterpret_cast<std::uintptr_t>(return_address);
        return (a != 0) ? a-1 : 0;
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

    // In the signal handler, the frames are written with print_frame(i) before
    // they are looked up, which can take long, in case the process ends
    // meanwhile. The final frames, which are not shorter, replace them.
    template <typename F>
    void write_frames_early (FILE* f, int nframes, F const& print_frame)
    {
        const long pos = std::ftell(f);
        for (int i = 0; i < nframes; ++i) {
            print_frame(i);
        }
        std::fflush(f);
        if (pos >= 0) { std::fseek(f, pos, SEEK_SET); }
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

#if defined(__APPLE__)
    // Start of the function that contains a in the image at header, from the
    // function starts that the linker records (LC_FUNCTION_STARTS), which
    // also stripped images have. Returns 0 if they are not found.
    std::uintptr_t function_start (void const* header, std::uintptr_t a)
    {
        auto const* mh = static_cast<mach_header_64 const*>(header);
        if (mh == nullptr || mh->magic != MH_MAGIC_64) { return 0; }
        const auto base = reinterpret_cast<std::uintptr_t>(mh);
        segment_command_64 const* text = nullptr;
        segment_command_64 const* linkedit = nullptr;
        linkedit_data_command const* starts = nullptr;
        auto cmd = base + sizeof(mach_header_64);
        for (std::uint32_t i = 0; i < mh->ncmds; ++i) {
            auto const* lc = reinterpret_cast<load_command const*>(cmd); // NOLINT(performance-no-int-to-ptr)
            if (lc->cmdsize < sizeof(load_command)) { return 0; }
            if (lc->cmd == LC_SEGMENT_64) {
                auto const* seg = reinterpret_cast<segment_command_64 const*>(lc);
                if (std::strncmp(seg->segname, SEG_TEXT, sizeof(seg->segname)) == 0) { text = seg; }
                if (std::strncmp(seg->segname, SEG_LINKEDIT, sizeof(seg->segname)) == 0) { linkedit = seg; }
            } else if (lc->cmd == LC_FUNCTION_STARTS) {
                starts = reinterpret_cast<linkedit_data_command const*>(lc);
            }
            cmd += lc->cmdsize;
        }
        if (text == nullptr || linkedit == nullptr || starts == nullptr) { return 0; }
        // ULEB128 offsets, the first from the start of __TEXT, which is base
        const std::uintptr_t slide = base - text->vmaddr;
        auto const* p = reinterpret_cast<unsigned char const*>( // NOLINT(performance-no-int-to-ptr)
            slide + linkedit->vmaddr - linkedit->fileoff + starts->dataoff);
        auto const* end = p + starts->datasize;
        std::uintptr_t addr = base;
        std::uintptr_t found = 0;
        while (p < end) {
            std::uintptr_t delta = 0;
            for (unsigned shift = 0; p < end && shift < 64; shift += 7) {
                delta |= static_cast<std::uintptr_t>(*p & 0x7fU) << shift;
                if ((*p++ & 0x80U) == 0) { break; }
            }
            if (delta == 0) { break; }
            addr += delta;
            if (addr > a) { break; }
            found = addr;
        }
        return found;
    }

    // Walks the stack of the code that a signal interrupted, starting from
    // the signal information info and the registers in the signal context
    // uctx: frames[0] is the PC that raised the signal, the others are
    // return addresses. Returns the number of frames, or 0 without a signal
    // context.
    //
    // backtrace() called in a signal handler misses the interrupted PC. Code
    // on macOS keeps a chain of frame records (the caller's fp and the return
    // address), which this follows from the interrupted fp. A function that
    // has not stored its return address in its frame record (e.g., a leaf, or
    // in its prologue or epilogue) still has it in lr on arm64, and on top of
    // the stack on x86_64 at its first instruction and at its ret.
    int unwind_signal_context (void* info, void* uctx, void** frames, int max_frames)
    {
#if defined(__aarch64__) || defined(__x86_64__)
        auto const* si = static_cast<siginfo_t const*>(info);
        auto const* uc = static_cast<ucontext_t const*>(uctx);
        if (uc == nullptr || uc->uc_mcontext == nullptr || max_frames < 2) { return 0; }
        auto const& ss = uc->uc_mcontext->__ss;

        // Do not follow a corrupted stack outside of this thread's stack: a
        // fault in this signal handler would not be handled.
        pthread_t self = pthread_self();
        const auto stack_hi = reinterpret_cast<std::uintptr_t>(pthread_get_stackaddr_np(self));
        const std::uintptr_t stack_lo = stack_hi - pthread_get_stacksize_np(self);
        auto on_stack = [=] (std::uintptr_t a) {
            return a % sizeof(void*) == 0 && a >= stack_lo && a < stack_hi - 2*sizeof(void*);
        };
        auto load = [] (std::uintptr_t a) {
            return *reinterpret_cast<std::uintptr_t const*>(a); // NOLINT(performance-no-int-to-ptr)
        };
        auto in_image = [] (std::uintptr_t a) {
            Dl_info dl;
            return dladdr(reinterpret_cast<void*>(a), &dl) != 0; // NOLINT(performance-no-int-to-ptr)
        };
        // Whether a is in a loaded image, and the start of its function, or 0
        // if unknown.
        struct Function { bool in_image = false; std::uintptr_t start = 0; };
        auto function_of = [] (std::uintptr_t a) {
            Function r;
            Dl_info dl;
            if (dladdr(reinterpret_cast<void*>(a), &dl) != 0) { // NOLINT(performance-no-int-to-ptr)
                r.in_image = true;
                r.start = function_start(dl.dli_fbase, a);
                if (r.start == 0) { r.start = reinterpret_cast<std::uintptr_t>(dl.dli_saddr); }
            }
            return r;
        };

        // ra is the return address of the interrupted function if it has not
        // stored it in its frame record (e.g., a leaf, or in its prologue).
        // Without the start of pc's function, this cannot be told, so ra is
        // then only used if pc is not in an image (e.g., after a call through
        // a bad pointer). It must follow a call instruction, which another
        // value in lr (used as a scratch register) or on the stack does not.
        std::uintptr_t ra = 0;
#if defined(__aarch64__)
        amrex::ignore_unused(si);
        const std::uintptr_t pc = __darwin_arm_thread_state64_get_pc(ss);
        const std::uintptr_t fp = __darwin_arm_thread_state64_get_fp(ss);
        const std::uintptr_t lr = __darwin_arm_thread_state64_get_lr(ss);
        const Function pc_function = function_of(pc);
        auto insn_at = [] (std::uintptr_t a) {
            return *reinterpret_cast<std::uint32_t const*>(a); // NOLINT(performance-no-int-to-ptr)
        };
        // Target of the branch (b or bl) at a
        auto branch_target = [] (std::uintptr_t a, std::uint32_t insn) {
            const auto offset = static_cast<std::intptr_t>(static_cast<std::int32_t>(insn << 6) >> 6) * 4;
            return a + static_cast<std::uintptr_t>(offset);
        };
        // Whether a follows a call instruction: bl, or blr (also with pointer
        // authentication)
        auto follows_call = [&] (std::uintptr_t a) {
            if (a % 4 != 0 || !in_image(a-4)) { return false; }
            const auto insn = insn_at(a-4);
            return (insn & 0xfc000000U) == 0x94000000U || (insn & 0xfefff000U) == 0xd63f0000U;
        };
        // Whether the interrupted function has set up its frame record: the
        // code before pc sets fp (mov or add x29, sp), and the code from pc on
        // restores the frame record (ldp x29, x30) before it returns.
        auto frame_set_up = [&] () {
            bool sets_fp = false;
            for (std::uintptr_t a = pc_function.start; a < pc && a < pc_function.start + 65536; a += 4) {
                if ((insn_at(a) & 0xff8003ffU) == 0x910003fdU) {
                    sets_fp = true;
                    break;
                }
            }
            if (!sets_fp) { return false; }
            std::uintptr_t a = pc;
            for (int i = 0; i < 64; ++i) {
                const auto insn = insn_at(a);
                if ((insn & 0xffc07fffU) == 0xa8c07bfdU || (insn & 0xffc07fffU) == 0xa9407bfdU) {
                    return true;
                }
                if (insn == 0xd65f03c0U || insn == 0xd65f0bffU || insn == 0xd65f0fffU) { // ret, retaa, retab
                    return false;
                }
                a = ((insn & 0xfc000000U) == 0x14000000U) ? branch_target(a, insn) : a+4; // b
                if ((a % 4096 == 0 || (insn & 0xfc000000U) == 0x14000000U) && !in_image(a)) {
                    break;
                }
            }
            return true;
        };
        if (lr != 0 && lr != pc && follows_call(lr)) {
            const Function lr_function = function_of(lr-1);
            if (!pc_function.in_image) {
                ra = lr;
            } else if (pc_function.start != 0 && lr_function.in_image) {
                if (lr_function.start != pc_function.start) {
                    // unless the interrupted function stored lr in its frame record
                    if (!(on_stack(fp) && load(fp + sizeof(std::uintptr_t)) == lr)) { ra = lr; }
                } else {
                    // lr in the same function is from an earlier call, unless
                    // it returns from a recursive call (bl to this function)
                    // and this function has not set up its frame record.
                    const auto insn = insn_at(lr-4);
                    if ((insn & 0xfc000000U) == 0x94000000U && // bl
                        branch_target(lr-4, insn) == pc_function.start && !frame_set_up()) {
                        ra = lr;
                    }
                }
            }
        }
#else
        auto opcode = [] (std::uintptr_t a) {
            return *reinterpret_cast<unsigned char const*>(a); // NOLINT(performance-no-int-to-ptr)
        };
        std::uintptr_t pc = ss.__rip;
        // After an int3 instruction (e.g., __builtin_debugtrap), which raises
        // SIGTRAP, pc is the next instruction, which can be in the next
        // function. macOS does not set si_code to TRAP_BRKPT.
        if (si != nullptr && si->si_signo == SIGTRAP && si->si_code != SI_USER &&
            si->si_code != SI_QUEUE && in_image(pc-1) && opcode(pc-1) == 0xcc) {
            --pc;
        }
        const std::uintptr_t fp = ss.__rbp;
        const Function pc_function = function_of(pc);
        // Whether a follows a call instruction: call rel32 (E8), or call
        // through a register or memory (FF /2), 2 to 7 bytes long
        auto follows_call = [&] (std::uintptr_t a) {
            if (!in_image(a-7) || !in_image(a-1)) { return false; }
            if (opcode(a-5) == 0xe8) { return true; }
            for (std::uintptr_t k : {2, 3, 4, 6, 7}) {
                if (opcode(a-k) == 0xff && ((opcode(a-k+1) >> 3) & 7) == 2) { return true; }
            }
            return false;
        };
        // The return address is on top of the stack at the first instruction,
        // after push rbp, and at ret.
        if (on_stack(ss.__rsp)) {
            // ret, ret imm16, rep ret
            const bool at_ret = pc_function.in_image &&
                (opcode(pc) == 0xc3 || opcode(pc) == 0xc2 || (opcode(pc) == 0xf3 && opcode(pc+1) == 0xc3));
            if (!pc_function.in_image) {
                ra = load(ss.__rsp);
            } else if (pc_function.start != 0) {
                if (pc == pc_function.start || at_ret) {
                    ra = load(ss.__rsp);
                } else if (pc == pc_function.start+1 && opcode(pc_function.start) == 0x55) {
                    ra = load(ss.__rsp + sizeof(std::uintptr_t));
                }
            }
        }
        if (ra != 0 && !follows_call(ra)) { ra = 0; }
#endif

        int n = 0;
        frames[n++] = reinterpret_cast<void*>(pc); // NOLINT(performance-no-int-to-ptr)
        if (ra != 0 && ra != pc) {
            frames[n++] = reinterpret_cast<void*>(ra); // NOLINT(performance-no-int-to-ptr)
        }
        // backtrace_from_fp cannot follow the chain here: it starts from the
        // current frame, which is on the alternate stack after an overflow.
        for (std::uintptr_t record = fp; n < max_frames && on_stack(record); ) {
            const std::uintptr_t ret = load(record + sizeof(std::uintptr_t));
            if (ret == 0) { break; }
            frames[n++] = reinterpret_cast<void*>(ret); // NOLINT(performance-no-int-to-ptr)
            const std::uintptr_t next = load(record);
            if (next <= record) { break; } // the stack grows down
            record = next;
        }
        return n;
#else
        amrex::ignore_unused(info, uctx, frames, max_frames);
        return 0;
#endif
    }
#endif
}
#endif

void
BLBackTrace::print_backtrace_info (FILE* f)
{
    print_backtrace_info_impl(f, nullptr, nullptr, false);
}

namespace {
void
print_backtrace_info_impl (FILE* f, void* info, void* uctx, bool in_handler)
{
#ifdef AMREX_BACKTRACE_SUPPORTED

    const int nbuf = 64;
    void *bt_buffer[nbuf];
    int nentries = 0;

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

#ifdef __linux__

    amrex::ignore_unused(info, uctx);
    nentries = backtrace(bt_buffer, nbuf);

    char **strings = backtrace_symbols(bt_buffer, nentries);
    if (strings != nullptr) {
        std::fprintf(f, "=== If no file names and line numbers are shown below, one can run\n");
        std::fprintf(f, "            llvm-addr2line -Cpfie my_exefile my_line_address\n");
        std::fprintf(f, "    (or addr2line, which may fail for executables built with clang)\n");
        std::fprintf(f, "    to convert `my_line_address` (e.g., 0x4a6b) into file name and line number.\n");
        std::fprintf(f, "    Or one can use amrex/Tools/Backtrace/parse_bt.py.\n\n");

        std::fprintf(f, "=== Please note that the line number reported by addr2line may not be accurate.\n");
        std::fprintf(f, "    One can use\n");
        std::fprintf(f, "            readelf -wl my_exefile | grep my_line_address'\n");
        std::fprintf(f, "    to find out the offset for that line.\n\n");

        if (in_handler) {
            write_frames_early(f, nentries, [&] (int i) {
                std::fprintf(f, "%2d: %s\n\n", i, strings[i]);
            });
        }

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

        Vector<std::string> results(nentries);

#if !defined(AMREX_USE_OMP) || !defined(__INTEL_COMPILER)
        if (amrex::system::call_addr2line) {
            // 0: nothing useful, 1: function only (has '?'), 2: complete
            auto rank = [] (std::string const& s) {
                if (s.empty() || s.starts_with("??")) { return 0; }
                return (s.find('?') == std::string::npos) ? 2 : 1;
            };

            // backtrace() returns the PC that raised the signal as it is, which
            // is looked up as it is (on x86_64, at the int3 instruction before
            // it if that raised SIGTRAP), and return addresses, which are
            // looked up at their call instruction.
            std::uintptr_t fault_pc = 0;
            std::uintptr_t fault_addr = 0;
            if (auto const* uc = static_cast<ucontext_t const*>(uctx)) {
#if defined(__x86_64__) && defined(REG_RIP)
                fault_pc = static_cast<std::uintptr_t>(uc->uc_mcontext.gregs[REG_RIP]);
                fault_addr = fault_pc;
                auto const* si = static_cast<siginfo_t const*>(info);
                if (si != nullptr && si->si_signo == SIGTRAP && si->si_code == SI_KERNEL &&
                    *reinterpret_cast<unsigned char const*>(fault_pc-1) == 0xcc) // NOLINT(performance-no-int-to-ptr)
                {
                    fault_addr = fault_pc-1;
                }
#elif defined(__aarch64__)
                fault_pc = static_cast<std::uintptr_t>(uc->uc_mcontext.pc);
                fault_addr = fault_pc;
#else
                amrex::ignore_unused(uc);
#endif
            }
            // the first frame with it, after the signal frame
            int fault_frame = -1;
            for (int i = 0; i < nentries && fault_pc != 0; ++i) {
                if (reinterpret_cast<std::uintptr_t>(bt_buffer[i]) == fault_pc) {
                    fault_frame = i;
                    break;
                }
            }
            auto lookup_addr = [&] (int i) {
                return (i == fault_frame) ? fault_addr : call_site(bt_buffer[i]);
            };

            // addr2line -e exe is used only for frames in the executable,
            // looked up by their offset there. Frames in shared libraries
            // are left to eu-addr2line.
            Vector<int> frames;
            Vector<std::string> addrs;
            for (int i = 0; i < nentries; ++i) {
                if (have_eu_addr2line && bt_buffer[i] != nullptr) {
                    results[i] = run_command(eu_cmd + " " + hex(lookup_addr(i)));
                }
                std::uintptr_t offset = 0;
                if (have_addr2line && !exe_file.empty() && rank(results[i]) < 2
                    && exe_offset(lookup_addr(i), offset))
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

    // Frames to print: the stack of the interrupted code if we are in a
    // signal handler, otherwise the current stack (e.g., amrex::Abort).
    void* sig_frames[nbuf];
    const int n_sig_frames = unwind_signal_context(info, uctx, sig_frames, nbuf);
    if (n_sig_frames == 0) {
        nentries = backtrace(bt_buffer, nbuf);
    }
    const bool has_fault_frame = n_sig_frames > 0;
    void** frames = has_fault_frame ? sig_frames : bt_buffer;
    const int nframes = has_fault_frame ? n_sig_frames : nentries;

    // The address to look up for each frame: the PC that raised the signal,
    // or the call instruction of a return address, and dladdr's result for it
    Vector<std::uintptr_t> addrs(nframes);
    Vector<Dl_info> dl(nframes);
    for (int i = 0; i < nframes; ++i) {
        addrs[i] = (i == 0 && has_fault_frame) ? reinterpret_cast<std::uintptr_t>(frames[i])
                                               : call_site(frames[i]);
        if (dladdr(reinterpret_cast<void*>(addrs[i]), &dl[i]) == 0) { dl[i] = Dl_info{}; } // NOLINT(performance-no-int-to-ptr)
    }
    // Writes the image, load address and address of frame i
    auto print_frame = [&] (int i) {
        std::fprintf(f, "%2d: [%s 0x%llx 0x%llx]\n", i,
                     dl[i].dli_fname ? dl[i].dli_fname : "??",
                     static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(dl[i].dli_fbase)),
                     static_cast<unsigned long long>(addrs[i]));
    };

    std::fprintf(f, "=== If no file names and line numbers are shown below, one can run\n");
    std::fprintf(f, "            atos -i -o my_image -l my_load_address my_line_address\n");
    std::fprintf(f, "    using the image, load address and address shown in [] for each frame.\n");
    std::fprintf(f, "    my_image must have debug info: either its object files still exist\n");
    std::fprintf(f, "    or a dSYM was made with `dsymutil my_image` (then pass my_image.dSYM).\n\n");

    if (in_handler) {
        write_frames_early(f, nframes, print_frame);
    }

    Vector<std::string> results(nframes);
    if (amrex::system::call_addr2line && file_exists("/usr/bin/atos")) {
        // One atos call for all frames. With -i, the entries of each address
        // (inlined frames first) are followed by a blank line. An empty
        // result is an empty line, followed by the blank line.
        std::string cmd = "/usr/bin/atos -i -p " + std::to_string(getpid());
        for (int i = 0; i < nframes; ++i) {
            cmd += " " + hex(addrs[i]);
        }
        std::istringstream iss(run_command(cmd + " </dev/null 2>/dev/null"));
        Vector<std::string> groups;
        std::string line, cur;
        bool empty_result = false;
        while (std::getline(iss, line)) {
            if (!line.empty()) {
                if (!cur.empty()) { cur += "    "; }
                cur += line + '\n';
            } else if (!cur.empty() || empty_result) {
                groups.push_back(cur);
                cur.clear();
                empty_result = false;
            } else {
                empty_result = true;
            }
        }
        if (!cur.empty() || empty_result) { groups.push_back(cur); }
        if (static_cast<int>(groups.size()) == nframes) {
            for (int i = 0; i < nframes; ++i) { results[i] = std::move(groups[i]); }
        }
    }

    for (int i = 0; i < nframes; ++i) {
        print_frame(i);
        std::string line = results[i];
        if (line.empty() && dl[i].dli_fbase != nullptr) {
            char const* sname = dl[i].dli_sname;
            int status = -1;
            char* demangled_name = sname
                ? abi::__cxa_demangle(sname, nullptr, 0, &status) : nullptr;
            line = (status == 0) ? demangled_name : (sname ? sname : "??");
            line += '\n';
            std::free(demangled_name);
        }
        std::fprintf(f, "    %s\n", line.c_str());
    }

#endif

#else
    amrex::ignore_unused(f, info, uctx, in_handler);
#endif
}
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
