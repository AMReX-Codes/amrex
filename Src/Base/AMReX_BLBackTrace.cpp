#include <AMReX_BLBackTrace.H>
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
#include <unistd.h>
#endif

namespace amrex {

std::stack<std::pair<std::string, std::string> >  BLBackTrace::bt_stack;

void
BLBackTrace::handler(int s)
{

    std::signal(s, SIG_DFL);

    AsyncOut::Finalize();

    switch (s) {
    case SIGSEGV:
        amrex::ErrorStream() << "Segfault\n";
        break;
    case SIGFPE:
        amrex::ErrorStream() << "Erroneous arithmetic operation\n";
        break;
    case SIGILL:
        amrex::ErrorStream() << "SIGILL Invalid, privileged, or ill-formed instruction\n";
        break;
    case SIGTERM:
        amrex::ErrorStream() << "SIGTERM\n";
        break;
    case SIGINT:
        amrex::ErrorStream() << "SIGINT\n";
        break;
    case SIGABRT:
        amrex::ErrorStream() << "SIGABRT\n";
        break;
    default: break;
    }

#if defined(AMREX_BACKTRACE_SUPPORTED) || defined(AMREX_TINY_PROFILING)

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
