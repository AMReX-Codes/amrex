#include <AMReX_CrtReport.H>
#include <AMReX.H>

#if defined(_MSC_VER) && defined(_DEBUG)
#include <AMReX_ParallelDescriptor.H>

#include <windows.h>
#include <crtdbg.h>

#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cwchar>

namespace {
    bool hook_installed = false;
    bool hook_w_installed = false;
    bool abort_behavior_changed = false;
    unsigned int prev_abort_behavior = 0;
    int rank = 0; // cached: reports can come after EndParallel

    /** Errors and assertions that the CRT would show in a dialog box */
    bool is_handled (int report_type)
    {
        if (report_type == _CRT_WARN) { return false; }
        int const mode = _CrtSetReportMode(report_type, _CRTDBG_REPORT_MODE);
        return mode != -1 && (mode & _CRTDBG_MODE_WNDW) != 0; // -1: unknown report type
    }

    /** A CRT report message with the rank, written to stderr at once
     *
     * This does not allocate on the heap, because the debug heap itself
     * reports, e.g., heap corruption, and avoids AMReX and CRT I/O functions,
     * which could fail and report again (the rank is cached). A single write
     * keeps reports from concurrent threads apart.
     */
    class Message
    {
        static constexpr std::size_t m_size = 4096;
        static constexpr std::size_t m_capacity = m_size - 16; // room for the end

    public:
        //! UTF-16 code units that fit at most (UTF-8 needs at most 3 bytes per code unit)
        static constexpr std::size_t max_utf16 = m_capacity / 3;

        explicit Message (int report_type)
        {
            append("\n");
            char rank_str[16];
            auto const r = std::to_chars(rank_str, rank_str + sizeof(rank_str), rank);
            append(rank_str, static_cast<std::size_t>(r.ptr - rank_str));
            append(report_type == _CRT_ASSERT ? "::MSVC debug CRT assertion: "
                                              : "::MSVC debug CRT error: ");
        }

        /** Append UTF-8 text (or ASCII), truncated only at a code point boundary */
        void append (char const* str, std::size_t len)
        {
            if (m_truncated) { return; }
            if (len > m_capacity - m_len) {
                len = m_capacity - m_len;
                m_truncated = true;
                // drop a code point that does not fit completely
                while (len > 0 && (static_cast<unsigned char>(str[len]) & 0xC0u) == 0x80u) { --len; }
            }
            std::memcpy(m_buf + m_len, str, len);
            m_len += len;
        }

        void append (char const* str) { append(str, std::strlen(str)); }

        void set_truncated () { m_truncated = true; }

        /** Append UTF-16 text as UTF-8, with U+FFFD for unpaired surrogates */
        void append (wchar_t const* str, std::size_t len)
        {
            if (m_truncated) { return; }
            std::size_t const max_len = (m_capacity - m_len) / 3;
            if (len > max_len) {
                len = max_len;
                m_truncated = true;
                auto const last = (len > 0) ? static_cast<std::uint16_t>(str[len-1]) : 0u;
                if (last >= 0xD800u && last < 0xDC00u) { --len; } // keep surrogate pairs
            }
            if (len == 0) { return; }
            int const n = WideCharToMultiByte(CP_UTF8, 0, str, static_cast<int>(len),
                                              m_buf + m_len, static_cast<int>(m_capacity - m_len),
                                              nullptr, nullptr);
            if (n > 0) { m_len += static_cast<std::size_t>(n); }
        }

        /** Write the message and request a break into the debugger
         *
         * An attached debugger stops at the failure. Without a debugger, the
         * system's post-mortem handling applies: by default, the process
         * terminates with STATUS_BREAKPOINT and Windows Error Reporting can
         * write a crash dump; a registered just-in-time debugger can start
         * instead.
         */
        int write_and_break (int* return_value)
        {
            char const* end = m_truncated ? " ... !!!\n" : " !!!\n";
            std::size_t const end_len = std::strlen(end);
            std::memcpy(m_buf + m_len, end, end_len); // fits: m_capacity < sizeof(m_buf)

            HANDLE const err = GetStdHandle(STD_ERROR_HANDLE);
            if (err != nullptr && err != INVALID_HANDLE_VALUE) {
                DWORD written = 0;
                WriteFile(err, m_buf, static_cast<DWORD>(m_len + end_len), &written, nullptr);
            }
            *return_value = 1; // the reporting function breaks into the debugger
            return 1; // report handled: no further processing by the CRT
        }

    private:
        char m_buf[m_size];
        std::size_t m_len = 0;
        bool m_truncated = false;
    };

    int __cdecl report_hook (int report_type, char* message, int* return_value)
    {
        if (!is_handled(report_type)) { return 0; }

        Message msg(report_type);
        std::size_t len = (message != nullptr) ? std::strlen(message) : 0;
        while (len > 0 && (message[len-1] == '\n' || message[len-1] == '\r')) { --len; }

        if (len == 0) { return msg.write_and_break(return_value); }

        // valid UTF-8 (incl. ASCII, e.g., __FILE__ with /utf-8): as is
        if (MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, message, static_cast<int>(len),
                                nullptr, 0) > 0)
        {
            msg.append(message, len);
            return msg.write_and_break(return_value);
        }

        // else the ANSI code page to UTF-16 (at most one code unit per byte), then UTF-8
        wchar_t wide[Message::max_utf16];
        std::size_t const wide_max = sizeof(wide) / sizeof(wide[0]);
        bool const truncated = len > wide_max;
        if (truncated) { len = wide_max; }
        int const n = MultiByteToWideChar(CP_ACP, 0, message, static_cast<int>(len),
                                          wide, static_cast<int>(wide_max));
        if (n > 0) { msg.append(wide, static_cast<std::size_t>(n)); }
        if (truncated) { msg.set_truncated(); }
        return msg.write_and_break(return_value);
    }

    int __cdecl report_hook_w (int report_type, wchar_t* message, int* return_value)
    {
        if (!is_handled(report_type)) { return 0; }

        Message msg(report_type);
        std::size_t len = (message != nullptr) ? std::wcslen(message) : 0;
        while (len > 0 && (message[len-1] == L'\n' || message[len-1] == L'\r')) { --len; }
        if (len > 0) { msg.append(message, len); }
        return msg.write_and_break(return_value);
    }
}
#endif

namespace amrex::detail
{

void
CrtReportInitialize ()
{
#if defined(_MSC_VER) && defined(_DEBUG)
    rank = amrex::ParallelDescriptor::MyProc();
    if (!hook_installed) {
        hook_installed = _CrtSetReportHook2(_CRT_RPTHOOK_INSTALL, report_hook) != -1;
    }
    if (!hook_w_installed) {
        hook_w_installed = _CrtSetReportHookW2(_CRT_RPTHOOK_INSTALL, report_hook_w) != -1;
    }
    if ((!hook_installed || !hook_w_installed) && amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Warning("amrex.handle_crt_reports: could not install CRT report hooks");
    }
#endif
}

void
CrtReportSkipAbortReport (bool skip)
{
#if defined(_MSC_VER) && defined(_DEBUG)
    if (skip && !abort_behavior_changed) {
        prev_abort_behavior = _set_abort_behavior(0, _WRITE_ABORT_MSG);
        abort_behavior_changed = true;
    } else if (!skip && abort_behavior_changed) {
        _set_abort_behavior(prev_abort_behavior, _WRITE_ABORT_MSG);
        abort_behavior_changed = false;
    }
#else
    amrex::ignore_unused(skip);
#endif
}

void
CrtReportFinalize ()
{
#if defined(_MSC_VER) && defined(_DEBUG)
    if (hook_installed) {
        _CrtSetReportHook2(_CRT_RPTHOOK_REMOVE, report_hook);
        hook_installed = false;
    }
    if (hook_w_installed) {
        _CrtSetReportHookW2(_CRT_RPTHOOK_REMOVE, report_hook_w);
        hook_w_installed = false;
    }
#endif
    CrtReportSkipAbortReport(false); // if not yet restored with the SIGABRT handler
}

}
