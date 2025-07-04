#ifndef AMREX_FILE_STREAM_H_
#define AMREX_FILE_STREAM_H_
#include <AMReX_Config.H>

#include <AMReX_INT.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#ifndef _WIN32

#include <cerrno>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <sys/types.h>

namespace amrex
{

class FileStream
{
public:

    using off_type = std::char_traits<char>::off_type;
    using pos_type = std::char_traits<char>::pos_type;

    explicit FileStream (char const* filename,
                         std::ios_base::openmode mode
                             = std::ios_base::in | std::ios_base::out);
    explicit FileStream (std::string const& filename,
                         std::ios_base::openmode mode
                             = std::ios_base::in | std::ios_base::out);

    FileStream () = default;
    FileStream (FileStream&& rhs) noexcept;
    FileStream (FileStream const& rhs) = delete;

    ~FileStream ();

    FileStream& operator= (FileStream&& rhs) noexcept;
    FileStream& operator= (FileStream const& rhs) = delete;

    void open (char const* filename,
               std::ios_base::openmode mode
                   = std::ios_base::in | std::ios_base::out);
    void open (std::string const& filename,
               std::ios_base::openmode mode
                   = std::ios_base::in | std::ios_base::out);

    void close ();

    FileStream& read (char* s, std::streamsize count);

    FileStream& write (char const* s, std::streamsize count);

    FileStream& flush ();

    FileStream& seekp (pos_type off);
    FileStream& seekp (off_type off, std::ios_base::seekdir dir);

    FileStream& seekg (pos_type off);
    FileStream& seekg (off_type off, std::ios_base::seekdir dir);

    [[nodiscard]] pos_type tellp () const;

    [[nodiscard]] pos_type tellg () const;

    [[nodiscard]] bool good () const;

    [[nodiscard]] bool bad () const;

    [[nodiscard]] bool fail () const;

private:

    double m_base_delay = 0.1;
    off_t m_pos = 0;
    off_t m_size = 0;
    int m_max_retries = 16;
    int m_fd = -1;
    bool m_good = false;
    bool m_binary = true; // not used
    bool m_append = false;

    std::unique_ptr<char[]> m_buffer;
    Long m_buffer_size = 32*1024*1024; // 32MB
    Long m_buffer_begin = 0;
    Long m_buffer_end = 0;
    enum class BufferMode { Read, Write, None};
    BufferMode m_buffer_mode = BufferMode::None;

    template <typename F>
    void execute_with_retry (F const& func, std::string const& op_name)
    {
        int attempt = 0;
        double delay = m_base_delay;

        while (attempt <= m_max_retries) {
            if (func()) { return; }

            // Retry for certain error
            if (errno == EINTR             // Interrupted function
                || errno == EAGAIN         // Resource unavailable, try again
                || errno == EWOULDBLOCK    // Operation would block
                || errno == EMFILE         // Too many files this process
                || errno == ENFILE         // Too many files sysetem wide
                || errno == EBUSY          // Busy
                || errno == ENOENT         // No such file or directory
                || errno == ENOBUFS        // No buffer
                || errno == ENODEV         // No such device
                || errno == ENOMEM         // Out of memory
                || errno == ENOTCONN       // Socket not connected
                || errno == ESTALE         // Stale file handle
                || errno == ETIMEDOUT      // Network timeout
                || errno == EIO)           // I/O error
            {
                attempt++;
                amrex::Sleep(delay);
                delay *= 2; // Exponential backoff
                if (amrex::Verbose()) {
                    amrex::AllPrint() << op_name << " failed after "
                                      << std::to_string(attempt)
                                      << ": " << std::strerror(errno) << "\n";
                }
            } else {
                break;
            }
        }

        m_good = false;
        throw std::runtime_error(op_name + " failed: " + std::strerror(errno));
    }

    void file_write (char const* s, Long count);

    void flush_write_buffer ();

    void fill_write_buffer (char const* s, Long count);

    void fill_read_buffer ();
};

}

#endif

#endif
