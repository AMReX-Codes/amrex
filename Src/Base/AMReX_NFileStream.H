// ABOUTME: This file defines AMReX_NFileStream, a custom stream class for NFilesIter
// ABOUTME: It uses POSIX I/O directly with proper error handling and exponential backoff

#ifndef AMREX_NFILESTREAM_H_
#define AMREX_NFILESTREAM_H_

#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_VisMFBuffer.H>
#include <AMReX_GpuQualifiers.H>
#include <string>
#include <ios>
#include <sstream>

namespace amrex {

/**
 * \brief A custom file stream class for NFilesIter that uses POSIX I/O directly
 *
 * This class provides the same interface as std::fstream but uses POSIX I/O calls
 * with proper error handling and exponential backoff for retryable errors.
 */
class NFileStream
{
public:
    NFileStream() = default;
    ~NFileStream();

    // File operations
    void open(const std::string& filename, std::ios_base::openmode mode);
    void close();
    [[nodiscard]] bool is_open() const { return m_fd >= 0; }
    [[nodiscard]] bool good() const { return m_good && is_open(); }

    // Stream state
    [[nodiscard]] bool fail() const { return m_fail; }
    [[nodiscard]] bool eof() const { return m_eof; }
    void clear() { m_fail = false; m_eof = false; m_good = true; }

    // Write operations
    NFileStream& write(const char* s, std::streamsize n);
    NFileStream& flush();

    // Read operations
    NFileStream& read(char* s, std::streamsize n);

    // Positioning
    NFileStream& seekp(std::streampos pos);
    NFileStream& seekp(std::streamoff off, std::ios_base::seekdir way);
    [[nodiscard]] std::streampos tellp();

    // Buffer management - for compatibility with VisMFBuffer
    std::streambuf* rdbuf() const { return m_stringbuf.rdbuf(); }

    // Conversion operators for compatibility
    explicit operator bool() const { return !fail(); }

    // Allow implicit conversion to std::istream& for reading operations
    operator std::istream&() { return m_stringbuf; }

    // Allow implicit conversion to std::ostream& for writing operations
    operator std::ostream&() { return m_stringbuf; }

private:
    int m_fd = -1;
    bool m_good = true;
    bool m_fail = false;
    bool m_eof = false;
    std::ios_base::openmode m_mode = {};
    mutable std::stringstream m_stringbuf;  // For rdbuf() compatibility

    // Helper methods
    int convert_openmode_to_flags(std::ios_base::openmode mode);
    void set_error_state(bool fail_state, bool eof_state = false);
};

} // namespace amrex

#endif // AMREX_NFILESTREAM_H_
