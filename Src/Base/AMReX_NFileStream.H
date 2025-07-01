// ABOUTME: This file defines AMReX_NFileStream, a custom stream class for NFilesIter
// ABOUTME: It replaces std::fstream with a more controlled file I/O interface

#ifndef AMREX_NFILESTREAM_H_
#define AMREX_NFILESTREAM_H_

#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_VisMFBuffer.H>
#include <AMReX_GpuQualifiers.H>
#include <fstream>
#include <string>
#include <ios>

namespace amrex {

/**
 * \brief A custom file stream class for NFilesIter that replaces std::fstream
 * 
 * This class provides the same interface as std::fstream but allows for
 * custom implementations and better control over file I/O operations.
 */
class NFileStream
{
public:
    NFileStream() = default;
    ~NFileStream();

    // File operations
    void open(const std::string& filename, std::ios_base::openmode mode);
    void close();
    [[nodiscard]] bool is_open() const { return m_file.is_open(); }
    [[nodiscard]] bool good() const { return m_file.good(); }
    
    // Stream state
    [[nodiscard]] bool fail() const { return m_file.fail(); }
    [[nodiscard]] bool eof() const { return m_file.eof(); }
    void clear() { m_file.clear(); }
    
    // Write operations
    NFileStream& write(const char* s, std::streamsize n);
    NFileStream& flush();
    
    // Read operations
    NFileStream& read(char* s, std::streamsize n);
    
    // Positioning
    NFileStream& seekp(std::streampos pos);
    NFileStream& seekp(std::streamoff off, std::ios_base::seekdir way);
    [[nodiscard]] std::streampos tellp();
    
    // Buffer management
    std::streambuf* rdbuf() const { return m_file.rdbuf(); }
    
    // Conversion operators for compatibility
    explicit operator bool() const { return !fail(); }
    
    // Allow implicit conversion to std::istream& for reading operations
    operator std::istream&() { return m_file; }
    
    // Allow implicit conversion to std::ostream& for writing operations
    operator std::ostream&() { return m_file; }

private:
    std::fstream m_file;
};

} // namespace amrex

#endif // AMREX_NFILESTREAM_H_