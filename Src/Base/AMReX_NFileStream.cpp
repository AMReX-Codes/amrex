// ABOUTME: This file implements AMReX_NFileStream, a custom stream class for NFilesIter
// ABOUTME: It provides controlled file I/O operations with POSIX error handling and exponential backoff

#include <AMReX_NFileStream.H>
#include <AMReX_Print.H>
#include <AMReX_BLassert.H>
#include <cerrno>
#include <chrono>
#include <thread>
#include <cstring>
#include <unistd.h>
#include <fcntl.h>

namespace amrex {

namespace {
    // Exponential backoff parameters
    constexpr int max_retries = 10;
    constexpr int initial_backoff_ms = 1;
    constexpr int max_backoff_ms = 1000;
    
    // Helper function to perform I/O with exponential backoff for retryable errors
    template <typename IOFunc>
    std::streamsize perform_io_with_retry(IOFunc&& io_func, const char* op_name)
    {
        int retry_count = 0;
        int backoff_ms = initial_backoff_ms;
        
        while (retry_count < max_retries) {
            errno = 0;
            std::streamsize result = io_func();
            
            if (result >= 0) {
                return result;  // Success
            }
            
            int err = errno;
            
            // Check for retryable errors
            if (err == EINTR || err == EAGAIN || err == EWOULDBLOCK) {
                if (retry_count > 0) {
                    // Exponential backoff
                    std::this_thread::sleep_for(std::chrono::milliseconds(backoff_ms));
                    backoff_ms = std::min(backoff_ms * 2, max_backoff_ms);
                }
                ++retry_count;
                continue;
            }
            
            // Non-retryable error
            std::string error_msg = std::string(op_name) + " failed: " + std::strerror(err);
            amrex::Error(error_msg);
        }
        
        // Max retries exceeded
        std::string error_msg = std::string(op_name) + " failed after " + 
                               std::to_string(max_retries) + " retries";
        amrex::Error(error_msg);
        return -1;  // Never reached due to Error() call
    }
}

NFileStream::~NFileStream()
{
    if (is_open()) {
        close();
    }
}

void NFileStream::open(const std::string& filename, std::ios_base::openmode mode)
{
    m_file.open(filename, mode);
    if (!m_file.good()) {
        int err = errno;
        std::string error_msg = "NFileStream::open failed to open file: " + filename;
        if (err != 0) {
            error_msg += " (errno: " + std::to_string(err) + " - " + std::strerror(err) + ")";
        }
        amrex::Error(error_msg);
    }
}

void NFileStream::close()
{
    if (is_open()) {
        // Flush before closing to ensure all data is written
        flush();
        
        // Close with retry logic
        int retry_count = 0;
        int backoff_ms = initial_backoff_ms;
        
        while (retry_count < max_retries) {
            errno = 0;
            m_file.close();
            
            if (!m_file.fail() || errno == 0) {
                return;  // Success
            }
            
            int err = errno;
            
            // Check for retryable errors
            if (err == EINTR) {
                if (retry_count > 0) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(backoff_ms));
                    backoff_ms = std::min(backoff_ms * 2, max_backoff_ms);
                }
                ++retry_count;
                m_file.clear();  // Clear error flags before retry
                continue;
            }
            
            // Non-retryable error - log but don't throw since we're in destructor path
            if (amrex::Verbose() > 0) {
                amrex::Print() << "Warning: NFileStream::close failed: " << std::strerror(err) << std::endl;
            }
            return;
        }
        
        if (amrex::Verbose() > 0) {
            amrex::Print() << "Warning: NFileStream::close failed after " << max_retries << " retries" << std::endl;
        }
    }
}

NFileStream& NFileStream::write(const char* s, std::streamsize n)
{
    if (n <= 0) return *this;
    
    std::streamsize total_written = 0;
    
    while (total_written < n) {
        std::streamsize to_write = n - total_written;
        
        std::streamsize written = perform_io_with_retry(
            [this, s, total_written, to_write]() -> std::streamsize {
                m_file.clear();  // Clear any error flags
                m_file.write(s + total_written, to_write);
                if (m_file.fail()) {
                    return -1;
                }
                return to_write;
            },
            "NFileStream::write"
        );
        
        total_written += written;
    }
    
    return *this;
}

NFileStream& NFileStream::flush()
{
    perform_io_with_retry(
        [this]() -> std::streamsize {
            m_file.clear();  // Clear any error flags
            m_file.flush();
            if (m_file.fail()) {
                return -1;
            }
            return 0;
        },
        "NFileStream::flush"
    );
    
    return *this;
}

NFileStream& NFileStream::read(char* s, std::streamsize n)
{
    if (n <= 0) return *this;
    
    std::streamsize total_read = 0;
    
    while (total_read < n) {
        std::streamsize to_read = n - total_read;
        
        std::streamsize bytes_read = perform_io_with_retry(
            [this, s, total_read, to_read]() -> std::streamsize {
                m_file.clear();  // Clear any error flags
                m_file.read(s + total_read, to_read);
                if (m_file.fail() && !m_file.eof()) {
                    return -1;
                }
                return m_file.gcount();
            },
            "NFileStream::read"
        );
        
        total_read += bytes_read;
        
        // Handle EOF
        if (m_file.eof() || bytes_read == 0) {
            break;
        }
    }
    
    return *this;
}

NFileStream& NFileStream::seekp(std::streampos pos)
{
    perform_io_with_retry(
        [this, pos]() -> std::streamsize {
            m_file.clear();  // Clear any error flags
            m_file.seekp(pos);
            if (m_file.fail()) {
                return -1;
            }
            return 0;
        },
        "NFileStream::seekp"
    );
    
    return *this;
}

NFileStream& NFileStream::seekp(std::streamoff off, std::ios_base::seekdir way)
{
    perform_io_with_retry(
        [this, off, way]() -> std::streamsize {
            m_file.clear();  // Clear any error flags
            m_file.seekp(off, way);
            if (m_file.fail()) {
                return -1;
            }
            return 0;
        },
        "NFileStream::seekp"
    );
    
    return *this;
}

std::streampos NFileStream::tellp()
{
    std::streampos pos = std::streampos(-1);
    
    perform_io_with_retry(
        [this, &pos]() -> std::streamsize {
            m_file.clear();  // Clear any error flags
            pos = m_file.tellp();
            if (pos == std::streampos(-1)) {
                return -1;
            }
            return 0;
        },
        "NFileStream::tellp"
    );
    
    return pos;
}

} // namespace amrex