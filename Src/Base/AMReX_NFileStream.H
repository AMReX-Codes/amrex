// ABOUTME: This file defines AMReX_NFileStream, a custom stream class for NFilesIter
// ABOUTME: It uses std::fstream with threaded writes and I/O hang detection

#ifndef AMREX_NFILESTREAM_H_
#define AMREX_NFILESTREAM_H_

#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_VisMFBuffer.H>
#include <AMReX_GpuQualifiers.H>
#include <fstream>
#include <string>
#include <ios>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <future>
#include <vector>
#include <memory>
#include <chrono>

namespace amrex {

/**
 * \brief A custom file stream class for NFilesIter that wraps std::fstream
 *
 * This class provides the same interface as std::fstream but uses threaded writes
 * with hang detection for better performance and reliability on network filesystems.
 */
class NFileStream
{
public:
    NFileStream() = default;
    ~NFileStream();

    // Disable copy and move operations (threading infrastructure cannot be safely moved)
    NFileStream(const NFileStream&) = delete;
    NFileStream& operator=(const NFileStream&) = delete;
    NFileStream(NFileStream&&) = delete;
    NFileStream& operator=(NFileStream&&) = delete;

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

    // Buffer management - for compatibility with VisMFBuffer
    std::streambuf* rdbuf() const { return m_file.rdbuf(); }

    // Conversion operators for compatibility
    explicit operator bool() const { return !fail(); }

    // Allow implicit conversion to std::istream& for reading operations
    operator std::istream&() { return m_file; }

    // Allow implicit conversion to std::ostream& for writing operations
    operator std::ostream&() { return m_file; }

private:
    std::fstream m_file;

    // Threading infrastructure for async writes
    struct WriteTask {
        std::vector<char> data;
        std::promise<bool> result;
        std::chrono::steady_clock::time_point start_time;
        WriteTask(const char* s, std::streamsize n)
            : data(s, s + n), start_time(std::chrono::steady_clock::now()) {}
    };

    std::thread m_write_thread;
    std::queue<std::unique_ptr<WriteTask>> m_write_queue;
    std::mutex m_queue_mutex;
    std::condition_variable m_queue_cv;
    std::atomic<bool> m_shutdown{false};
    std::atomic<bool> m_thread_running{false};

    // Hang detection configuration
    static constexpr std::chrono::seconds default_write_timeout{30};
    std::chrono::seconds m_write_timeout{default_write_timeout};

    // Helper methods
    void start_write_thread();
    void stop_write_thread();
    void write_thread_worker();
    bool perform_write_operation(const char* data, std::streamsize size);

public:
    // Hang detection configuration
    void set_write_timeout(std::chrono::seconds timeout) { m_write_timeout = timeout; }
    [[nodiscard]] std::chrono::seconds get_write_timeout() const { return m_write_timeout; }

    // Initialize timeout from runtime parameters
    void init_timeout_from_parmparse();

private:
};

} // namespace amrex

#endif // AMREX_NFILESTREAM_H_