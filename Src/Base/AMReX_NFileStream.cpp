// ABOUTME: This file implements AMReX_NFileStream using std::fstream with threaded writes
// ABOUTME: It provides controlled file I/O operations with hang detection and exponential backoff

#include <AMReX_NFileStream.H>
#include <AMReX_Print.H>
#include <AMReX_BLassert.H>
#include <AMReX_ParmParse.H>
#include <cerrno>
#include <chrono>
#include <thread>
#include <cstring>

namespace amrex {

namespace {
    // Exponential backoff parameters
    constexpr int max_retries = 10;
    constexpr int initial_backoff_ms = 1;
    constexpr int max_backoff_ms = 1000;

    // Helper function to perform I/O with exponential backoff for retryable errors
    template <typename IOFunc>
    bool perform_io_with_retry(IOFunc&& io_func, const char* op_name)
    {
        int retry_count = 0;
        int backoff_ms = initial_backoff_ms;

        while (retry_count < max_retries) {
            errno = 0;
            bool result = io_func();

            if (result) {
                return true;  // Success
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

            // Non-retryable error - log but don't throw during destruction
            if (amrex::Verbose() > 0) {
                std::string error_msg = std::string(op_name) + " failed";
                if (err != 0) {
                    error_msg += ": " + std::string(std::strerror(err));
                }
                amrex::Print() << "Error: " << error_msg << std::endl;
            }
            return false;
        }

        // Max retries exceeded - log but don't throw during destruction
        if (amrex::Verbose() > 0) {
            std::string error_msg = std::string(op_name) + " failed after " +
                                   std::to_string(max_retries) + " retries";
            amrex::Print() << "Error: " << error_msg << std::endl;
        }
        return false;
    }
}

NFileStream::~NFileStream()
{
    try {
        // Stop write thread first if it's running
        if (m_thread_running) {
            stop_write_thread();
        }
        
        if (is_open()) {
            close();
        }
    } catch (...) {
        // Suppress all exceptions in destructor to prevent segfaults
        // Log the error if verbose mode is enabled
        if (amrex::Verbose() > 0) {
            amrex::Print() << "Warning: Exception caught in NFileStream destructor" << std::endl;
        }
    }
}

void NFileStream::open(const std::string& filename, std::ios_base::openmode mode)
{
    if (is_open()) {
        close();
    }

    m_file.open(filename, mode);
    if (!m_file.good()) {
        int err = errno;
        std::string error_msg = "NFileStream::open failed to open file: " + filename;
        if (err != 0) {
            error_msg += " (errno: " + std::to_string(err) + " - " + std::strerror(err) + ")";
        }
        amrex::Error(error_msg);
    }

    // Initialize timeout from runtime parameters
    init_timeout_from_parmparse();
    
    // Start write thread if we're opening for writing
    if (mode & std::ios_base::out) {
        start_write_thread();
    }
}

void NFileStream::close()
{
    if (is_open()) {
        // Flush before closing if we were writing
        flush();

        // Close with retry logic
        perform_io_with_retry(
            [this]() -> bool {
                m_file.close();
                return !m_file.fail();
            },
            "NFileStream::close"
        );
    }
}

NFileStream& NFileStream::write(const char* s, std::streamsize n)
{
    if (!is_open() || n <= 0) {
        return *this;
    }

    // If write thread is running, queue the write operation
    if (m_thread_running) {
        auto task = std::make_unique<WriteTask>(s, n);
        auto future = task->result.get_future();

        {
            std::lock_guard<std::mutex> lock(m_queue_mutex);
            m_write_queue.push(std::move(task));
        }
        m_queue_cv.notify_one();

        // Wait for the write to complete with timeout
        auto future_status = future.wait_for(m_write_timeout);

        if (future_status == std::future_status::timeout) {
            // Write operation timed out - this indicates a hung filesystem
            if (amrex::Verbose() > 0) {
                amrex::Print() << "Warning: NFileStream write operation timed out after "
                              << m_write_timeout.count() << " seconds. Filesystem may be hung." << std::endl;
            }
            m_file.setstate(std::ios::failbit);
        } else if (!future.get()) {
            m_file.setstate(std::ios::failbit);
        }
    } else {
        // Fall back to synchronous write if no thread is running
        if (!perform_write_operation(s, n)) {
            m_file.setstate(std::ios::failbit);
        }
    }

    return *this;
}

NFileStream& NFileStream::flush()
{
    if (!is_open()) {
        return *this;
    }

    perform_io_with_retry(
        [this]() -> bool {
            m_file.flush();
            return !m_file.fail();
        },
        "NFileStream::flush"
    );

    return *this;
}

NFileStream& NFileStream::read(char* s, std::streamsize n)
{
    if (!is_open() || n <= 0) {
        return *this;
    }

    perform_io_with_retry(
        [this, s, n]() -> bool {
            m_file.read(s, n);
            return !m_file.fail() || m_file.eof();
        },
        "NFileStream::read"
    );

    return *this;
}

NFileStream& NFileStream::seekp(std::streampos pos)
{
    if (!is_open()) {
        return *this;
    }

    perform_io_with_retry(
        [this, pos]() -> bool {
            m_file.seekp(pos);
            return !m_file.fail();
        },
        "NFileStream::seekp"
    );

    return *this;
}

NFileStream& NFileStream::seekp(std::streamoff off, std::ios_base::seekdir way)
{
    if (!is_open()) {
        return *this;
    }

    perform_io_with_retry(
        [this, off, way]() -> bool {
            m_file.seekp(off, way);
            return !m_file.fail();
        },
        "NFileStream::seekp"
    );

    return *this;
}

std::streampos NFileStream::tellp()
{
    if (!is_open()) {
        return std::streampos(-1);
    }

    std::streampos pos = std::streampos(-1);

    perform_io_with_retry(
        [this, &pos]() -> bool {
            pos = m_file.tellp();
            return pos != std::streampos(-1);
        },
        "NFileStream::tellp"
    );

    return pos;
}

void NFileStream::start_write_thread()
{
    if (!m_thread_running) {
        m_shutdown = false;
        m_thread_running = true;
        m_write_thread = std::thread(&NFileStream::write_thread_worker, this);
    }
}

void NFileStream::stop_write_thread()
{
    if (m_thread_running) {
        try {
            // Signal shutdown
            m_shutdown = true;
            m_queue_cv.notify_all();

            // Wait for thread to finish with timeout to prevent hanging
            if (m_write_thread.joinable()) {
                // Use detach() if join() might hang during destruction
                auto future = std::async(std::launch::async, [this]() {
                    m_write_thread.join();
                });
                
                auto status = future.wait_for(std::chrono::seconds(5));
                if (status == std::future_status::timeout) {
                    // Thread didn't finish in time, detach it to prevent hanging
                    if (amrex::Verbose() > 0) {
                        amrex::Print() << "Warning: Write thread did not finish cleanly, detaching" << std::endl;
                    }
                    m_write_thread.detach();
                }
            }

            m_thread_running = false;

            // Clear any remaining tasks
            std::lock_guard<std::mutex> lock(m_queue_mutex);
            while (!m_write_queue.empty()) {
                auto task = std::move(m_write_queue.front());
                m_write_queue.pop();
                try {
                    task->result.set_value(false);  // Mark as failed
                } catch (...) {
                    // Ignore promise exceptions during cleanup
                }
            }
        } catch (...) {
            // Ensure thread state is reset even if cleanup fails
            m_thread_running = false;
            if (m_write_thread.joinable()) {
                m_write_thread.detach();
            }
        }
    }
}

void NFileStream::write_thread_worker()
{
    while (!m_shutdown) {
        std::unique_ptr<WriteTask> task;

        // Wait for a task or shutdown signal
        {
            std::unique_lock<std::mutex> lock(m_queue_mutex);
            m_queue_cv.wait(lock, [this] { return !m_write_queue.empty() || m_shutdown; });

            if (m_shutdown && m_write_queue.empty()) {
                break;
            }

            if (!m_write_queue.empty()) {
                task = std::move(m_write_queue.front());
                m_write_queue.pop();
            }
        }

        // Perform the write operation if we have a task
        if (task) {
            bool success = perform_write_operation(task->data.data(), task->data.size());
            task->result.set_value(success);
        }
    }
}

bool NFileStream::perform_write_operation(const char* data, std::streamsize size)
{
    return perform_io_with_retry(
        [this, data, size]() -> bool {
            m_file.write(data, size);
            return !m_file.fail();
        },
        "NFileStream::perform_write_operation"
    );
}

void NFileStream::init_timeout_from_parmparse()
{
    amrex::ParmParse pp("nfilestream");
    int timeout_seconds = static_cast<int>(default_write_timeout.count());
    pp.queryWithParser("write_timeout", timeout_seconds);

    if (timeout_seconds > 0) {
        m_write_timeout = std::chrono::seconds(timeout_seconds);
    } else {
        // Use default timeout if invalid value provided
        m_write_timeout = default_write_timeout;
    }
}

} // namespace amrex