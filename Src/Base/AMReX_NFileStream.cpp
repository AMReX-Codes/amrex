// ABOUTME: This file implements AMReX_NFileStream using POSIX I/O directly
// ABOUTME: It provides controlled file I/O operations with proper error handling and exponential backoff

#include <AMReX_NFileStream.H>
#include <AMReX_Print.H>
#include <AMReX_BLassert.H>
#include <AMReX_ParmParse.H>
#include <cerrno>
#include <chrono>
#include <thread>
#include <cstring>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

namespace amrex {

namespace {
    // Exponential backoff parameters
    constexpr int max_retries = 10;
    constexpr int initial_backoff_ms = 1;
    constexpr int max_backoff_ms = 1000;

    // Helper function to perform I/O with exponential backoff for retryable errors
    template <typename IOFunc>
    auto perform_io_with_retry(IOFunc&& io_func, const char* op_name) -> decltype(io_func())
    {
        int retry_count = 0;
        int backoff_ms = initial_backoff_ms;

        while (retry_count < max_retries) {
            errno = 0;
            auto result = io_func();

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
        return decltype(io_func())(-1);  // Never reached due to Error() call
    }
}

NFileStream::~NFileStream()
{
    // Stop write thread first if it's running
    if (m_thread_running) {
        stop_write_thread();
    }
    
    if (is_open()) {
        close();
    }
}


int NFileStream::convert_openmode_to_flags(std::ios_base::openmode mode)
{
    int flags = 0;

    bool read_mode = (mode & std::ios_base::in) != 0;
    bool write_mode = (mode & std::ios_base::out) != 0;
    bool append_mode = (mode & std::ios_base::app) != 0;
    bool trunc_mode = (mode & std::ios_base::trunc) != 0;

    if (read_mode && write_mode) {
        flags = O_RDWR;
    } else if (write_mode) {
        flags = O_WRONLY;
    } else {
        flags = O_RDONLY;
    }

    if (write_mode) {
        flags |= O_CREAT;

        if (append_mode) {
            flags |= O_APPEND;
        } else if (trunc_mode) {
            flags |= O_TRUNC;
        }
    }

    if (mode & std::ios_base::binary) {
        // Binary mode is default in POSIX, no special flag needed
    }

    return flags;
}

void NFileStream::set_error_state(bool fail_state, bool eof_state)
{
    m_fail = fail_state;
    m_eof = eof_state;
    m_good = !fail_state && !eof_state;
}

void NFileStream::open(const std::string& filename, std::ios_base::openmode mode)
{
    if (is_open()) {
        close();
    }

    m_mode = mode;
    int flags = convert_openmode_to_flags(mode);
    mode_t file_mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;  // 644 permissions

    m_fd = perform_io_with_retry(
        [&filename, flags, file_mode]() -> int {
            return ::open(filename.c_str(), flags, file_mode);
        },
        "NFileStream::open"
    );

    if (m_fd >= 0) {
        clear();  // Reset error flags on successful open
        // Initialize timeout from runtime parameters
        init_timeout_from_parmparse();
        // Start write thread if we're opening for writing
        if (mode & std::ios_base::out) {
            start_write_thread();
        }
    } else {
        set_error_state(true);
    }
}

void NFileStream::close()
{
    if (is_open()) {
        // Flush before closing if we were writing
        if (m_mode & std::ios_base::out) {
            flush();
        }

        // Close with retry logic
        int retry_count = 0;
        int backoff_ms = initial_backoff_ms;

        while (retry_count < max_retries) {
            errno = 0;
            int result = ::close(m_fd);

            if (result == 0) {
                m_fd = -1;
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
                continue;
            }

            // Non-retryable error - log but don't throw since we're in destructor path
            if (amrex::Verbose() > 0) {
                amrex::Print() << "Warning: NFileStream::close failed: " << std::strerror(err) << std::endl;
            }
            m_fd = -1;  // Mark as closed even if close() failed
            return;
        }

        if (amrex::Verbose() > 0) {
            amrex::Print() << "Warning: NFileStream::close failed after " << max_retries << " retries" << std::endl;
        }
        m_fd = -1;  // Mark as closed even if close() failed
    }
}

NFileStream& NFileStream::write(const char* s, std::streamsize n)
{
    if (!is_open() || n <= 0) {
        set_error_state(true);
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
            set_error_state(true);
        } else if (!future.get()) {
            set_error_state(true);
        }
    } else {
        // Fall back to synchronous write if no thread is running
        if (!perform_write_syscall(s, n)) {
            set_error_state(true);
        }
    }

    return *this;
}

NFileStream& NFileStream::flush()
{
    if (!is_open()) {
        set_error_state(true);
        return *this;
    }

    // For POSIX, we use fsync to ensure data is written to storage
    int result = perform_io_with_retry(
        [this]() -> int {
            return ::fsync(m_fd);
        },
        "NFileStream::flush"
    );

    if (result != 0) {
        set_error_state(true);
    }

    return *this;
}

NFileStream& NFileStream::read(char* s, std::streamsize n)
{
    if (!is_open() || n <= 0) {
        set_error_state(true);
        return *this;
    }

    std::streamsize total_read = 0;

    while (total_read < n) {
        std::streamsize to_read = n - total_read;

        ssize_t bytes_read = perform_io_with_retry(
            [this, s, total_read, to_read]() -> ssize_t {
                return ::read(m_fd, s + total_read, static_cast<size_t>(to_read));
            },
            "NFileStream::read"
        );

        if (bytes_read > 0) {
            total_read += bytes_read;
        } else if (bytes_read == 0) {
            // EOF reached
            set_error_state(false, true);
            break;
        } else {
            set_error_state(true);
            break;
        }
    }

    return *this;
}

NFileStream& NFileStream::seekp(std::streampos pos)
{
    if (!is_open()) {
        set_error_state(true);
        return *this;
    }

    off_t result = perform_io_with_retry(
        [this, pos]() -> off_t {
            return ::lseek(m_fd, static_cast<off_t>(pos), SEEK_SET);
        },
        "NFileStream::seekp"
    );

    if (result == static_cast<off_t>(-1)) {
        set_error_state(true);
    }

    return *this;
}

NFileStream& NFileStream::seekp(std::streamoff off, std::ios_base::seekdir way)
{
    if (!is_open()) {
        set_error_state(true);
        return *this;
    }

    int whence;
    switch (way) {
        case std::ios_base::beg:
            whence = SEEK_SET;
            break;
        case std::ios_base::cur:
            whence = SEEK_CUR;
            break;
        case std::ios_base::end:
            whence = SEEK_END;
            break;
        default:
            set_error_state(true);
            return *this;
    }

    off_t result = perform_io_with_retry(
        [this, off, whence]() -> off_t {
            return ::lseek(m_fd, static_cast<off_t>(off), whence);
        },
        "NFileStream::seekp"
    );

    if (result == static_cast<off_t>(-1)) {
        set_error_state(true);
    }

    return *this;
}

std::streampos NFileStream::tellp()
{
    if (!is_open()) {
        set_error_state(true);
        return std::streampos(-1);
    }

    off_t pos = perform_io_with_retry(
        [this]() -> off_t {
            return ::lseek(m_fd, 0, SEEK_CUR);
        },
        "NFileStream::tellp"
    );

    if (pos == static_cast<off_t>(-1)) {
        set_error_state(true);
        return std::streampos(-1);
    }

    return std::streampos(pos);
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
        // Signal shutdown
        m_shutdown = true;
        m_queue_cv.notify_all();

        // Wait for thread to finish
        if (m_write_thread.joinable()) {
            m_write_thread.join();
        }

        m_thread_running = false;

        // Clear any remaining tasks
        std::lock_guard<std::mutex> lock(m_queue_mutex);
        while (!m_write_queue.empty()) {
            auto task = std::move(m_write_queue.front());
            m_write_queue.pop();
            task->result.set_value(false);  // Mark as failed
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
            bool success = perform_write_syscall(task->data.data(), task->data.size());
            task->result.set_value(success);
        }
    }
}

bool NFileStream::perform_write_syscall(const char* data, std::streamsize size)
{
    std::streamsize total_written = 0;

    while (total_written < size) {
        std::streamsize to_write = size - total_written;

        ssize_t written = perform_io_with_retry(
            [this, data, total_written, to_write]() -> ssize_t {
                return ::write(m_fd, data + total_written, static_cast<size_t>(to_write));
            },
            "NFileStream::perform_write_syscall"
        );

        if (written > 0) {
            total_written += written;
        } else {
            return false;
        }
    }

    return true;
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
