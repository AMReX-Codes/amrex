// ABOUTME: This file implements AMReX_NFileStream using std::fstream with threaded writes
// ABOUTME: It provides a simple wrapper around std::fstream with I/O hang detection

#include <AMReX_NFileStream.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

namespace amrex {

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
        std::string error_msg = "NFileStream::open failed to open file: " + filename;
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
        m_file.close();
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
            std::string error_msg = "NFileStream write operation timed out after " +
                                   std::to_string(m_write_timeout.count()) +
                                   " seconds. Filesystem may be hung.";
            amrex::Abort(error_msg);
        } else if (!future.get()) {
            m_file.setstate(std::ios::failbit);
        }
    } else {
        // Fall back to synchronous write if no thread is running
        m_file.write(s, n);
    }

    return *this;
}

NFileStream& NFileStream::flush()
{
    if (is_open()) {
        m_file.flush();
    }
    return *this;
}

NFileStream& NFileStream::read(char* s, std::streamsize n)
{
    if (is_open() && n > 0) {
        m_file.read(s, n);
    }
    return *this;
}

NFileStream& NFileStream::seekp(std::streampos pos)
{
    if (is_open()) {
        m_file.seekp(pos);
    }
    return *this;
}

NFileStream& NFileStream::seekp(std::streamoff off, std::ios_base::seekdir way)
{
    if (is_open()) {
        m_file.seekp(off, way);
    }
    return *this;
}

std::streampos NFileStream::tellp()
{
    if (is_open()) {
        return m_file.tellp();
    }
    return std::streampos(-1);
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
            m_file.write(task->data.data(), task->data.size());
            bool success = !m_file.fail();
            task->result.set_value(success);
        }
    }
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