#include <AMReX_BackgroundThread.H>

#ifndef _WIN32
#include <csignal>
#include <pthread.h>
#endif

namespace amrex {

BackgroundThread::BackgroundThread ()
{
    m_thread = std::make_unique<std::thread>(&BackgroundThread::do_job, this);
}

BackgroundThread::~BackgroundThread ()
{
    if (m_thread) {
        Submit([this] () { m_finalizing = true; });
        m_thread->join();
        m_thread.reset();
    }
}

std::thread::id BackgroundThread::ThreadId () const
{
    return m_thread ? m_thread->get_id() : std::thread::id{};
}

void BackgroundThread::do_job ()
{
#ifndef _WIN32
    // SIGINT and SIGTERM go to other threads, so that the signal handler of
    // AMReX can let the queued jobs finish here.
    sigset_t set;
    sigemptyset(&set);
    sigaddset(&set, SIGINT);
    sigaddset(&set, SIGTERM);
    pthread_sigmask(SIG_BLOCK, &set, nullptr);
#endif

    while (true)
    {
        std::unique_lock<std::mutex> lck(m_mutx);
        m_job_cond.wait(lck, [this] () -> bool { return !m_func.empty(); });
        auto f = m_func.front();
        m_func.pop();
        lck.unlock();
        f();
        if (m_clearing) { // All jobs before this have finished.
            m_done_cond.notify_one();
        }
        if (m_finalizing) {
            break;
        }
    }
}

void BackgroundThread::Submit (std::function<void()>&& a_f)
{
    std::scoped_lock lck(m_mutx);
    m_func.emplace(std::move(a_f));
    m_job_cond.notify_one();
}

void BackgroundThread::Submit (std::function<void()> const& a_f)
{
    std::scoped_lock lck(m_mutx);
    m_func.emplace(a_f);
    m_job_cond.notify_one();
}

void BackgroundThread::Finish ()
{
    if (m_thread) {
        Submit([this] () { m_clearing = true; });
        std::unique_lock<std::mutex> lck(m_mutx);
        m_done_cond.wait(lck, [this] () -> bool { return m_func.empty(); });
        m_clearing = false;
        lck.unlock();
    }
}

}
