#include <AMReX_BackgroundThread.H>

namespace amrex {

BackgroundThread::BackgroundThread ()
{
    m_thread.reset(new std::thread(&BackgroundThread::do_job, this));
}

BackgroundThread::~BackgroundThread ()
{
    if (m_thread) {
        Submit([this] () { m_finalizing = true; });
        m_thread->join();
        m_thread.reset();
    }
}

void BackgroundThread::do_job ()
{
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
    std::lock_guard<std::mutex> lck(m_mutx);
    m_func.emplace(std::move(a_f));
    m_job_cond.notify_one();
}

void BackgroundThread::Submit (std::function<void()> const& a_f)
{
    std::lock_guard<std::mutex> lck(m_mutx);
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
