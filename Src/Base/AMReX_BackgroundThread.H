#ifndef AMREX_BACKGROUND_THREAD_H_
#define AMREX_BACKGROUND_THREAD_H_

#include <algorithm>
#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>

namespace amrex {

class BackgroundThread
{
public:
    BackgroundThread ();
    ~BackgroundThread ();

    void Submit (std::function<void()>&& a_f);
    void Submit (std::function<void()> const& a_f);

    void Finish (); // Not required. Call this if you want all jobs to finish.

private:
    void do_job ();

    std::unique_ptr<std::thread> m_thread;
    std::mutex m_mutx;
    std::condition_variable m_job_cond;
    std::condition_variable m_done_cond;
    std::queue<std::function<void()> > m_func;
    bool m_clearing = false;
    bool m_finalizing = false;
};

}

#endif
