#ifndef AMREX_BACKGROUND_THREAD_H_
#define AMREX_BACKGROUND_THREAD_H_
#include <AMReX_Config.H>

#include <algorithm>
#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>

/**
 * \file AMReX_BackgroundThread.H
 * \brief Simple job queue serviced by a dedicated background std::thread.
 */

namespace amrex {

/**
 * \class amrex::BackgroundThread
 * \brief Utility class running deferred work on a single background thread.
 * \ingroup amrex_utilities
 */
class BackgroundThread
{
public:
    /**
     * \brief Launch the worker thread and initialize synchronization primitives.
     */
    BackgroundThread ();
    /**
     * \brief Stop the worker thread and drain any queued work.
     */
    ~BackgroundThread ();
    BackgroundThread (BackgroundThread const&) = delete;
    BackgroundThread (BackgroundThread &&) = delete;
    BackgroundThread& operator= (BackgroundThread const&) = delete;
    BackgroundThread& operator= (BackgroundThread &&) = delete;

    /**
     * \brief Enqueue a job that will run asynchronously on the worker thread.
     * \param a_f Functor representing the work to execute.
     */
    void Submit (std::function<void()>&& a_f);
    /**
     * \overload
     */
    void Submit (std::function<void()> const& a_f);

    /**
     * \brief Block the caller until all queued jobs have completed.
     */
    void Finish ();

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
