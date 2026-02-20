#ifndef BL_LAZY_H
#define BL_LAZY_H
#include <AMReX_Config.H>

/**
 * \file AMReX_Lazy.H
 * \brief Deferred-evaluation reduction queue for lazy AMReX operations.
 */

#include <vector>
#include <functional>
#include <algorithm>

namespace amrex {
/**
 * \ingroup amrex_utilities
 * \brief Namespace for deferred (lazy) reduction evaluation.
 *
 * Reduction operations can be queued with \c QueueReduction and executed
 * later in bulk with \c EvalReduction.
 */
namespace Lazy
{
    using Func = std::function<void()>;
    using FuncQue = std::vector<Func>;

    extern FuncQue reduction_queue;

    //! Enqueue a reduction function \p f to be evaluated later.
    void QueueReduction (Func);
    //! Execute all queued reduction functions and clear the queue.
    void EvalReduction ();

    //! Finalize the Lazy subsystem (called at AMReX finalization).
    void Finalize ();
}
}

#endif
